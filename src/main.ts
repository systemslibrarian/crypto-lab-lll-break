import './style.css';
import { bkzReduce, buildLWELattice, generateLWE, recoverFromReducedBasis } from './bkz';
import { determinant, lllReduce, orthogonalityDefect, shortestVectorInBasis } from './lll';
import { dot, gramSchmidt, norm, type Basis, type Vec } from './lattice';
import { getSeed, randomIntInclusive, setSeed } from './rng';

const app = document.querySelector<HTMLDivElement>('#app');
if (!app) {
  throw new Error('Missing app root');
}

const theme = document.documentElement.getAttribute('data-theme') ?? 'dark';

app.innerHTML = `
<a class="skip-link" href="#main-content">Skip to main content</a>
<div class="topbar">
  <h1>crypto-lab-lll-break</h1>
  <button id="theme-toggle" type="button" class="theme-toggle" aria-label="${theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode'}">${theme === 'dark' ? '🌙' : '☀️'}</button>
</div>

<main id="main-content">
<section class="exhibit" id="lab-controls" aria-labelledby="lab-controls-title">
  <h2 id="lab-controls-title">Reproducible Labs</h2>
  <p class="objective"><strong>Goal:</strong> make any run repeatable. Enter a numeric seed to make all randomness deterministic - the same seed and parameters always produce the same instance - or leave it blank for fresh crypto-random exploration. "Copy lab link" puts the whole setup in the URL so every student opens the identical lab.</p>
  <div class="controls-row">
    <label>Seed (blank = random)
      <input id="lab-seed" type="text" inputmode="numeric" autocomplete="off" placeholder="e.g. 1234" />
    </label>
    <label>Canonical lab
      <select id="lab-preset">
        <option value="">-- choose a lab --</option>
        <option value="swap-cascade">LLL swap cascade</option>
        <option value="delta-sensitivity">Delta sensitivity</option>
        <option value="toy-lwe-success">Toy LWE success</option>
        <option value="bkz-helps">BKZ helps</option>
        <option value="kyber-failure">Kyber-like failure</option>
      </select>
    </label>
    <button id="lab-apply-seed" type="button" class="btn">Apply seed</button>
    <button id="lab-copy-link" type="button" class="btn">Copy lab link</button>
  </div>
  <p class="mono" id="lab-status" role="status" aria-live="polite">Mode: crypto-random (not reproducible).</p>
</section>
<section class="exhibit" id="exhibit-1" aria-labelledby="exhibit-1-title">
  <h2 id="exhibit-1-title">Exhibit 1 - What Is a Lattice?</h2>
  <p class="objective"><strong>Goal:</strong> see that many different bases generate the same lattice, and that the determinant (fundamental-cell area) is invariant under unimodular basis changes.</p>
  <div class="grid two">
    <div>
      <canvas id="lattice-canvas" width="500" height="400" role="img" aria-label="2D lattice canvas"></canvas>
      <p class="mono" id="det-label"></p>
      <p class="mono" id="same-lattice-label"></p>
    </div>
    <div>
      <label>b1 angle <input id="b1-angle" type="range" min="0" max="180" value="25" /></label>
      <label>b1 length <input id="b1-len" type="range" min="1" max="10" value="5" /></label>
      <label>b2 angle <input id="b2-angle" type="range" min="0" max="180" value="95" /></label>
      <label>b2 length <input id="b2-len" type="range" min="1" max="10" value="4" /></label>
      <button id="same-lattice-btn" type="button" class="btn">Same lattice, different basis</button>
      <p>
        A lattice is all integer combinations of two basis vectors. Determinant equals area of the
        fundamental parallelogram, an invariant under unimodular basis changes.
      </p>
    </div>
  </div>
  <details class="selfcheck">
    <summary>Can you answer this now?</summary>
    <p>Q: You click "Same lattice, different basis" and the determinant label does not change. Why must that be the case? A: The new basis is the old one multiplied by an integer matrix with determinant +-1 (unimodular). Multiplying basis vectors by such a matrix re-labels the same set of lattice points, so the fundamental-cell area - the absolute determinant - is preserved even though the vectors look different.</p>
  </details>
</section>

<section class="exhibit" id="exhibit-2" aria-labelledby="exhibit-2-title">
  <h2 id="exhibit-2-title">Exhibit 2 - Gram-Schmidt Orthogonalization</h2>
  <p class="objective"><strong>Goal:</strong> compute the Gram-Schmidt vectors b*, read the mu coefficients, and use them to decide the Lovasz condition.</p>
  <div class="grid two">
    <div>
      <div class="matrix-grid">
        <label class="sr-only" for="gso-a11">GSO matrix row 1 column 1</label>
        <input id="gso-a11" type="number" inputmode="numeric" value="3" />
        <label class="sr-only" for="gso-a12">GSO matrix row 1 column 2</label>
        <input id="gso-a12" type="number" inputmode="numeric" value="1" />
        <label class="sr-only" for="gso-a21">GSO matrix row 2 column 1</label>
        <input id="gso-a21" type="number" inputmode="numeric" value="1" />
        <label class="sr-only" for="gso-a22">GSO matrix row 2 column 2</label>
        <input id="gso-a22" type="number" inputmode="numeric" value="2" />
      </div>
      <button id="gso-next" type="button" class="btn">Next</button>
      <pre id="gso-text" class="mono panel"></pre>
      <p class="panel">
        Gram-Schmidt vectors and mu coefficients are floating point. Basis vectors remain exact integers.
      </p>
    </div>
    <div>
      <canvas id="gso-canvas" width="500" height="360" role="img" aria-label="Gram-Schmidt visualization"></canvas>
    </div>
  </div>
  <details class="selfcheck">
    <summary>Can you answer this now?</summary>
    <p>Q: mu21 measures what, and why does LLL round it to the nearest integer during size reduction? A: mu21 = &lt;b2, b1*&gt; / &lt;b1*, b1*&gt; is how much of b1* sits inside b2. Subtracting round(mu21)*b1 from b2 removes that overlap using an integer multiple, which keeps b2 a lattice vector while shrinking it - you cannot subtract the exact real multiple without leaving the lattice.</p>
  </details>
</section>

<section class="exhibit" id="exhibit-3" aria-labelledby="exhibit-3-title">
  <h2 id="exhibit-3-title">Exhibit 3 - LLL Step-by-Step</h2>
  <p class="objective"><strong>Goal:</strong> predict whether LLL will size-reduce, swap, or advance at each step - and watch the metrics confirm the lattice never changes (det fixed, transform U unimodular).</p>
  <div class="controls-row">
    <label>Dimension
      <select id="lll-dim">
        <option value="2">2D</option>
        <option value="3">3D</option>
      </select>
    </label>
    <label>Delta <input id="lll-delta" type="range" min="0.5" max="0.999" step="0.001" value="0.75" /></label>
    <button type="button" class="btn preset" data-preset="random2d">Random bad 2D</button>
    <button type="button" class="btn preset" data-preset="classic">Classic example</button>
    <button type="button" class="btn preset" data-preset="near">Near-orthogonal</button>
    <button type="button" class="btn preset" data-preset="challenge3d">3D challenge</button>
    <button type="button" class="btn preset" data-preset="qary">LWE-style q-ary</button>
  </div>
  <label class="sr-only" for="lll-matrix">LLL input basis matrix, one row of space-separated integers per line</label>
  <textarea id="lll-matrix" class="mono" rows="4" aria-label="LLL input basis matrix, one row of space-separated integers per line">19 2
7 1</textarea>
  <div class="controls-row">
    <button id="lll-step" type="button" class="btn">Step</button>
    <button id="lll-auto" type="button" class="btn">Auto</button>
    <button id="lll-reset" type="button" class="btn">Reset</button>
  </div>
  <div class="grid two">
    <div>
      <canvas id="lll-canvas" width="420" height="420" role="img" aria-label="LLL basis evolution canvas"></canvas>
      <canvas id="defect-chart" width="420" height="140" role="img" aria-label="Orthogonality defect chart"></canvas>
    </div>
    <div>
      <div id="lll-log" class="panel mono log" aria-live="polite"></div>
      <pre id="lll-metrics" class="panel mono"></pre>
    </div>
  </div>
  <details class="selfcheck">
    <summary>Can you answer this now?</summary>
    <p>Q: After many size-reductions and swaps the basis vectors are completely different numbers, yet "Lattice det |det B|" never changes. What does that prove? A: Every LLL step is a unimodular row operation, so the whole reduction is reducedBasis = U x original with det(U) = +-1 (shown live as the matrix U). The basis is rewritten but the lattice - the set of points and its determinant - is identical. LLL improves basis quality, it does not change the lattice.</p>
  </details>
</section>

<section class="exhibit" id="exhibit-4" aria-labelledby="exhibit-4-title">
  <h2 id="exhibit-4-title">Exhibit 4 - Break a Toy LWE Instance</h2>
  <p class="objective"><strong>Goal:</strong> derive why (-s, e, &plusmn;1) is short in the primal embedding, then recover s only when reduction actually exposes that vector - never confusing it with the brute-force baseline.</p>
  <p>
    The instance b = A s + e (mod q) is placed in a primal (Kannan) embedding whose rows are
    [I | A&#7488; | 0], [0 | qI | 0], [0 | b | 1]. The vector (-s, e, &plusmn;1) is short, so LLL/BKZ
    surfaces it and the secret is read straight off the first block. Recovery here is genuine
    lattice reduction; the exhaustive-search baseline only appears, clearly labelled, when reduction
    fails.
  </p>
  <details class="model-note">
    <summary>What this toy BKZ models / omits</summary>
    <p>
      Models: LLL preprocessing, local block enumeration (exact SVP on small blocks),
      lattice-preserving insert-and-reduce, tours, and convergence logging - the shape of BKZ
      reasoning. Omits: projected-block handling, pruning, floating-point/precision engineering,
      progressive and modern BKZ variants, and real cost models. Consequence: it teaches how BKZ
      thinks, not production cryptanalysis - use a lattice estimator and a real reduction library for
      security claims.
    </p>
  </details>
  <div class="controls-row">
    <label>n <input id="lwe-n" type="range" min="2" max="12" value="4" /></label>
    <label>q <input id="lwe-q" type="range" min="7" max="257" value="71" /></label>
    <label>sigma <input id="lwe-sigma" type="range" min="1" max="10" value="2" /></label>
    <label>beta <input id="lwe-beta" type="range" min="2" max="8" value="2" /></label>
    <button id="lwe-generate" type="button" class="btn">Generate LWE Instance</button>
    <button id="lwe-attack" type="button" class="btn">Run LLL/BKZ Attack</button>
    <button id="kyber-try" type="button" class="btn warn">Try Kyber-512 parameters</button>
  </div>
  <progress id="lwe-progress" max="100" value="0" aria-label="LWE attack progress"></progress>
  <div class="lwe-meter panel" role="status" aria-live="polite">
    <div class="lwe-meter-head mono" id="lwe-meter-text">Norm-gap confidence: waiting for attack run.</div>
    <div class="lwe-meter-track" aria-hidden="true">
      <div id="lwe-meter-fill" class="lwe-meter-fill"></div>
    </div>
  </div>
  <pre id="lwe-output" class="panel mono" aria-live="polite"></pre>
  <details class="selfcheck">
    <summary>Can you answer this now?</summary>
    <p>Q: The output says "Lattice attack result: FAILED" but "Teaching baseline: brute force found the secret". Did you break LWE? A: No. Brute force just enumerated every candidate secret - it is exhaustive search, not lattice reduction, and it only works because the toy secret space is tiny. A genuine lattice break requires reduction to expose the (-s, e, +-1) vector. Counting the baseline as a break is exactly the mistake this exhibit guards against.</p>
  </details>
</section>

<section class="exhibit" id="exhibit-5" aria-labelledby="exhibit-5-title">
  <h2 id="exhibit-5-title">Exhibit 5 - Parameter Explorer</h2>
  <p class="objective"><strong>Goal:</strong> explain why tiny parameters fall to LLL/BKZ while Kyber-like dimensions push the required block size - and the cost - out of reach.</p>
  <label>n <input id="explore-n" type="range" min="4" max="256" value="8" /></label>
  <div class="grid two">
    <pre id="explore-left" class="panel mono"></pre>
    <pre id="explore-right" class="panel mono"></pre>
  </div>
  <canvas id="threshold-chart" width="900" height="260" role="img" aria-label="Security threshold chart"></canvas>
  <p>
    ML-KEM (Kyber), FrodoKEM, ML-DSA (Dilithium), and FALCON all rely on LLL/BKZ failing at chosen parameters.
  </p>
  <details class="selfcheck">
    <summary>Can you answer this now?</summary>
    <p>Q: Moving n upward, the required block size beta and the cost exponent climb steeply. Which single fact makes real lattice schemes secure? A: BKZ cost grows roughly as 2^(0.292*beta), and the beta needed to find the short vector grows with the dimension. At Kyber-like parameters the required beta puts the cost far beyond practical resources - the attack shape still exists, it is just computationally out of reach. (This is intuition, not a certified security estimate.)</p>
  </details>
</section>

<section class="exhibit" id="challenges" aria-labelledby="challenges-title">
  <h2 id="challenges-title">Challenges</h2>
  <p class="objective"><strong>Goal:</strong> predict first, then check. Each challenge loads a reproducible (seeded) setup - commit to an answer out loud before you reveal it.</p>

  <div class="challenge">
    <p><strong>1. Will it swap?</strong> Load the setup. Before pressing Step, predict: does LLL's first Lovasz check at k=1 pass (advance) or fail (swap)?</p>
    <button type="button" class="btn challenge-setup" data-lab="swap-cascade">Set up this challenge</button>
    <details class="selfcheck"><summary>Reveal answer</summary><p>It swaps - repeatedly. With b1 and b2 nearly parallel and almost equal length, after size reduction b2* is much shorter than b1*, so delta*||b1*||^2 &gt; ||b2* + mu*b1*||^2 and the Lovasz condition fails. Each swap then exposes a new violation, producing a cascade until the basis is reduced.</p></details>
  </div>

  <div class="challenge">
    <p><strong>2. Delta sensitivity.</strong> The setup starts at delta=0.5. Predict whether raising delta toward 0.99 (then Reset) produces more or fewer swaps, and why.</p>
    <button type="button" class="btn challenge-setup" data-lab="delta-sensitivity">Set up this challenge</button>
    <details class="selfcheck"><summary>Reveal answer</summary><p>Higher delta is a stricter Lovasz bound, so it demands a better-ordered basis and generally triggers MORE swaps and a shorter final first vector. Lower delta (down to 0.5) accepts weaker bases, terminating faster with fewer swaps but a worse basis. Delta trades reduction quality against work.</p></details>
  </div>

  <div class="challenge">
    <p><strong>3. Same lattice, ugly basis.</strong> In Exhibit 1, find a basis that generates the same lattice as the current one but looks much worse (longer, more skewed vectors). What operation guarantees you stay on the same lattice?</p>
    <button type="button" class="btn challenge-setup" data-lab="swap-cascade">Jump to Exhibit 1 idea</button>
    <details class="selfcheck"><summary>Reveal answer</summary><p>Apply any unimodular transform - e.g. b2 &lt;- b2 + k*b1 for an integer k (the "Same lattice, different basis" button does exactly this). Adding an integer multiple of one basis vector to another, or swapping them, has determinant +-1, so the lattice (and its determinant) is unchanged while the vectors get longer and more skewed.</p></details>
  </div>

  <div class="challenge">
    <p><strong>4. Read the embedding vector.</strong> Run the toy LWE setup and attack. In the printed short vector, identify which coordinates are the secret block, which are the error block, and what the final coordinate means.</p>
    <button type="button" class="btn challenge-setup" data-lab="toy-lwe-success">Set up this challenge</button>
    <details class="selfcheck"><summary>Reveal answer</summary><p>The vector is (-s, e, +-1) over [n secret coords | m error coords | 1 embed coord]. The first n entries are the secret (negated by the sign of the last coordinate); the next m are the LWE errors; the final +-1 is the embedding coordinate that anchors the shift. The app prints exactly this decomposition under "Structure (-s, e, +-1)".</p></details>
  </div>

  <div class="challenge">
    <p><strong>5. Break it until it breaks.</strong> Starting from the toy success setup, raise sigma and/or n until the lattice attack fails (not the brute-force baseline). Which parameter degrades the attack fastest, and why?</p>
    <button type="button" class="btn challenge-setup" data-lab="toy-lwe-success">Set up this challenge</button>
    <details class="selfcheck"><summary>Reveal answer</summary><p>Raising sigma shrinks the gap between the target vector (-s, e, 1) and ordinary lattice vectors, so reduction stops isolating it - the norm-gap meter falls and "Lattice attack result" flips to FAILED. Raising n both enlarges the embedding and raises the beta needed. Noise (sigma) usually bites fastest at fixed small n because it directly lengthens the target vector relative to q.</p></details>
  </div>

  <div class="challenge">
    <p><strong>6. BKZ-2 vs BKZ-6.</strong> Load the BKZ setup (seeded), run the attack at beta=6, then set beta=2 and re-run on the same seed. Compare tours, block improvements, and the norm gap.</p>
    <button type="button" class="btn challenge-setup" data-lab="bkz-helps">Set up this challenge</button>
    <details class="selfcheck"><summary>Reveal answer</summary><p>Because the seed is fixed, both runs attack the identical instance. BKZ-6 enumerates larger blocks, so it typically logs more block improvements and produces a shorter first vector (smaller norm gap) than BKZ-2 (which is essentially LLL). On these toy parameters both often still recover, but BKZ-6 reaches a better-reduced basis - the point real attacks exploit as dimensions grow.</p></details>
  </div>
</section>

</main>

`;

function byId<T extends HTMLElement>(id: string): T {
  const el = document.getElementById(id) as T | null;
  if (!el) throw new Error(`Missing element: ${id}`);
  return el;
}

function clamp(v: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, v));
}

function vec2FromPolar(angleDeg: number, len: number): Vec {
  const rad = (angleDeg * Math.PI) / 180;
  return [len * Math.cos(rad), len * Math.sin(rad)];
}

function randInt(min: number, max: number): number {
  return randomIntInclusive(min, max);
}

function determinant2(b1: Vec, b2: Vec): number {
  return b1[0] * b2[1] - b1[1] * b2[0];
}

function drawAxes(ctx: CanvasRenderingContext2D, w: number, h: number): void {
  ctx.strokeStyle = '#334';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(w / 2, 0);
  ctx.lineTo(w / 2, h);
  ctx.moveTo(0, h / 2);
  ctx.lineTo(w, h / 2);
  ctx.stroke();
}

function drawArrow(
  ctx: CanvasRenderingContext2D,
  origin: Vec,
  tip: Vec,
  color: string,
  dashed = false,
  width = 2,
): void {
  ctx.save();
  if (dashed) ctx.setLineDash([6, 4]);
  ctx.strokeStyle = color;
  ctx.fillStyle = color;
  ctx.lineWidth = width;
  ctx.beginPath();
  ctx.moveTo(origin[0], origin[1]);
  ctx.lineTo(tip[0], tip[1]);
  ctx.stroke();

  const ang = Math.atan2(tip[1] - origin[1], tip[0] - origin[0]);
  const ah = 8;
  ctx.beginPath();
  ctx.moveTo(tip[0], tip[1]);
  ctx.lineTo(tip[0] - ah * Math.cos(ang - 0.3), tip[1] - ah * Math.sin(ang - 0.3));
  ctx.lineTo(tip[0] - ah * Math.cos(ang + 0.3), tip[1] - ah * Math.sin(ang + 0.3));
  ctx.closePath();
  ctx.fill();
  ctx.restore();
}

function project(v: Vec): Vec {
  if (v.length === 2) return [v[0], v[1]];
  const x = v[0] - 0.45 * v[2];
  const y = v[1] - 0.25 * v[2];
  return [x, y];
}

function drawBasisScene(
  canvas: HTMLCanvasElement,
  original: Basis,
  current: Basis,
  gso?: Vec[],
  label?: string,
): void {
  const ctx = canvas.getContext('2d');
  if (!ctx) return;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  drawAxes(ctx, canvas.width, canvas.height);

  const center: Vec = [canvas.width / 2, canvas.height / 2];
  const scale = 22;

  const b1 = project(current[0]);
  const b2 = project(current[1]);
  for (let a = -5; a <= 5; a += 1) {
    for (let b = -5; b <= 5; b += 1) {
      const p: Vec = [
        center[0] + scale * (a * b1[0] + b * b2[0]),
        center[1] - scale * (a * b1[1] + b * b2[1]),
      ];
      ctx.fillStyle = '#1a3a6a';
      ctx.beginPath();
      ctx.arc(p[0], p[1], 2, 0, 2 * Math.PI);
      ctx.fill();
    }
  }

  ctx.fillStyle = '#ffd700';
  ctx.beginPath();
  ctx.arc(center[0], center[1], 4, 0, 2 * Math.PI);
  ctx.fill();

  for (const v of original) {
    const p = project(v);
    drawArrow(ctx, center, [center[0] + scale * p[0], center[1] - scale * p[1]], '#ff3366', false, 2);
  }
  for (const v of current) {
    const p = project(v);
    drawArrow(ctx, center, [center[0] + scale * p[0], center[1] - scale * p[1]], '#00ff88', false, 3);
  }
  if (gso) {
    for (const v of gso) {
      const p = project(v);
      drawArrow(ctx, center, [center[0] + scale * p[0], center[1] - scale * p[1]], '#00d4ff', true, 2);
    }
  }

  if (label) {
    ctx.fillStyle = '#ddd';
    ctx.font = '13px "JetBrains Mono", monospace';
    ctx.fillText(label, 10, 20);
  }
}

const themeToggle = byId<HTMLButtonElement>('theme-toggle');

function setupThemeToggle(button: HTMLButtonElement): void {
  const sync = () => {
    const mode = document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
    button.textContent = mode === 'dark' ? '🌙' : '☀️';
    button.setAttribute('aria-label', mode === 'dark' ? 'Switch to light mode' : 'Switch to dark mode');
    button.setAttribute('aria-pressed', mode === 'light' ? 'true' : 'false');
    button.title = mode === 'dark' ? 'Switch to light mode' : 'Switch to dark mode';
  };

  button.addEventListener('click', () => {
    const next = document.documentElement.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem('theme', next);
    sync();
  });

  sync();
}

setupThemeToggle(themeToggle);

const b1Angle = byId<HTMLInputElement>('b1-angle');
const b1Len = byId<HTMLInputElement>('b1-len');
const b2Angle = byId<HTMLInputElement>('b2-angle');
const b2Len = byId<HTMLInputElement>('b2-len');
const detLabel = byId<HTMLParagraphElement>('det-label');
const sameLatticeLabel = byId<HTMLParagraphElement>('same-lattice-label');
const latticeCanvas = byId<HTMLCanvasElement>('lattice-canvas');
let latticeBasis: Basis = [vec2FromPolar(25, 5), vec2FromPolar(95, 4)];

function updateSliderAria(el: HTMLInputElement, label: string): void {
  el.setAttribute('aria-valuenow', el.value);
  el.setAttribute('aria-valuetext', `${label} ${el.value}`);
}

function bindRangeAria(el: HTMLInputElement, label: string): void {
  const update = () => updateSliderAria(el, label);
  el.addEventListener('input', update);
  update();
}

function renderEx1(): void {
  latticeBasis = [
    vec2FromPolar(Number(b1Angle.value), Number(b1Len.value)),
    vec2FromPolar(Number(b2Angle.value), Number(b2Len.value)),
  ];
  const det = Math.abs(determinant2(latticeBasis[0], latticeBasis[1]));
  detLabel.textContent = `det(Lambda) = ${det.toFixed(3)}`;
  latticeCanvas.setAttribute('aria-label', `Lattice canvas with determinant ${det.toFixed(3)}`);
  drawBasisScene(latticeCanvas, latticeBasis, latticeBasis, undefined, 'Original basis in red');
  updateSliderAria(b1Angle, 'b1 angle');
  updateSliderAria(b1Len, 'b1 length');
  updateSliderAria(b2Angle, 'b2 angle');
  updateSliderAria(b2Len, 'b2 length');
}

[b1Angle, b1Len, b2Angle, b2Len].forEach((el) => el.addEventListener('input', renderEx1));
bindRangeAria(b1Angle, 'b1 angle');
bindRangeAria(b1Len, 'b1 length');
bindRangeAria(b2Angle, 'b2 angle');
bindRangeAria(b2Len, 'b2 length');

const unimodularChoices = [
  [
    [1, 1],
    [0, 1],
  ],
  [
    [1, 0],
    [1, 1],
  ],
  [
    [1, -1],
    [0, 1],
  ],
  [
    [0, 1],
    [1, 0],
  ],
  [
    [-1, 0],
    [0, 1],
  ],
] as const;

byId<HTMLButtonElement>('same-lattice-btn').addEventListener('click', () => {
  const t = unimodularChoices[randInt(0, unimodularChoices.length - 1)]!;
  const b1 = latticeBasis[0];
  const b2 = latticeBasis[1];
  const nb1: Vec = [t[0][0] * b1[0] + t[0][1] * b2[0], t[0][0] * b1[1] + t[0][1] * b2[1]];
  const nb2: Vec = [t[1][0] * b1[0] + t[1][1] * b2[0], t[1][0] * b1[1] + t[1][1] * b2[1]];
  drawBasisScene(latticeCanvas, [nb1, nb2], [nb1, nb2], undefined, 'Same lattice. Same det. Uglier basis.');
  sameLatticeLabel.textContent = 'Same lattice. Same det. Uglier basis.';
});

const gsoInputs = ['gso-a11', 'gso-a12', 'gso-a21', 'gso-a22'].map((id) => byId<HTMLInputElement>(id));
const gsoNext = byId<HTMLButtonElement>('gso-next');
const gsoText = byId<HTMLPreElement>('gso-text');
const gsoCanvas = byId<HTMLCanvasElement>('gso-canvas');
let gsoStep = 0;

function readGsoBasis(): Basis {
  const a11 = Number(gsoInputs[0].value);
  const a12 = Number(gsoInputs[1].value);
  const a21 = Number(gsoInputs[2].value);
  const a22 = Number(gsoInputs[3].value);
  return [
    [a11, a12],
    [a21, a22],
  ];
}

function renderEx2(): void {
  const basis = readGsoBasis();
  const { gso, mu } = gramSchmidt(basis);
  const mu21 = mu[1][0];
  const lhs = 0.75 * dot(gso[0], gso[0]);
  const rhs = dot([gso[1][0] + mu21 * gso[0][0], gso[1][1] + mu21 * gso[0][1]], [gso[1][0] + mu21 * gso[0][0], gso[1][1] + mu21 * gso[0][1]]);
  const lovaszOk = lhs <= rhs;

  let text = '';
  if (gsoStep === 0) {
    text = `Step 1\nb1~ = [${gso[0].map((x) => x.toFixed(3)).join(', ')}]\n||b1~|| = ${norm(gso[0]).toFixed(3)}`;
  } else if (gsoStep === 1) {
    text = [
      'Step 2',
      `mu21 = ${mu21.toFixed(6)}`,
      `b2~ = [${gso[1].map((x) => x.toFixed(3)).join(', ')}]`,
      `||b2~|| = ${norm(gso[1]).toFixed(3)}`,
    ].join('\n');
  } else {
    text = [
      'Step 3 - Lovasz check',
      `0.75 * ||b1~||^2 = ${lhs.toFixed(4)}`,
      `||b2~ + mu21*b1~||^2 = ${rhs.toFixed(4)}`,
      lovaszOk ? 'Condition holds: advance.' : 'Condition fails: swap needed.',
    ].join('\n');
  }
  gsoText.textContent = text;
  drawBasisScene(gsoCanvas, basis, basis, gso, lovaszOk ? 'Lovasz pass' : 'Lovasz fail');
}

gsoNext.addEventListener('click', () => {
  gsoStep = (gsoStep + 1) % 3;
  renderEx2();
});
gsoInputs.forEach((i) => i.addEventListener('input', () => {
  gsoStep = 0;
  renderEx2();
}));

const lllDim = byId<HTMLSelectElement>('lll-dim');
const lllDelta = byId<HTMLInputElement>('lll-delta');
const lllMatrix = byId<HTMLTextAreaElement>('lll-matrix');
const lllLog = byId<HTMLDivElement>('lll-log');
const lllMetrics = byId<HTMLPreElement>('lll-metrics');
const lllCanvas = byId<HTMLCanvasElement>('lll-canvas');
const defectCanvas = byId<HTMLCanvasElement>('defect-chart');
let lllSteps = lllReduce([
  [19, 2],
  [7, 1],
]).steps;
let lllCursor = 0;
let lllAutoTimer: number | null = null;
let lllOriginal: Basis = [
  [19, 2],
  [7, 1],
];

function parseMatrix(text: string): Basis {
  const rows = text
    .trim()
    .split('\n')
    .map((row) => row.trim())
    .filter((row) => row.length > 0)
    .map((row) => row.split(/\s+/).map((x) => Number(x)));
  if (rows.length === 0) throw new Error('Empty matrix');
  const n = rows.length;
  if (!rows.every((r) => r.length === n)) throw new Error('Matrix must be square');
  if (!rows.every((r) => r.every((x) => Number.isFinite(x) && Number.isInteger(x)))) {
    throw new Error('Matrix values must be integers');
  }
  return rows;
}

function setPreset(name: string): void {
  if (name === 'random2d') {
    lllDim.value = '2';
    lllMatrix.value = `${randInt(8, 30)} ${randInt(1, 7)}\n${randInt(3, 12)} ${randInt(1, 5)}`;
  } else if (name === 'classic') {
    lllDim.value = '2';
    lllMatrix.value = '19 2\n7 1';
  } else if (name === 'near') {
    lllDim.value = '2';
    lllMatrix.value = '10 0\n1 9';
  } else if (name === 'challenge3d') {
    lllDim.value = '3';
    lllMatrix.value = '19 2 3\n7 1 6\n2 9 5';
  } else {
    lllDim.value = '3';
    lllMatrix.value = '71 0 0\n17 1 0\n63 0 1';
  }
  resetLLL();
}

function drawDefectChart(values: number[]): void {
  const ctx = defectCanvas.getContext('2d');
  if (!ctx) return;
  ctx.clearRect(0, 0, defectCanvas.width, defectCanvas.height);
  if (values.length === 0) return;
  const finite = values.filter((v) => Number.isFinite(v));
  if (finite.length === 0) return;
  const minV = Math.min(...finite);
  const maxV = Math.max(...finite);
  const span = Math.max(maxV - minV, 1e-9);
  ctx.strokeStyle = '#00d4ff';
  ctx.lineWidth = 2;
  ctx.beginPath();
  values.forEach((v, i) => {
    const x = (i / Math.max(values.length - 1, 1)) * (defectCanvas.width - 20) + 10;
    const y = defectCanvas.height - 10 - ((v - minV) / span) * (defectCanvas.height - 20);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  });
  ctx.stroke();
}

function renderLLLState(): void {
  const step = lllSteps[Math.min(lllCursor, lllSteps.length - 1)]!;
  drawBasisScene(lllCanvas, lllOriginal, step.after, step.gsoAfter, step.description);

  lllLog.textContent = lllSteps
    .slice(0, lllCursor + 1)
    .map((s, idx) => `${idx + 1}. ${s.description}`)
    .join('\n');

  const initial = orthogonalityDefect(lllOriginal);
  const current = orthogonalityDefect(step.after);
  const shortest = shortestVectorInBasis(step.after);
  const defects = lllSteps.slice(0, lllCursor + 1).map((s) => orthogonalityDefect(s.after));
  drawDefectChart(defects);

  // The invariant, made visible: det of the lattice is unchanged and the running
  // transform U (with reducedBasis = U * original) is unimodular at every step.
  const detOriginal = Math.abs(determinant(lllOriginal));
  const detNow = Math.abs(determinant(step.after));
  const detU = Math.round(determinant(step.transformAfter));
  const uRows = step.transformAfter.map((r) => `[${r.join(', ')}]`).join('  ');

  lllMetrics.textContent = [
    `Orthogonality defect: ${initial.toFixed(4)} -> ${current.toFixed(4)}`,
    `Basis norms: [${step.after.map((v) => norm(v).toFixed(3)).join(', ')}]`,
    `GSO norms: [${step.gsoAfter.map((v) => norm(v).toFixed(3)).join(', ')}]`,
    `Swap count: ${lllSteps.filter((s) => s.type === 'swap').length}`,
    `Step count: ${lllCursor + 1} / ${lllSteps.length}`,
    `Shortest basis vector: [${shortest.join(', ')}], norm=${norm(shortest).toFixed(4)}`,
    `Lattice det |det B|: ${detOriginal.toFixed(3)} -> ${detNow.toFixed(3)} (unchanged: same lattice)`,
    `Transform U (reduced = U x original), det(U)=${detU} ${Math.abs(detU) === 1 ? '(unimodular)' : ''}`,
    `U = ${uRows}`,
  ].join('\n');
}

function resetLLL(): void {
  try {
    if (lllAutoTimer !== null) {
      window.clearInterval(lllAutoTimer);
      lllAutoTimer = null;
    }
    lllOriginal = parseMatrix(lllMatrix.value);
    const delta = Number(lllDelta.value);
    lllSteps = lllReduce(lllOriginal, delta).steps;
    lllCursor = 0;
    renderLLLState();
  } catch (err) {
    lllLog.textContent = `Input error: ${(err as Error).message}`;
  }
}

byId<HTMLButtonElement>('lll-step').addEventListener('click', () => {
  if (lllCursor < lllSteps.length - 1) {
    lllCursor += 1;
    renderLLLState();
  }
});

byId<HTMLButtonElement>('lll-auto').addEventListener('click', () => {
  if (lllAutoTimer !== null) {
    window.clearInterval(lllAutoTimer);
    lllAutoTimer = null;
    return;
  }
  lllAutoTimer = window.setInterval(() => {
    if (lllCursor >= lllSteps.length - 1) {
      if (lllAutoTimer !== null) window.clearInterval(lllAutoTimer);
      lllAutoTimer = null;
      return;
    }
    lllCursor += 1;
    renderLLLState();
  }, 500);
});

byId<HTMLButtonElement>('lll-reset').addEventListener('click', resetLLL);
lllDelta.addEventListener('input', () => {
  updateSliderAria(lllDelta, 'delta');
});
bindRangeAria(lllDelta, 'delta');
lllDim.addEventListener('change', () => {
  if (lllDim.value === '2') {
    lllMatrix.value = '19 2\n7 1';
  } else {
    lllMatrix.value = '19 2 3\n7 1 6\n2 9 5';
  }
  resetLLL();
});

document.querySelectorAll<HTMLButtonElement>('.preset').forEach((btn) => {
  btn.addEventListener('click', () => setPreset(btn.dataset.preset ?? 'classic'));
});

const lweN = byId<HTMLInputElement>('lwe-n');
const lweQ = byId<HTMLInputElement>('lwe-q');
const lweSigma = byId<HTMLInputElement>('lwe-sigma');
const lweBeta = byId<HTMLInputElement>('lwe-beta');
const lweOutput = byId<HTMLPreElement>('lwe-output');
const lweProgress = byId<HTMLProgressElement>('lwe-progress');
const lweMeterText = byId<HTMLDivElement>('lwe-meter-text');
const lweMeterFill = byId<HTMLDivElement>('lwe-meter-fill');
// Keep the parameters the instance was generated with. The attack must use these,
// not the live slider values -- moving a slider after Generate would otherwise
// build the embedding with a mismatched q/n while b still encodes the old ones,
// silently breaking (or falsifying) recovery.
type LWEInstance = Awaited<ReturnType<typeof generateLWE>> & {
  n: number;
  q: number;
  sigma: number;
};
let lweInstance: LWEInstance | null = null;

function updateLWEMeter(
  score: number,
  message: string,
  color: 'good' | 'warn' | 'bad',
): void {
  const clamped = clamp(score, 0, 100);
  lweMeterFill.style.width = `${clamped.toFixed(1)}%`;
  if (color === 'good') {
    lweMeterFill.style.background = 'linear-gradient(90deg, #00ff88, #66ffb5)';
  } else if (color === 'warn') {
    lweMeterFill.style.background = 'linear-gradient(90deg, #ffaa00, #ffd166)';
  } else {
    lweMeterFill.style.background = 'linear-gradient(90deg, #ff3366, #ff6b8f)';
  }
  lweMeterText.textContent = `Norm-gap confidence: ${clamped.toFixed(1)}% - ${message}`;
}

function estimateTargetNorm(n: number, m: number, sigma: number): number {
  // The target short vector is (-s, e, 1). generateLWE draws each secret
  // coordinate uniformly from {0..4}, so E[s_i^2] = (0+1+4+9+16)/5 = 6; the
  // secret block contributes ~6n, not n. Using n here made the gap meter
  // systematically pessimistic.
  const secretEnergy = 6 * n;
  return Math.sqrt(secretEnergy + m * sigma * sigma + 1);
}

function scoreFromNormGap(shortNorm: number | null, targetNorm: number): number {
  if (shortNorm === null || !Number.isFinite(shortNorm) || targetNorm <= 0) {
    return 0;
  }
  const gap = shortNorm / targetNorm;
  if (!Number.isFinite(gap) || gap <= 1) {
    return 100;
  }
  const decay = 18;
  const score = 100 * Math.exp(-(gap - 1) / decay);
  return clamp(score, 0, 100);
}

function summarizeBKZImpact(
  tours: number,
  improvements: number,
  success: boolean,
  recoverMethod: 'short-vector' | 'bruteforce' | 'none',
): string {
  if (improvements === 0) {
    return 'Why this tour mattered: BKZ made no block improvements, so the basis stayed too coarse to expose a useful secret vector.';
  }
  if (!success) {
    return `Why this tour mattered: BKZ improved ${improvements} block(s) over ${tours} tour(s), but not enough to isolate a target-shaped short vector.`;
  }
  if (recoverMethod === 'short-vector') {
    return `Why this tour mattered: BKZ improvements compressed the basis enough to expose a short vector consistent with (-s, e, +-1).`;
  }
  return `Why this tour mattered: BKZ improved the basis structure, then the toy fallback recovered s from the low-noise residual pattern.`;
}

function fmtMatrix(m: number[][], maxRows = 6): string {
  const rows = m.slice(0, maxRows).map((r) => `  [${r.join(', ')}]`);
  if (m.length > maxRows) rows.push('  ...');
  return rows.join('\n');
}

byId<HTMLButtonElement>('lwe-generate').addEventListener('click', async () => {
  const n = Number(lweN.value);
  const q = Number(lweQ.value);
  const sigma = Number(lweSigma.value);
  updateSliderAria(lweN, 'n');
  updateSliderAria(lweQ, 'q');
  updateSliderAria(lweSigma, 'sigma');
  updateSliderAria(lweBeta, 'beta');
  lweProgress.value = 15;
  // When a seed is active, rewind the deterministic stream before generating so
  // the same seed + parameters always reproduce the exact same instance.
  const seed = getSeed();
  if (seed !== null) setSeed(seed);
  lweInstance = { ...(await generateLWE(n, q, sigma)), n, q, sigma };
  lweProgress.value = 35;
  const m = lweInstance.A.length;
  const lattice = buildLWELattice(lweInstance.A, lweInstance.b, q);
  lweOutput.textContent = [
    `Secret s = [${lweInstance.secret.join(', ')}]`,
    `Samples m = ${m}`,
    `A (${m}x${n}):`,
    fmtMatrix(lweInstance.A, 5),
    `b = [${lweInstance.b.join(', ')}]`,
    `hidden errors e = [${lweInstance.errors.join(', ')}]`,
    `Embedding lattice (${lattice.length}x${lattice[0].length}):`,
    fmtMatrix(lattice, 5),
  ].join('\n');
  updateLWEMeter(0, 'instance generated, run attack to measure gap.', 'warn');
  lweProgress.value = 40;
});

byId<HTMLButtonElement>('lwe-attack').addEventListener('click', () => {
  if (!lweInstance) {
    lweOutput.textContent = 'Generate an LWE instance first.';
    return;
  }
  try {
    // n/q/sigma come from the generated instance, not the live sliders -- see the
    // LWEInstance note. beta is the one genuine attack-time choice, so read it live.
    const { n, q, sigma } = lweInstance;
    const beta = Number(lweBeta.value);
    lweProgress.value = 55;
    const lattice = buildLWELattice(lweInstance.A, lweInstance.b, q);
    const bkz = bkzReduce(lattice, beta);
    lweProgress.value = 80;
    // Recover from the SAME basis we just reduced and report, so the beta slider
    // genuinely drives the attack instead of a hidden second reduction.
    const result =
      n >= 20
        ? { success: false, recovered: null, shortVec: null, recoverMethod: 'none' as const }
        : recoverFromReducedBasis(bkz.reducedBasis, lweInstance.A, lweInstance.b, q, n, sigma);
    const ok = result.success && result.recovered !== null;
    const shortNorm = result.shortVec ? norm(result.shortVec) : null;
    const targetNorm = estimateTargetNorm(n, lweInstance.A.length, sigma);
    const gapScore = scoreFromNormGap(shortNorm, targetNorm);
    lweOutput.textContent += `\n\nRunning BKZ-${beta} on ${lattice.length}-dim lattice...`;
    lweOutput.textContent += `\nTours: ${bkz.tours}, improvements: ${bkz.improvements}`;
    if (bkz.tourLogs.length > 0) {
      lweOutput.textContent += '\nTour log:';
      for (const t of bkz.tourLogs) {
        lweOutput.textContent += `\n  tour ${t.tour}: ${t.improvementsInTour} improvements${t.converged ? ' (converged)' : ''}`;
      }
    }
    if (bkz.blockImprovements.length > 0) {
      lweOutput.textContent += '\nBlock improvements:';
      for (const imp of bkz.blockImprovements.slice(0, 8)) {
        lweOutput.textContent += `\n  [${imp.blockStart}..${imp.blockEnd}] ||b0|| ${imp.headNormBefore.toFixed(3)} -> ${imp.insertedNorm.toFixed(3)}`;
      }
      if (bkz.blockImprovements.length > 8) {
        lweOutput.textContent += `\n  ... ${bkz.blockImprovements.length - 8} more`;
      }
    }
    if (result.shortVec) {
      lweOutput.textContent += `\nShortest basis vector found: [${result.shortVec.join(', ')}]`;
      lweOutput.textContent += `\nNorm: ${shortNorm!.toFixed(4)}`;
      lweOutput.textContent += `\nTarget-like norm estimate: ${targetNorm.toFixed(4)} (gap ${(shortNorm! / targetNorm).toFixed(2)}x)`;
      // Reveal the primal-embedding structure of the recovered vector so the
      // lesson is visible, not asserted: the LWE vector is (-s, e, +-1). When
      // reduction surfaces that vector, its three blocks line up with the
      // secret, the error, and the embedding coordinate.
      const v = result.shortVec;
      const m = lweInstance.A.length;
      const embed = v[v.length - 1];
      if (Math.abs(embed) === 1 && v.length === n + m + 1) {
        const sign = -embed;
        const sBlock = v.slice(0, n).map((x) => sign * x);
        const eBlock = v.slice(n, n + m).map((x) => sign * x);
        lweOutput.textContent += '\nStructure (-s, e, +-1) read off this vector:';
        lweOutput.textContent += `\n  secret block -> s = [${sBlock.join(', ')}]`;
        lweOutput.textContent += `\n  error  block -> e = [${eBlock.join(', ')}]`;
        lweOutput.textContent += `\n  embed coord  -> ${embed} (anchors the shift)`;
      }
    }
    // Keep the lattice attack and the teaching baseline strictly separated, so a
    // baseline recovery is never mistaken for a lattice break.
    const latticeOk = ok && result.recoverMethod === 'short-vector';
    lweOutput.textContent += `\n\nLattice attack result: ${latticeOk ? 'SUCCESS - secret read off a reduced short vector' : 'FAILED - reduction did not expose the (-s, e, +-1) vector'}`;
    if (result.recoverMethod === 'bruteforce' && result.recovered) {
      lweOutput.textContent += '\nTeaching baseline: brute force found the secret (exhaustive search over the tiny toy secret space)';
      lweOutput.textContent += '\n  -> do NOT count baseline recovery as a lattice break';
    } else if (!latticeOk) {
      lweOutput.textContent +=
        n > 8
          ? '\nTeaching baseline: skipped (n too large for brute force)'
          : '\nTeaching baseline: brute force did not recover the secret either';
    }
    if (ok) {
      const recovered = result.recovered!;
      const exact = recovered.join(',') === lweInstance.secret.join(',');
      lweOutput.textContent += `\nRecovered secret: [${recovered.join(', ')}] ${exact ? 'EXACT MATCH' : 'close'}`;
    }
    lweOutput.textContent += `\n${summarizeBKZImpact(bkz.tours, bkz.improvements, ok, result.recoverMethod)}`;
    if (ok && result.recoverMethod === 'short-vector') {
      updateLWEMeter(gapScore, 'short vector is close to target regime.', 'good');
    } else if (ok) {
      updateLWEMeter(gapScore, 'basis improved, but recovery leaned on toy fallback.', 'warn');
    } else {
      updateLWEMeter(gapScore, 'short vector remains far from target regime.', 'bad');
    }
    lweProgress.value = 100;
  } catch (err) {
    const message = err instanceof Error ? err.message : 'Unknown attack error';
    lweOutput.textContent += `\nAttack error: ${message}`;
    updateLWEMeter(0, 'attack run errored before convergence.', 'bad');
    lweProgress.value = 0;
  }
});

byId<HTMLButtonElement>('kyber-try').addEventListener('click', () => {
  const n = 256;
  const q = 3329;
  const sigma = 1;
  const dim = n + 1 + 256;
  const reqBeta = 400;
  const costExp = 0.292 * reqBeta;
  lweOutput.textContent += `\n\nKyber-like test:`;
  lweOutput.textContent += `\nn=${n}, q=${q}, sigma=${sigma}, embedding dimension~${dim}`;
  lweOutput.textContent += '\nBKZ-2 (LLL) does not recover target-length vectors in this regime.';
  lweOutput.textContent += `\nRequired block size beta~${reqBeta}, cost~2^${costExp.toFixed(0)} operations.`;
  lweOutput.textContent += '\nThis demonstrates why real Kyber parameters resist LLL/BKZ at practical resources.';
  updateLWEMeter(1, 'Kyber-scale regime: norm gap is overwhelmingly large.', 'bad');
});

const exploreN = byId<HTMLInputElement>('explore-n');
const exploreLeft = byId<HTMLPreElement>('explore-left');
const exploreRight = byId<HTMLPreElement>('explore-right');
const thresholdCanvas = byId<HTMLCanvasElement>('threshold-chart');

function drawThresholdChart(nNow: number): void {
  const ctx = thresholdCanvas.getContext('2d');
  if (!ctx) return;
  ctx.clearRect(0, 0, thresholdCanvas.width, thresholdCanvas.height);
  ctx.strokeStyle = '#345';
  ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i += 1) {
    const y = 20 + i * 55;
    ctx.beginPath();
    ctx.moveTo(40, y);
    ctx.lineTo(thresholdCanvas.width - 20, y);
    ctx.stroke();
  }

  ctx.strokeStyle = '#00ff88';
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let n = 4; n <= 256; n += 1) {
    const beta = n / (2 * Math.log2(3329 / 1));
    const cost = 0.292 * beta;
    const x = 40 + ((n - 4) / (256 - 4)) * (thresholdCanvas.width - 60);
    const y = thresholdCanvas.height - 20 - clamp(cost / 140, 0, 1) * (thresholdCanvas.height - 40);
    if (n === 4) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  const xNow = 40 + ((nNow - 4) / (256 - 4)) * (thresholdCanvas.width - 60);
  ctx.strokeStyle = '#ffd700';
  ctx.beginPath();
  ctx.moveTo(xNow, 20);
  ctx.lineTo(xNow, thresholdCanvas.height - 20);
  ctx.stroke();

  const xThreshold = 40 + ((50 - 4) / (256 - 4)) * (thresholdCanvas.width - 60);
  ctx.strokeStyle = '#ffaa00';
  ctx.beginPath();
  ctx.moveTo(xThreshold, 20);
  ctx.lineTo(xThreshold, thresholdCanvas.height - 20);
  ctx.stroke();
  ctx.fillStyle = '#ccc';
  ctx.font = '12px "JetBrains Mono", monospace';
  ctx.fillText('Security threshold ~ n=50', xThreshold + 8, 34);
}

function renderEx5(): void {
  const n = Number(exploreN.value);
  const q = 3329;
  const sigma = 1;
  const approxFactorExp = n / 2;
  const beta = n / (2 * Math.log2(q / sigma));
  const costExp = 0.292 * beta;
  const insecure = n <= 20;
  exploreLeft.textContent = [
    'INSECURE (toy-ish)',
    `n=${Math.min(n, 16)}, q=101, sigma=3`,
    `SVP approx factor 2^(n/2)=2^${(Math.min(n, 16) / 2).toFixed(1)}`,
    `Required beta~${(Math.min(n, 16) / (2 * Math.log2(101 / 3))).toFixed(1)}`,
    'Status: BROKEN by LLL/BKZ-2 for tiny dimensions',
  ].join('\n');
  exploreRight.textContent = [
    'SECURE (Kyber-512 style)',
    `n=${n}, q=${q}, sigma=${sigma}`,
    `SVP approx factor: 2^${approxFactorExp.toFixed(1)}`,
    `Estimated beta~${beta.toFixed(2)}`,
    `Attack cost~2^${costExp.toFixed(2)}`,
    insecure ? 'Status: still toy scale' : 'Status: practical attacks out of reach',
  ].join('\n');
  exploreN.style.accentColor = n < 50 ? '#ff3366' : n < 100 ? '#ffaa00' : '#00ff88';
  drawThresholdChart(n);
  updateSliderAria(exploreN, 'parameter n');
}

exploreN.addEventListener('input', renderEx5);
bindRangeAria(lweN, 'n');
bindRangeAria(lweQ, 'q');
bindRangeAria(lweSigma, 'sigma');
bindRangeAria(lweBeta, 'beta');
bindRangeAria(exploreN, 'parameter n');

// ---- Reproducible labs: seeding, shareable links, canonical labs ----
const labSeed = byId<HTMLInputElement>('lab-seed');
const labPreset = byId<HTMLSelectElement>('lab-preset');
const labStatus = byId<HTMLParagraphElement>('lab-status');

function updateLabStatus(extra?: string): void {
  const seed = getSeed();
  const base =
    seed !== null
      ? `Mode: seeded (seed=${seed}) - reproducible: same seed + parameters give the same instance.`
      : 'Mode: crypto-random (not reproducible). Enter a seed for repeatable runs.';
  labStatus.textContent = extra ? `${base} ${extra}` : base;
}

function applySeedFromInput(): void {
  const raw = labSeed.value.trim();
  if (raw === '') {
    setSeed(null);
  } else {
    const parsed = Number(raw);
    setSeed(Number.isFinite(parsed) ? parsed : null);
    if (!Number.isFinite(parsed)) labSeed.value = '';
  }
  updateLabStatus();
}

function refreshLWEAria(): void {
  updateSliderAria(lweN, 'n');
  updateSliderAria(lweQ, 'q');
  updateSliderAria(lweSigma, 'sigma');
  updateSliderAria(lweBeta, 'beta');
}

function buildLabQuery(): string {
  const p = new URLSearchParams();
  const seed = getSeed();
  if (seed !== null) p.set('seed', String(seed));
  p.set('n', lweN.value);
  p.set('q', lweQ.value);
  p.set('sigma', lweSigma.value);
  p.set('beta', lweBeta.value);
  p.set('dim', lllDim.value);
  p.set('delta', lllDelta.value);
  p.set('matrix', lllMatrix.value);
  p.set('en', exploreN.value);
  return p.toString();
}

function applyLabQuery(query: string): void {
  const p = new URLSearchParams(query);
  if (p.has('seed')) {
    const s = Number(p.get('seed'));
    if (Number.isFinite(s)) {
      labSeed.value = String(Math.trunc(s));
      setSeed(s);
    }
  }
  const setIf = (key: string, el: HTMLInputElement | HTMLSelectElement | HTMLTextAreaElement) => {
    const v = p.get(key);
    if (v !== null) el.value = v;
  };
  setIf('n', lweN);
  setIf('q', lweQ);
  setIf('sigma', lweSigma);
  setIf('beta', lweBeta);
  setIf('dim', lllDim);
  setIf('delta', lllDelta);
  setIf('matrix', lllMatrix);
  setIf('en', exploreN);
  refreshLWEAria();
  updateSliderAria(lllDelta, 'delta');
  resetLLL();
  renderEx5();
  updateLabStatus();
}

byId<HTMLButtonElement>('lab-apply-seed').addEventListener('click', applySeedFromInput);
labSeed.addEventListener('change', applySeedFromInput);

byId<HTMLButtonElement>('lab-copy-link').addEventListener('click', async () => {
  const url = `${location.origin}${location.pathname}#${buildLabQuery()}`;
  try {
    location.hash = buildLabQuery();
  } catch {
    /* ignore hash write failures */
  }
  try {
    await navigator.clipboard.writeText(url);
    updateLabStatus('Lab link copied to clipboard.');
  } catch {
    updateLabStatus(`Copy failed - link is in the address bar.`);
  }
});

interface CanonicalLab {
  seed: number;
  focus: string;
  note: string;
  apply: () => void;
}

const CANONICAL_LABS: Record<string, CanonicalLab> = {
  'swap-cascade': {
    seed: 7,
    focus: 'exhibit-3',
    note: 'Use Step to watch each swap; the lattice determinant never changes.',
    apply: () => {
      lllDim.value = '2';
      lllDelta.value = '0.99';
      lllMatrix.value = '233 144\n144 89';
      updateSliderAria(lllDelta, 'delta');
      resetLLL();
    },
  },
  'delta-sensitivity': {
    seed: 1,
    focus: 'exhibit-3',
    note: 'Run as-is (delta=0.5), then raise delta toward 0.99 and Reset to compare swap counts.',
    apply: () => {
      lllDim.value = '2';
      lllDelta.value = '0.5';
      lllMatrix.value = '40 1\n39 1';
      updateSliderAria(lllDelta, 'delta');
      resetLLL();
    },
  },
  'toy-lwe-success': {
    seed: 1234,
    focus: 'exhibit-4',
    note: 'Instance generated. Click "Run LLL/BKZ Attack" - the secret is read off the short vector.',
    apply: () => {
      lweN.value = '4';
      lweQ.value = '71';
      lweSigma.value = '2';
      lweBeta.value = '2';
      refreshLWEAria();
      byId<HTMLButtonElement>('lwe-generate').click();
    },
  },
  'bkz-helps': {
    seed: 2024,
    focus: 'exhibit-4',
    note: 'Generated at beta=6. Run the attack, then lower beta to 2 and re-run on the same seed to compare.',
    apply: () => {
      lweN.value = '6';
      lweQ.value = '101';
      lweSigma.value = '2';
      lweBeta.value = '6';
      refreshLWEAria();
      byId<HTMLButtonElement>('lwe-generate').click();
    },
  },
  'kyber-failure': {
    seed: 0,
    focus: 'exhibit-5',
    note: 'In the Parameter Explorer, drag n toward 256 - required beta and attack cost explode out of reach.',
    apply: () => {
      exploreN.value = '256';
      renderEx5();
    },
  },
};

labPreset.addEventListener('change', () => {
  const lab = CANONICAL_LABS[labPreset.value];
  if (!lab) return;
  labSeed.value = String(lab.seed);
  setSeed(lab.seed);
  lab.apply();
  updateLabStatus(lab.note);
  document.getElementById(lab.focus)?.scrollIntoView({ block: 'start' });
});

// Challenge "Set up" buttons reuse the canonical-lab machinery: select the lab
// and fire its change handler so each challenge loads a reproducible setup.
document.querySelectorAll<HTMLButtonElement>('.challenge-setup').forEach((btn) => {
  btn.addEventListener('click', () => {
    const lab = btn.dataset.lab;
    if (!lab || !(lab in CANONICAL_LABS)) return;
    labPreset.value = lab;
    labPreset.dispatchEvent(new Event('change'));
  });
});

renderEx1();
renderEx2();
updateSliderAria(lllDelta, 'delta');
resetLLL();
renderEx5();

// Apply a shared lab link if present, otherwise show the default status.
if (location.hash.length > 1) {
  applyLabQuery(location.hash.slice(1));
} else {
  updateLabStatus();
}
