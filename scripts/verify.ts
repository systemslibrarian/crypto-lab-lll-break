import { dot, gramSchmidt, norm, type Basis } from '../src/lattice';
import { applyTransform, determinant, lllReduce } from '../src/lll';
import { attackLWE, buildLWELattice, generateLWE } from '../src/bkz';
import { setSeed } from '../src/rng';

function assertTrue(condition: boolean, message: string): void {
  if (!condition) {
    throw new Error(message);
  }
}

function approxEqual(a: number, b: number, eps = 1e-9): boolean {
  return Math.abs(a - b) <= eps;
}

function centeredMod(x: number, q: number): number {
  const r = ((x % q) + q) % q;
  return r > q / 2 ? r - q : r;
}

async function run(): Promise<void> {
  const gs = gramSchmidt([
    [3, 1],
    [1, 2],
  ]);
  assertTrue(approxEqual(gs.mu[1][0], 0.5), `Expected mu21=0.5, got ${gs.mu[1][0]}`);

  assertTrue(dot([1, 2], [3, 4]) === 11, 'dot product sanity check failed');
  assertTrue(approxEqual(norm([3, 4]), 5), 'norm sanity check failed');

  const lll = lllReduce([
    [19, 2],
    [7, 1],
  ]);
  const shortestNormSq = Math.min(...lll.reducedBasis.map((v) => dot(v, v)));
  assertTrue(shortestNormSq <= 5, `Expected short vector norm^2 <= 5, got ${shortestNormSq}`);

  // The headline LLL invariant, made provable: the reduction is a UNIMODULAR
  // change of basis, so reducedBasis = U * inputBasis exactly, |det(U)| = 1, and
  // the lattice determinant is unchanged. Check across several bases (including
  // a 3D one and an already-reduced one).
  const transformCases: Basis[] = [
    [
      [19, 2],
      [7, 1],
    ],
    [
      [1, 0],
      [0, 1],
    ],
    [
      [19, 2, 3],
      [7, 1, 6],
      [2, 9, 5],
    ],
  ];
  for (const original of transformCases) {
    const res = lllReduce(original);
    const reconstructed = applyTransform(res.transform, original);
    for (let r = 0; r < reconstructed.length; r += 1) {
      for (let c = 0; c < reconstructed[r].length; c += 1) {
        assertTrue(
          reconstructed[r][c] === res.reducedBasis[r][c],
          `U*B mismatch at [${r}][${c}]: ${reconstructed[r][c]} != ${res.reducedBasis[r][c]}`,
        );
      }
    }
    const detU = Math.round(determinant(res.transform));
    assertTrue(Math.abs(detU) === 1, `Transform must be unimodular, |det(U)|=${Math.abs(detU)}`);
    const detBefore = Math.abs(determinant(original));
    const detAfter = Math.abs(determinant(res.reducedBasis));
    assertTrue(
      approxEqual(detBefore, detAfter, 1e-6),
      `Lattice determinant changed: ${detBefore} -> ${detAfter}`,
    );
  }

  // The primal embedding must genuinely contain the (-s, e, 1) vector. Build it
  // by the exact integer combination the construction promises and confirm it is
  // a real lattice point with that shape -- if this fails, the demo is teaching a
  // broken embedding no matter what the attack happens to recover.
  {
    const q = 71;
    const inst = await generateLWE(4, q, 2);
    const L = buildLWELattice(inst.A, inst.b, q);
    const n = 4;
    const m = inst.A.length;
    const dim = n + m + 1;
    // coeffs: +1 on the last (embedding) row, -s_i on row i.
    const coeffs = new Array<number>(L.length).fill(0);
    coeffs[L.length - 1] = 1;
    for (let i = 0; i < n; i += 1) coeffs[i] = -inst.secret[i];
    const combo = new Array<number>(dim).fill(0);
    for (let r = 0; r < L.length; r += 1) {
      for (let c = 0; c < dim; c += 1) combo[c] += coeffs[r] * L[r][c];
    }
    // First n coords must be -s; embed coord must be 1; middle block must reduce
    // (mod q) to the small error vector e.
    for (let i = 0; i < n; i += 1) {
      assertTrue(combo[i] === -inst.secret[i], `Embedding secret coord ${i} wrong: ${combo[i]} != ${-inst.secret[i]}`);
    }
    assertTrue(combo[dim - 1] === 1, `Embedding coord must be 1, got ${combo[dim - 1]}`);
    for (let r = 0; r < m; r += 1) {
      assertTrue(
        centeredMod(combo[n + r] - inst.errors[r], q) === 0,
        `Embedding error coord ${r} not congruent to e mod q`,
      );
    }
  }

  // The toy attack must succeed BY LATTICE REDUCTION, not by the exhaustive-search
  // fallback. We assert on the short-vector method specifically so a regression to
  // "brute force disguised as an attack" fails CI.
  const trials = 60;
  let shortVectorWins = 0;
  let bruteWins = 0;
  let exact = 0;
  for (let i = 0; i < trials; i += 1) {
    const inst = await generateLWE(4, 71, 2);
    const atk = attackLWE(inst.A, inst.b, 71, 4, 2);
    const isExact = atk.recovered !== null && atk.recovered.join(',') === inst.secret.join(',');
    if (isExact) exact += 1;
    if (atk.recoverMethod === 'short-vector' && isExact) shortVectorWins += 1;
    if (atk.recoverMethod === 'bruteforce' && isExact) bruteWins += 1;
  }
  const successRate = exact / trials;
  const reductionRate = shortVectorWins / trials;
  // Overall success can dip on the occasional high-noise instance, so allow margin.
  assertTrue(successRate >= 0.85, `Expected toy attack success rate >= 0.85, got ${successRate}`);
  // The load-bearing invariant: the toy LWE break must come from LATTICE REDUCTION.
  // If recoveries ever start leaning on the exhaustive-search baseline, the demo is
  // back to teaching a facade -- fail CI loudly.
  assertTrue(
    exact > 0 && shortVectorWins / exact >= 0.95,
    `Expected >=95% of recoveries via short-vector reduction, got ${exact > 0 ? shortVectorWins / exact : 0} (short=${shortVectorWins}, brute=${bruteWins})`,
  );

  // Regression guard for the sigma-aware acceptance threshold: a genuine recovery
  // at higher noise must be REPORTED as success, not rejected by a fixed bound.
  let noisyExact = 0;
  const noisyTrials = 30;
  for (let i = 0; i < noisyTrials; i += 1) {
    const inst = await generateLWE(4, 257, 5);
    const atk = attackLWE(inst.A, inst.b, 257, 4, 5);
    if (atk.success && atk.recovered !== null && atk.recovered.join(',') === inst.secret.join(',')) {
      noisyExact += 1;
    }
  }
  const noisyRate = noisyExact / noisyTrials;
  assertTrue(
    noisyRate >= 0.7,
    `Expected sigma=5 recoveries to be reported as success >= 0.7, got ${noisyRate} -- threshold not noise-aware?`,
  );

  // Seeded determinism: the same seed must reproduce the exact same instance, so
  // shareable lab links and canonical labs are byte-for-byte repeatable. Restore
  // crypto randomness afterwards so the rest of the suite stays non-deterministic.
  setSeed(42);
  const seededA = await generateLWE(5, 97, 2);
  setSeed(42);
  const seededB = await generateLWE(5, 97, 2);
  assertTrue(
    JSON.stringify(seededA) === JSON.stringify(seededB),
    'Seeded RNG is not reproducible: same seed produced different instances',
  );
  setSeed(7);
  const seededC = await generateLWE(5, 97, 2);
  assertTrue(
    JSON.stringify(seededA) !== JSON.stringify(seededC),
    'Different seeds produced identical instances (seed is being ignored)',
  );
  setSeed(null);

  const hard = await generateLWE(20, 3329, 1);
  const hardAttack = attackLWE(hard.A, hard.b, 3329, 20);
  assertTrue(!hardAttack.success && hardAttack.recovered === null, 'Expected n=20 hard case to fail');

  console.log('verify:algo passed');
  console.log(`toy_success_rate=${successRate.toFixed(3)}`);
  console.log(`reduction_recovery_rate=${reductionRate.toFixed(3)}`);
  console.log(`noisy_sigma5_success_rate=${noisyRate.toFixed(3)}`);
}

await run();
