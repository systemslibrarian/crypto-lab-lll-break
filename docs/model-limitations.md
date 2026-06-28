# Model Limitations: What This Toy Models vs. Omits

**Read this before teaching with the demo.** crypto-lab-lll-break is built to teach
the *shape* of LLL/BKZ reduction and the primal LWE attack. It performs **genuine**
lattice reduction and **genuine** short-vector recovery on small instances, but it
is deliberately simplified and is **not** production cryptanalysis. Nothing in this
demo shows that LLL or BKZ breaks real Kyber.

The single most important honesty point: the brute-force fallback in Exhibit 4 is
an **exhaustive search over a tiny secret space**, not a lattice break. It must
never be counted as one. The UI labels it "Teaching baseline" precisely so this
distinction is impossible to miss.

---

## LLL (`src/lll.ts`, `src/lattice.ts`)

**What it models faithfully:**

- Real size reduction: `b_k <- b_k - round(mu[k][j]) * b_j` until `|mu[k][j]| <= 1/2`.
- The real Lovasz condition: `delta * ||b*_{k-1}||^2 <= ||b*_k + mu * b*_{k-1}||^2`,
  with swap-on-failure / advance-on-success, default `delta = 0.75`.
- A complete, inspectable step trace (size-reduce / swap / advance / done) with
  before/after basis, GSO, and `mu` at every step.
- A running **unimodular transform `U`** with `reducedBasis = U * original` and
  `det(U) = +-1`, plus live determinant before/after, making the "same lattice"
  invariant auditable rather than asserted.
- Orthogonality-defect tracking so basis-quality improvement is visible.

**What it omits or simplifies:**

- Floating-point engineering. GSO and `mu` are plain JavaScript doubles; there is
  no use of the careful FP / exact-arithmetic strategies (e.g. `fplll`'s
  proved/heuristic variants) that real implementations need for stability at
  scale. A hard operation cap (`maxOps = 20000`) guards against pathological runs.
- Scale. It is meant for 2-3 dimensional teaching bases, not the hundreds of
  dimensions of real instances.
- Deep LLL, segment/recursive LLL, and other performance variants.

---

## BKZ (`src/bkz.ts`)

**What it models faithfully:**

- **LLL preprocessing** before and between tours.
- **Exact SVP on each block** by full enumeration (`enumerateSVP`), restricted to
  blocks of size up to 8 (the demo's `beta` slider goes 2-8).
- **Lattice-preserving insert-and-reduce**: the enumerated short vector is added
  as an extra generator and the `(block+1)` set is LLL-reduced, dropping the
  resulting zero vector. This preserves the block lattice exactly (unlike the
  naive "overwrite `b_0`" shortcut, which is only valid when the short vector's
  coefficient on `b_0` is `+-1`).
- **Tours** across all block positions, **convergence detection** (stop when a
  tour makes no improvement), and per-tour / per-block improvement logging.

**What it omits or simplifies:**

- **Projected-block handling.** Real BKZ enumerates in the lattice *projected*
  orthogonally to the preceding basis vectors; this toy enumerates each block as a
  standalone sublattice. This is the most significant theoretical simplification.
- **Pruning.** No (extreme) pruning of the enumeration tree; the demo does plain
  exhaustive enumeration, which is why blocks must stay tiny.
- **Floating-point / precision engineering** for enumeration radii and GSO.
- **Progressive BKZ and modern variants** (BKZ 2.0, self-dual BKZ, G6K-style
  sieving, etc.).
- **Real cost models.** Tours are capped at 20; there is no realistic
  time/memory accounting. Exhibit 5's `cost ~ 2^(0.292*beta)` is rule-of-thumb
  intuition, not a calibrated estimate.

**Consequence:** the demo teaches *how BKZ reasons* (reduce, enumerate a block,
reinsert, repeat, converge), not how to mount a real attack.

---

## The LWE Attack (Exhibit 4, `src/bkz.ts`)

**What it models faithfully:**

- The **textbook primal (Kannan) embedding**, rows `[I | A^T | 0]`,
  `[0 | qI | 0]`, `[0 | b | 1]`, so the short vector is `(-s, e, +-1)`
  (see the long comment on `buildLWELattice`).
- A **genuine reduction attack**: it builds the embedding, runs LLL or toy BKZ
  on it, and reads `s` off the embedding-shaped reduced vector
  (`recoverFromShortVector`: last coordinate `+-1`, secret is the first `n`
  coordinates negated by the embedding sign).
- Honest instance generation (`generateLWE`): secret entries drawn uniformly from
  `{0..4}`, `m = n + 4` samples, discrete Gaussian errors with standard deviation
  `sigma`.
- A norm-gap meter comparing the shortest vector found against the estimated
  target norm `sqrt(6n + m*sigma^2 + 1)` (the `6n` reflects `E[s_i^2] = 6` for a
  uniform `{0..4}` secret).
- Strict separation of "Lattice attack result" from "Teaching baseline" in the
  output and the meter coloring.

**What it omits or simplifies — and what the fallback is NOT:**

- **The brute-force fallback is exhaustive search, not a lattice break.**
  `bruteForceRecover` enumerates every candidate in `{0..4}^n` (only when
  `n <= 8`) and accepts one whose residual looks like genuine LWE noise. It only
  works because the toy secret space is microscopic. The UI explicitly prints
  "do NOT count baseline recovery as a lattice break."
- **Small dimensions only.** `attackLWE` refuses `n >= 20`; enumeration and brute
  force are gated at block/secret size 8. Real LWE lives far beyond this.
- **No real cost or success-probability modeling**, no parameter-driven failure
  analysis beyond the visible norm gap.
- **The Kyber-512 button is illustrative, not computational.** It does not run a
  reduction; it prints the *regime* (`n=256`, `q=3329`, embedding dim ~513,
  required `beta ~ 400`, cost `~2^117`) to show the attack is out of reach. It is
  hard-wired never to report recovery.

**Do not claim:** that LLL or this toy BKZ breaks real Kyber/ML-KEM, or any
real-world parameter set. The demo's own README and guardrails say the same.

---

## Further Reading — A Tiered Path

Read top to bottom; each tier assumes the previous one. One sentence per item says
why it is worth reading.

**Tier 1 — Intuitive lattices and LLL**

- Regev's lecture notes / many "Intro to Lattices" course pages — the gentlest way
  to build the geometric picture (bases, determinant, fundamental cell) this demo
  draws.
- Lenstra, Lenstra & Lovasz, *Factoring polynomials with rational coefficients*
  (1982) — the original LLL paper; worth skimming to see size reduction and the
  Lovasz condition exactly as Exhibits 2-3 implement them.
- Galbraith, *Mathematics of Public Key Cryptography*, lattice chapters — a clear,
  modern textbook treatment of LLL and reduced bases with worked examples.

**Tier 2 — Regev-style LWE background**

- Regev, *On lattices, learning with errors, random linear codes, and
  cryptography* (2005/2009) — the paper that defines LWE and its hardness; read it
  to understand what `b = A*s + e (mod q)` is and why it is hard.
- Peikert, *A Decade of Lattice Cryptography* — a survey tying LWE to real schemes;
  good for seeing where ML-KEM/Kyber fits.

**Tier 3 — Kannan embedding and primal attacks**

- Kannan, *Minkowski's convex body theorem and integer programming* (1987) — the
  source of the embedding technique Exhibit 4 uses to turn LWE into an SVP
  instance.
- Albrecht & Ducas, *Lattice Attacks on NTRU and LWE: A History of Refinements*
  — explains the primal (uSVP) and dual attacks and why `(-s, e, 1)` being short is
  the whole game.

**Tier 4 — BKZ, root-Hermite factor, estimators, and libraries**

- Chen & Nguyen, *BKZ 2.0: Better Lattice Security Estimates* — the reference for
  how real BKZ (pruning, tours, cost) behaves, i.e. exactly what this toy omits.
- Albrecht, Player & Scott, *On the concrete hardness of Learning with Errors*
  (and the maintained `lattice-estimator`) — the standard way to estimate the
  `beta` and cost to break given LWE parameters; use this, not Exhibit 5, for any
  security claim.
- **fplll** — a production floating-point LLL/BKZ library; shows what real
  reduction (with proper FP and pruning) actually requires.
- **G6K (General Sieve Kernel)** — state-of-the-art lattice sieving used in record
  SVP/LWE challenge solves; illustrates how modern attacks really scale.

This demo is a teaching scaffold that points at these tools; it does not replace
them.
