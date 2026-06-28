# Glossary

Concise reference for the terms used in **crypto-lab-lll-break**. Definitions are
kept to one or two sentences and matched to how each term is actually used in the
demo's five exhibits.

- **Lattice** — The set of all integer linear combinations of a set of basis
  vectors; in Exhibit 1 it is the grid of dots you see, and it does not change
  when you pick a different basis for it.

- **Basis** — A set of linearly independent vectors whose integer combinations
  generate the lattice; many different bases generate the *same* lattice, which
  is the whole point of Exhibit 1.

- **Determinant** — For a square basis matrix, the (signed) volume of one
  fundamental cell; its absolute value is an invariant of the lattice, shown in
  Exhibit 1 as `det(Lambda)` and in Exhibit 3 as `|det B|`, unchanged by
  reduction.

- **Fundamental parallelogram** — In 2D, the region spanned by the two basis
  vectors; tiling the plane with copies of it covers every point exactly once,
  and its area equals the absolute determinant.

- **Unimodular matrix** — An integer square matrix with determinant `+-1`;
  multiplying a basis by one produces a different basis for the *same* lattice,
  which is exactly what the "Same lattice, different basis" button and the
  running transform `U` in Exhibit 3 demonstrate.

- **Gram-Schmidt orthogonalization (GSO)** — A procedure that turns a basis
  `b_1..b_n` into mutually orthogonal vectors `b*_1..b*_n` by subtracting the
  projections onto earlier vectors; Exhibit 2 computes and draws these `b*`.

- **mu coefficients** — The Gram-Schmidt projection coefficients
  `mu[i][j] = <b_i, b*_j> / <b*_j, b*_j>`, measuring how much of an earlier
  orthogonal direction sits inside a later basis vector; Exhibit 2 reports
  `mu21` directly.

- **Size reduction** — Replacing `b_k` with `b_k - round(mu[k][j]) * b_j` so that
  every `|mu[k][j]| <= 1/2`; it shrinks a vector using an *integer* multiple,
  keeping it a lattice vector (the rounding is why an integer, not the exact real
  multiple, is used).

- **Lovasz condition** — The test `delta * ||b*_{k-1}||^2 <= ||b*_k + mu[k][k-1] * b*_{k-1}||^2`;
  if it holds, LLL advances, otherwise it swaps `b_{k-1}` and `b_k`. Exhibit 2
  evaluates both sides and Exhibit 3 acts on the result.

- **delta** — The Lovasz parameter in `(1/4, 1)` controlling how aggressively LLL
  reorders vectors; larger `delta` means a stronger (but slower) reduction. The
  demo defaults to `delta = 0.75` and exposes a slider in Exhibit 3.

- **LLL (Lenstra-Lenstra-Lovasz)** — A polynomial-time basis-reduction algorithm
  that interleaves size reduction and Lovasz-condition swaps to produce a
  reasonably short, nearly orthogonal basis; Exhibit 3 runs it step by step.

- **Orthogonality defect** — The ratio `(prod ||b_i||) / |det B|`, equal to 1 for a
  perfectly orthogonal basis and larger for a skewed one; Exhibit 3 charts it
  falling as LLL improves the basis.

- **SVP (Shortest Vector Problem)** — Find the shortest non-zero vector in a
  lattice; the toy BKZ in this demo solves SVP *exactly* on small blocks by
  enumeration (`enumerateSVP`, blocks of size up to 8).

- **CVP (Closest Vector Problem)** — Find the lattice point closest to a given
  target point; the primal/Kannan embedding turns an LWE (approximate-CVP)
  instance into an SVP instance so that reduction can attack it.

- **LWE (Learning With Errors)** — The hardness assumption `b = A*s + e (mod q)`
  with secret `s` and small noise `e`; recovering `s` from `(A, b)` is hard at
  large parameters and underpins ML-KEM/Kyber, Dilithium, and others.

- **Primal / Kannan embedding** — A construction that places an LWE instance into
  a lattice (rows `[I | A^T | 0]`, `[0 | qI | 0]`, `[0 | b | 1]`) so that the
  vector `(-s, e, +-1)` is unusually short; reading its first block recovers `s`.
  This is the core of Exhibit 4.

- **q-ary lattice** — A lattice that contains `q * Z^n` (defined by relations mod
  `q`); the `qI` block in the LWE embedding is what makes it q-ary, and the
  "LWE-style q-ary" preset in Exhibit 3 shows a small one.

- **BKZ (Block Korkine-Zolotarev)** — A stronger reduction that repeatedly solves
  SVP on local blocks of size `beta` and reinserts the short vector; this demo's
  toy BKZ does LLL preprocessing, exact block enumeration, lattice-preserving
  insert-and-reduce, and tours.

- **Block size beta** — The local block dimension BKZ works on; `beta = 2` is
  essentially LLL, and larger `beta` gives stronger reduction at exponentially
  higher cost. Exhibit 4 exposes a `beta` slider (2-8).

- **Tour** — One full sweep of BKZ across all block positions; the demo logs how
  many improvements each tour makes and stops when a tour produces none
  (convergence), capped at 20 tours.

- **Root-Hermite factor** — The per-dimension quality measure `delta_0` of a
  reduced basis (roughly `(||b_1|| / det^{1/n})^{1/n}`); smaller is better, and
  it is the standard way to summarize how strong a reduction (and thus how large
  a `beta`) you need.

- **Lattice estimator** — A tool (e.g. the Albrecht et al. LWE estimator) that
  predicts the `beta` and cost needed to break given LWE parameters; this demo is
  *not* one, and Exhibit 5 only gives intuition (`cost ~ 2^(0.292*beta)`), not a
  certified estimate.
