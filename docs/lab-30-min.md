# Lab: 30-Minute Guided Walkthrough

**Audience:** students or engineers meeting LLL/BKZ/LWE for the first time.
**Goal:** a quick, correct tour of all five exhibits and the key takeaways.
**Setup:** open `https://systemslibrarian.github.io/crypto-lab-lll-break/` (or run
`npm install && npm run dev`). No accounts or keys needed.

Timings are guidance; the bracketed `[min]` marks are cumulative.

---

## 0. Frame the session (2 min) `[2]`

Say up front what this is and is not:

- It teaches the **shape** of lattice reduction and the primal LWE attack.
- It does **real** reduction and **real** short-vector recovery on tiny instances.
- It does **not** break real Kyber. The brute-force fallback is exhaustive search,
  not a lattice break — see `docs/model-limitations.md`.

## 1. Exhibit 1 — What Is a Lattice? (4 min) `[6]`

1. Move the b1/b2 angle and length sliders; note the dot grid (the lattice) and the
   `det(Lambda)` label.
2. Click **"Same lattice, different basis"** two or three times.

**Point to make:** the green basis vectors change, but `det(Lambda)` does not — the
button multiplied the basis by a unimodular (`det = +-1`) matrix, which is a new
basis for the *same* lattice. Determinant = area of the fundamental parallelogram,
an invariant.

## 2. Exhibit 2 — Gram-Schmidt (5 min) `[11]`

1. Keep the default matrix `3 1 / 1 2`. Click **"Next"** through all three steps.
2. Step 1 shows `b1*`; Step 2 shows `mu21` and `b2*`; Step 3 shows the Lovasz check
   with both sides computed.

**Points to make:** `mu21` is how much of `b1*` lives inside `b2`. LLL will later
subtract `round(mu)` times a vector (an *integer*, to stay on the lattice). The
Lovasz check is the decision rule LLL uses to swap or advance.

## 3. Exhibit 3 — LLL Step-by-Step (7 min) `[18]`

1. Load the **"Classic example"** preset (`19 2 / 7 1`).
2. Click **"Step"** a few times, then **"Auto"** to finish. Watch the canvas
   (red = original, green = current, dashed blue = GSO) and the log.
3. Read the metrics panel, especially:
   - `Lattice det |det B|: ... -> ...` (unchanged),
   - `Transform U` with `det(U) = +-1` (unimodular),
   - falling orthogonality defect and shortest basis vector.

**Point to make:** every step is a unimodular operation, so `reducedBasis = U *
original` with `det(U) = +-1`. LLL changes the *basis*, never the *lattice*. Quickly
nudge the **delta** slider higher and reset to show it drives more swaps / stronger
reduction.

## 4. Exhibit 4 — Break a Toy LWE Instance (8 min) `[26]`

1. Set `n = 4`, `q = 71`, `sigma = 2`, `beta = 2`. Click **"Generate LWE
   Instance"**. Read the printed secret `s`, samples `m = n + 4`, matrix `A`,
   vector `b`, and the embedding lattice (rows `[I | A^T | 0]`, `[0 | qI | 0]`,
   `[0 | b | 1]`).
2. Click **"Run LLL/BKZ Attack"**. Find:
   - the "Structure (-s, e, +-1)" block — the recovered secret matches `s`,
   - "Lattice attack result: SUCCESS" when reduction exposes the target vector.
3. Raise `sigma` to ~9-10, regenerate, and rerun. Show the attack now often FAILS,
   and if `n <= 8` the **"Teaching baseline: brute force"** may still find `s`.
   Stress: that is exhaustive search, **not** a lattice break.
4. Click **"Try Kyber-512 parameters"**: it prints `n=256, q=3329, beta~400,
   cost~2^117` and reports **no** recovery.

**Point to make:** the secret rides on the short `(-s, e, +-1)` vector; reduction
either exposes it or it doesn't. Brute force is a separate, honest baseline.

## 5. Exhibit 5 — Parameter Explorer (3 min) `[29]`

1. Drag the `n` slider and watch the SECURE panel's `Estimated beta` and
   `Attack cost`, plus the threshold chart (orange line near `n = 50`).

**Point to make:** `cost ~ 2^(0.292*beta)` and required `beta` grows with `n`, so
real parameters push attacks out of reach. This is intuition, not a certified
estimate.

## 6. Wrap-up (1 min) `[30]`

Have students state these in their own words:

- A lattice is unchanged by unimodular basis operations; basis *quality* changes.
- Gram-Schmidt data tells LLL when to size-reduce, swap, or advance.
- LLL/BKZ search for short vectors; they do not magically decrypt.
- Toy LWE breaks when `(-s, e, +-1)` is short *and* reduction is strong enough.
- Real parameters keep this attack shape computationally out of reach.
- This is a teaching model, not a replacement for real estimators/libraries.
