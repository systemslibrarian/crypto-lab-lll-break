# Lab: 60-Minute LLL + Toy LWE Attack

**Audience:** students who have seen the 30-minute tour, or a class with one full
period.
**Goal:** work LLL thoroughly (Exhibits 1-3) and run the toy primal LWE attack
end-to-end (Exhibit 4), with checkpoints you can grade or discuss.
**Setup:** `https://systemslibrarian.github.io/crypto-lab-lll-break/` or
`npm run dev`. Pair this with `docs/worksheet.md` (Parts A-D).

Cumulative timings in brackets.

---

## Part 0 — Orientation (3 min) `[3]`

State the contract (see `docs/model-limitations.md`): real reduction and real
short-vector recovery on tiny instances; **not** a break of real Kyber; the
brute-force fallback is exhaustive search, never counted as a lattice break.

## Part 1 — Lattices and invariance (8 min) `[11]`

1. Exhibit 1: move sliders, read `det(Lambda)`. Establish that the dot grid is the
   lattice and the determinant is the fundamental-cell area.
2. Click **"Same lattice, different basis"** several times.

**Checkpoint 1:** *Why does the determinant stay fixed when the basis changes?*
Expected: the new basis = old basis times a unimodular (`det = +-1`) integer
matrix, so it is a different basis for the same lattice; area is preserved.

## Part 2 — Gram-Schmidt and the Lovasz test (10 min) `[21]`

1. Exhibit 2, default `3 1 / 1 2`. Step through all three panels.
2. Record `mu21` and the Lovasz left/right values.
3. Change the matrix to a skewed one, e.g. `1 9 / 10 0`, and step again.

**Checkpoint 2:** *What does `mu21` measure, and what does a failed Lovasz check
make LLL do?* Expected: `mu21 = <b2, b1*>/<b1*, b1*>` is the projection of `b2` onto
`b1*`; a failed check triggers a **swap** of `b_{k-1}` and `b_k`.

## Part 3 — LLL step by step (18 min) `[39]`

1. Exhibit 3, **"Classic example"** (`19 2 / 7 1`). Use **"Step"** so students see
   each size-reduce / swap / advance entry in the log; finish with **"Auto"**.
2. Track in the metrics panel:
   - `Lattice det |det B|` (unchanged),
   - `Transform U` and `det(U) = +-1`,
   - orthogonality defect (and its chart) falling,
   - shortest basis vector and swap count.
3. **delta experiment:** Reset; set delta near 0.5, run, note swaps. Set delta near
   0.999, run, note swaps.
4. **3D:** load **"3D challenge"** (`19 2 3 / 7 1 6 / 2 9 5`) and run to show LLL
   beyond 2D. Optionally load **"LWE-style q-ary"** (`71 0 0 / 17 1 0 / 63 0 1`) to
   preview the embedding structure used in Part 4.

**Checkpoint 3:** *The numbers in the basis change completely but `|det B|` never
does. What does that prove?* Expected: every step is unimodular, so
`reducedBasis = U * original` with `det(U) = +-1`; the basis is rewritten but the
lattice is identical. LLL improves quality, not the lattice.

**Checkpoint 4:** *Which delta gave more swaps, and what does delta control?*
Expected: larger delta -> more swaps / stronger reduction; delta is the
quality/work knob.

## Part 4 — Toy LWE attack end to end (18 min) `[57]`

1. Exhibit 4: set `n = 4`, `q = 71`, `sigma = 2`, `beta = 2`.
2. **Generate LWE Instance.** Walk through the printout: secret `s`, samples
   `m = n + 4 = 8`, matrix `A`, vector `b = A*s + e (mod q)`, hidden errors `e`, and
   the embedding lattice rows `[I | A^T | 0]`, `[0 | qI | 0]`, `[0 | b | 1]`.
3. Derive on the board (or read the `buildLWELattice` comment): take the last row
   once and subtract `s_i` times row `i`; modulo `q` (via the `qI` rows) this yields
   `(-s, e, 1)`, a short vector.
4. **Run LLL/BKZ Attack.** Show the "Structure (-s, e, +-1)" block, the matched
   recovered secret, and "Lattice attack result: SUCCESS".
5. **Make it fail honestly:** raise `sigma` to ~9-10, regenerate, rerun a few times.
   The lattice attack flips to FAILED as the norm gap widens. Because `n = 4 <= 8`,
   the **brute-force teaching baseline** may still print the secret.

**Checkpoint 5:** *On a run where the lattice attack FAILED but brute force found
`s`, did we break LWE?* Expected: **No.** Brute force enumerated the tiny
`{0..4}^n` secret space — exhaustive search, not reduction. Only short-vector
extraction counts as a lattice break.

**Checkpoint 6 (interpretation):** *Which single parameter did raising `sigma`
change, and why did it hurt the attack?* Expected: it enlarged the error `e`, so the
target `(-s, e, +-1)` is no longer unusually short and reduction can't isolate it.

6. Click **"Try Kyber-512 parameters"**: `n=256, q=3329, beta~400, cost~2^117`, no
   recovery. Use this to bridge to the parameter intuition.

## Wrap-up (3 min) `[60]`

Round-robin: each student states one takeaway.

- Unimodular operations preserve the lattice; quality is what LLL improves.
- Gram-Schmidt `mu` data drives the size-reduce / swap / advance decision.
- The primal embedding makes `(-s, e, +-1)` short; reduction recovers `s` only when
  it exposes that vector.
- Brute force is a labeled baseline, never a lattice break.
- Real parameters keep the cost (`~2^(0.292*beta)`) out of reach — see the
  90-minute lab for BKZ block-size experiments and `docs/model-limitations.md` for
  the boundaries of this toy.
