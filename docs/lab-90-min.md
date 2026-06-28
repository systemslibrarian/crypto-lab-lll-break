# Lab: 90-Minute Deep Dive — BKZ Block Size and Parameters

**Audience:** students who have completed the 60-minute lab (or equivalent) and are
ready for BKZ and parameter exploration.
**Goal:** understand how BKZ block size `beta` strengthens reduction, run
controlled `beta` experiments on the toy LWE attack, and build honest parameter
intuition (Exhibit 5) — with explicit attention to what the toy model omits.
**Setup:** `https://systemslibrarian.github.io/crypto-lab-lll-break/` or
`npm run dev`. Read alongside `docs/model-limitations.md`; finish `docs/worksheet.md`
including Part E and the Stretch question.

Cumulative timings in brackets.

---

## Part 0 — Recap and contract (5 min) `[5]`

Quickly restate the same-lattice invariant and the primal embedding from the
60-minute lab. Re-emphasize the honesty boundaries: the toy BKZ teaches the *shape*
of BKZ; the brute-force fallback is exhaustive search, never a lattice break; the
Kyber button is illustrative, not computational.

## Part 1 — From LLL to BKZ: what a block adds (15 min) `[20]`

Explain, then demonstrate, what this demo's BKZ actually does (`src/bkz.ts`):

1. **LLL preprocessing** runs first (and between tours).
2. For each block position, it takes a window of `beta` vectors and solves **exact
   SVP** on that block by full enumeration (`enumerateSVP`, blocks up to size 8).
3. It does **lattice-preserving insert-and-reduce**: the enumerated short vector is
   added as an extra generator and the set is LLL-reduced; the resulting zero
   vector is dropped, leaving a true basis of the same block lattice with the short
   vector pulled to the front.
4. One sweep across all positions is a **tour**; tours repeat until one makes no
   improvement (convergence), capped at 20.

Tie `beta = 2` to "essentially LLL" and larger `beta` to "stronger, costlier
reduction". Have students keep `docs/model-limitations.md` open: note the omissions
(projected blocks, pruning, FP engineering, modern BKZ, real cost models).

## Part 2 — Controlled beta experiment (25 min) `[45]`

Goal: see `beta` change the *outcome*, not just the runtime.

1. In Exhibit 4 set a regime that is borderline at `beta = 2`. A good starting
   point: `n = 6`, `q = 71`, `sigma = 4-6`. (Larger `sigma` lengthens `e` and
   widens the norm gap, making the attack harder so `beta` has room to matter.)
2. With `beta = 2`, **Generate** then **Run** the attack. Record:
   - Tours and improvements,
   - Tour log and any block improvements (`[start..end] ||b0|| X -> Y`),
   - Shortest vector, target-norm estimate, and the `gap` multiplier,
   - Lattice attack result and the norm-gap meter color.
3. **Without regenerating**, raise `beta` to 6, rerun **Run LLL/BKZ Attack** on the
   same instance (the demo reduces and recovers from the *same* basis, so `beta` is
   the only change). Record the same fields.
4. Repeat with `beta = 8`.

**Checkpoint A:** *What changed between `beta = 2`, 6, and 8 on the same instance?*
Expected: larger `beta` tends to make more / deeper block improvements and a shorter
head vector, sometimes flipping a borderline run from FAILED to SUCCESS. Bigger
blocks mean a stronger reduction more likely to expose `(-s, e, +-1)`.

**Checkpoint B:** *Why does the demo recover from the exact basis it reduced, rather
than re-running plain LLL?* Expected: otherwise the `beta` slider would be cosmetic —
a hidden LLL re-run would decide success regardless of the chosen block size.

> Note on randomness: each **Generate** makes a fresh instance, so compare `beta`
> values *on the same generated instance* (regenerate only when you want a new
> case). Run several instances to see tendencies, not a single anecdote.

## Part 3 — Mapping the success/failure frontier (20 min) `[65]`

Sweep parameters to find where the toy attack lives. For each setting, generate and
run 2-3 instances at `beta = 2` and again at `beta = 6`, recording success rate and
the typical norm gap:

| n | q | sigma | beta | typical result | norm gap |
|---|---|-------|------|----------------|----------|
| 4 | 71 | 2 | 2 | | |
| 4 | 71 | 8 | 2 | | |
| 4 | 71 | 8 | 6 | | |
| 8 | 71 | 4 | 2 | | |
| 8 | 71 | 4 | 8 | | |
| 12 | 257 | 2 | 2 | | |
| 12 | 257 | 2 | 8 | | |

Discuss the levers:

- **sigma up** -> longer `e` -> larger target norm -> wider gap -> harder.
- **q up** (relative to noise) -> generally easier to separate the short vector, but
  watch the interplay with `n`.
- **n up** -> larger embedding dimension (`n + m + 1`, with `m = n + 4`) -> harder,
  and recall `attackLWE` refuses `n >= 20` and enumeration/brute force are gated at
  size 8.
- **beta up** -> stronger reduction -> can rescue borderline cases.

**Checkpoint C:** *Which single parameter most reliably moved runs from SUCCESS to
FAILED in your table, and why?* Expected (typically): `sigma`, because it directly
lengthens the target vector the attack must isolate.

## Part 4 — Parameter intuition and the wall (15 min) `[80]`

1. Exhibit 5: drag `n`. Read the SECURE panel (`SVP approx factor 2^(n/2)`,
   `Estimated beta`, `Attack cost ~ 2^(0.292*beta)`) and the threshold chart (orange
   line near `n = 50`, gold line = current `n`).
2. Back in Exhibit 4, click **"Try Kyber-512 parameters"**: `n=256, q=3329,
   beta~400, cost~2^117`, and **no** recovery is reported.

**Checkpoint D:** *Using `cost ~ 2^(0.292*beta)` and the fact that required `beta`
grows with dimension, explain in two sentences why real parameters resist LLL/BKZ.*
Expected: the block size needed to surface the short vector grows with `n`, and cost
is exponential in `beta`, so at Kyber-like sizes the exponent is astronomically
large; the attack shape still exists but is computationally out of reach.

**Honesty checkpoint:** *Name two things this toy BKZ omits that a real attack
needs.* Expected (any two): projected-block enumeration, pruning, floating-point
engineering, progressive/modern BKZ, realistic cost models. Therefore Exhibit 5 is
intuition, not a certified estimate — real numbers come from an LWE estimator.

## Wrap-up and bridge to real tools (10 min) `[90]`

Have students articulate:

- `beta` is the strength knob of BKZ; larger blocks = stronger, costlier reduction
  that can expose a target vector LLL misses.
- The toy faithfully shows LLL preprocessing, exact block enumeration,
  lattice-preserving insertion, tours, and convergence — but omits the engineering
  that real attacks require.
- Security comes from the cost wall (`~2^(0.292*beta)`) at parameters where the
  required `beta` is huge.

Point serious follow-up at the real tooling (see `docs/model-limitations.md` and the
references in the README): the **LWE estimator** (Albrecht et al.) for concrete
hardness, **fplll** for production LLL/BKZ, and **G6K** for state-of-the-art
sieving. None of this demo's output is a security claim about real schemes.
