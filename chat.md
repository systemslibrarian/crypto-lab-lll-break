# Making This Demo the Gold Standard

This project already has the hard part many lattice demos miss: it teaches a real attack workflow instead of hand-waving. The Kannan embedding is explicit, toy LWE recovery is checked through a short-vector path, BKZ improvements are logged, and the README/guardrails avoid the false claim that LLL breaks real Kyber. To make it the gold standard for teaching LLL/BKZ/LWE attacks, I would focus on turning that correctness into a guided learning system.

## 1. Add a Clear Learning Contract

Give every exhibit a small, explicit learning objective and a success check.

Suggested objectives:

- Exhibit 1: Explain that many bases can generate the same lattice and determinant is invariant under unimodular changes.
- Exhibit 2: Compute Gram-Schmidt vectors and interpret `mu` coefficients.
- Exhibit 3: Predict whether LLL will size-reduce, swap, or advance from the Lovasz condition.
- Exhibit 4: Derive why `(-s, e, +-1)` is short in the primal embedding and recover `s` only when reduction exposes that vector.
- Exhibit 5: Explain why toy parameters fail while Kyber-like dimensions remain out of reach.

Gold-standard move: add a small “Can you answer this now?” prompt at the end of each exhibit. It does not need grading at first; even a revealable answer is enough to shift the app from demo to lesson.

## 2. Show the Invariants, Not Just the Results

The best teaching upgrade is to make every algorithmic step auditable.

For LLL, display these per step:

- `|mu[k][j]| <= 1/2` for size reduction, with pass/fail coloring.
- Lovasz left side vs right side, already partly shown in GSO, but integrated directly into the LLL trace.
- The exact row operation performed, such as `b_k <- b_k - c*b_j`.
- A running unimodular transform matrix `U` where `B_reduced = U * B_original`.
- Determinant before/after each operation.

This would make the core promise visible: LLL changes the basis without changing the lattice.

## 3. Add a Guided Mode and an Expert Mode

Right now the app is rich, but all learners see the same density. Split the experience into two modes.

Guided mode:

- One primary action per exhibit.
- Short explanations next to the exact number currently changing.
- “Try this next” presets that create predictable teaching moments.
- Revealable math derivations instead of always-visible dense output.

Expert mode:

- Full matrices, logs, BKZ tours, norm gaps, embedding basis, and raw vectors.
- Copy/export experiment output.
- Adjustable seeds and URL-shareable state.

This keeps beginners oriented without taking power away from instructors or advanced learners.

## 4. Make Experiments Reproducible

The app currently uses browser crypto randomness, which is good for avoiding fake randomness but hard for teaching repeatable examples. Add a seedable deterministic mode for classroom runs.

Recommended features:

- Optional seed input.
- “Copy lab link” that serializes parameters, seed, selected exhibit, matrix, and generated instance.
- “Copy transcript” that exports the LWE instance, reduction settings, result, and interpretation.
- A few named canonical labs: `LLL swap cascade`, `delta sensitivity`, `toy LWE success`, `BKZ helps`, `Kyber-like failure`.

Keep crypto randomness as the default for fresh exploration, but make deterministic mode available when a teacher needs every student to see the same run.

## 5. Upgrade the LWE Attack Explanation

Exhibit 4 is the signature teaching moment. It should get the most polish.

Add a step-by-step derivation panel:

1. Start with `b = A*s + e mod q`.
2. Build the embedding rows.
3. Combine the last row with `-s_i` times the secret rows.
4. Reduce the middle block modulo `q` using the `qI` rows.
5. Arrive at `(-s, e, 1)`.
6. Explain why its norm is small compared with typical lattice vectors.

Then pair the derivation with the actual recovered vector from the current run, highlighting:

- secret block
- error block
- embedding coordinate
- residual norm
- whether recovery came from short-vector extraction or the brute-force baseline

The app already reports much of this. The improvement is to make the derivation and live vector line up visually so learners can see the theory land.

## 6. Be More Precise About “Toy BKZ”

The BKZ implementation is intentionally educational. That is fine, but gold-standard teaching should label exactly what is modeled and what is simplified.

Add a “What this BKZ model includes / omits” section:

- Includes: LLL preprocessing, local block enumeration, block improvements, tours, convergence logs.
- Simplifies: projected block handling, pruning strategies, floating-point precision concerns, modern BKZ variants, high-performance enumeration, and real cost modeling.
- Consequence: it teaches the shape of BKZ reasoning, not production cryptanalysis.

This prevents advanced learners from leaving with a distorted idea of real lattice-reduction tooling.

## 7. Add Parameter Intuition Beyond a Single Threshold

Exhibit 5 should become a parameter intuition lab, not just a threshold chart.

Add side-by-side sliders for:

- dimension `n`
- modulus `q`
- noise `sigma`
- sample count `m`
- block size `beta`

Show derived quantities:

- embedding dimension
- estimated target norm `sqrt(||s||^2 + ||e||^2 + 1)`
- shortest vector found
- norm gap
- rough root-Hermite factor intuition
- plain-English status: “visible to LLL”, “needs stronger BKZ”, “outside this toy model”

Important: avoid presenting the estimate as a real Kyber security claim. Frame it as intuition and point serious readers to lattice estimators and production reduction libraries.

## 8. Add Challenges and Assessments

Gold-standard teaching tools let learners test themselves.

Suggested challenges:

- Predict whether a 2D basis will swap under LLL.
- Choose a delta that causes fewer/more swaps, then explain why.
- Find a bad basis that generates the same lattice as a good one.
- Given an embedding vector, identify the secret and error blocks.
- Tune toy LWE parameters until LLL fails, then explain which parameter mattered most.
- Compare BKZ-2 and BKZ-6 on the same seeded instance.

Add revealable answers and a short explanation for each. This turns the app into a lab exercise rather than a passive visualization.

## 9. Add Instructor Materials

If you want this to be adopted by teachers, include classroom-ready assets.

Recommended files:

- `docs/lab-30-min.md`: quick interactive walkthrough.
- `docs/lab-60-min.md`: full LLL plus toy LWE lab.
- `docs/lab-90-min.md`: deeper BKZ/parameter exploration.
- `docs/worksheet.md`: student questions.
- `docs/answer-key.md`: instructor answers.
- `docs/glossary.md`: lattice, basis, determinant, GSO, mu, Lovasz, SVP, CVP, LWE, embedding, BKZ, beta.

These assets are often the difference between “cool demo” and “course-ready module.”

## 10. Strengthen Visual Explanations

The canvases are useful. The next level is better annotation.

Visual upgrades:

- Highlight the exact vector being size-reduced.
- Animate row operations as vector subtraction.
- Show the fundamental parallelogram area in Exhibit 1.
- In the LLL chart, mark swaps as vertical ticks.
- In the LWE attack, show a compact block diagram of the embedding matrix.
- Use consistent colors for original basis, current basis, GSO vectors, target vector, and recovered vector across all exhibits.

The goal is not more decoration. It is visual continuity: once a color means “GSO vector,” it should mean that everywhere.

## 11. Add Better Failure Stories

A great cryptography lesson teaches failure clearly.

Add curated failed runs:

- LLL fails but brute-force fallback succeeds.
- BKZ improves the basis but does not expose an embedding-shaped vector.
- Noise too high for the current toy settings.
- Dimension too high for browser reduction.
- Kyber-like parameters are intentionally blocked from fake recovery.

For each failure, explain what failed: reduction quality, target norm gap, embedding shape, or computational budget. This will prevent learners from reading every failure as “the code broke.”

## 12. Separate “Attack Result” From “Teaching Baseline”

The brute-force fallback is useful pedagogically, but it should be visually and semantically separated from the lattice attack.

Suggested UI language:

- “Lattice attack result: success/failure.”
- “Teaching baseline: brute force found the secret / did not run / skipped because n is too large.”
- “Do not count baseline recovery as a lattice break.”

The code already labels the fallback. The gold-standard version should make that distinction impossible to miss.

## 13. Add References and Further Reading

Add a short, curated reading path instead of a huge bibliography.

Suggested tiers:

- First: intuitive lattice and LLL notes.
- Next: Regev-style LWE background.
- Then: Kannan embedding and primal attacks.
- Advanced: BKZ, root-Hermite factor, lattice estimators, and real reduction libraries.

Each reference should include one sentence saying why it is worth reading.

## 14. Improve Verification Coverage

The existing verification script is already unusually strong for a browser demo. I would extend it in three directions.

Algorithm checks:

- Assert determinant preservation for LLL row operations.
- Assert `B_reduced = U * B_original` once the transform matrix exists.
- Add seeded fixtures for known LLL traces.
- Add known embedding fixtures where `(-s, e, 1)` is constructed exactly.

UI checks:

- Add Playwright smoke tests for all buttons and sliders.
- Verify canvases render nonblank pixels.
- Verify the Kyber-like button never reports recovery.
- Verify moving sliders after generating an LWE instance does not silently mismatch `n/q/sigma`.

Accessibility checks:

- Run an automated accessibility pass.
- Confirm keyboard operation for all controls.
- Confirm canvas information has meaningful textual equivalents.

## 15. Polish the Project Surface

Small project-level changes would make the repository feel more authoritative.

Recommended additions:

- README table of exhibits and learning goals.
- Screenshot or animated GIF of the main lesson flow.
- “Accuracy boundaries” section: what is real, what is simplified, what not to infer.
- `docs/model-limitations.md` for LLL/BKZ/LWE caveats.
- Remove unused Vite starter residue such as `src/counter.ts` if it is not used.
- Add `npm run test:ui` once browser checks exist.

## Suggested Priority Order

1. Add exhibit learning objectives and revealable checks.
2. Add LLL invariant display: row operations, determinant, Lovasz, size-reduction pass/fail.
3. Add reproducible seeded labs and shareable URLs.
4. Upgrade Exhibit 4 with a live derivation-to-vector explanation.
5. Separate lattice attack result from brute-force teaching baseline in the UI.
6. Add curated challenge mode.
7. Add instructor docs and worksheets.
8. Add UI/accessibility verification.
9. Expand parameter intuition with better derived metrics.
10. Add references and model limitations.

## North Star

The gold-standard version should let a learner leave with these claims in their own words:

- A lattice is unchanged by unimodular basis operations, but the basis quality can change dramatically.
- Gram-Schmidt data tells LLL when to size-reduce, swap, or advance.
- LLL and BKZ search for short vectors; they do not magically decrypt.
- Toy LWE breaks when the embedding makes `(-s, e, 1)` unusually short and reduction is strong enough to expose it.
- Real post-quantum parameters are chosen so this attack shape is computationally out of reach.
- A browser demo can teach the mechanics, but it is not a replacement for production lattice estimators or reduction libraries.

That combination of visible math, honest limitations, reproducible experiments, and self-checking exercises is what would make this a reference-quality teaching module.