# Instructor Answer Key

Answers and explanations for `docs/worksheet.md`. Exact numbers vary because the
demo uses real browser randomness for the LWE instance (Exhibit 4) and the
"Random bad 2D" preset; the *behavior* and reasoning are what to assess.

---

## Part A — What Is a Lattice?

**A1.** The determinant **does** change when the b2 angle changes:
`det = |b1| * |b2| * sin(theta)`, so rotating b2 changes the angle `theta` between
the vectors, which changes the parallelogram's area and therefore the determinant.
The determinant equals the **area of the fundamental parallelogram**. *(Teaching
point: the sliders change the actual geometry, so the determinant should change.
The determinant is invariant only under unimodular basis changes — which is exactly
what Q2 isolates with the "Same lattice, different basis" button.)*

**A2.** Clicking the button multiplies the basis by an integer matrix with
determinant `+-1` (a unimodular matrix). That relabels exactly the same set of
lattice points, so the fundamental-cell area — the absolute determinant — is
preserved even though the vectors look different. *(This is the core invariant of
the whole demo.)*

## Part B — Gram-Schmidt

**A3.** For `3 1 / 1 2`: `b1* = b1 = [3, 1]`, and
`mu21 = <b2, b1*> / <b1*, b1*> = (1*3 + 2*1) / (9 + 1) = 5/10 = 0.5`.
`mu21` measures how much of the orthogonal direction `b1*` is contained in `b2`
(the projection coefficient).

**A4.** `b1* = [3,1]`, `||b1*||^2 = 10`, so left `= 0.75 * 10 = 7.5`.
`b2* = b2 - mu21*b1* = [1,2] - 0.5*[3,1] = [-0.5, 1.5]`;
`b2* + mu21*b1* = b2 = [1,2]`, so right `= 1 + 4 = 5`.
Left `7.5 <= 5`? **No** — the condition fails for this basis. *(Note: the demo's
Step-3 right-hand side reconstructs `||b2* + mu21*b1*||^2 = ||b2||^2`.)*

**A5.** Entering `1 9 / 10 0` is a deliberately skewed basis; the Lovasz condition
fails. A failed condition tells LLL to **swap** `b_{k-1}` and `b_k` (then continue
reducing), because the later Gram-Schmidt vector is shorter than `delta` times the
earlier one and the order should be improved.

## Part C — LLL Step-by-Step

**A6.** For `19 2 / 7 1`, LLL reduces to a much shorter basis (e.g. vectors like
`[-2, 1]` / `[5, -1]` up to sign/order; `det = |19*1 - 2*7| = 5` is preserved).
The swap count is small (typically 1-2). Accept any correct short reduced basis
with `|det| = 5`.

**A7.** `|det B|` is identical before the first step and after the last step (e.g.
`5.000 -> 5.000`). This proves LLL **changes the basis but not the lattice** — it
improves basis quality (shorter, more orthogonal vectors) while the lattice, and
hence its determinant, is unchanged.

**A8.** `U` is the accumulated integer transform with `reducedBasis = U * original`.
`det(U) = +-1` (unimodular) because every LLL step is a unimodular row operation
(integer size reduction `b_k -= c*b_j`, which has det 1, or a row swap, which has
det -1); the product of such operations is still `+-1`. This is the live proof of
the same-lattice invariant.

**A9.** A **larger delta produces more swaps** (and a stronger reduction); a small
delta (near 0.5) tolerates a worse-ordered basis and swaps less. delta controls how
aggressively LLL insists each Gram-Schmidt step not shrink too much — it is the
quality/work knob of the reduction. *(Default 0.75; range in the app 0.5-0.999.)*

## Part D — Break a Toy LWE Instance

**A10.** Accept the printed secret (entries in `{0..4}`). For `n = 4`, `m = 8`, so
`m = n + 4`. *(`generateLWE` always uses `m = n + 4` samples.)*

**A11.** On a successful lattice run the recovered secret block matches the printed
`s` exactly ("EXACT MATCH"), and the embedding coordinate is `+-1`. The secret is
read off the first `n` coordinates of the reduced `(-s, e, +-1)` vector, negated by
the sign of the embedding coordinate. *(If it instead failed and only brute force
recovered s, that is Q13's situation, not a lattice break.)*

**A12.** Raising `sigma` enlarges the error vector `e`, so the target
`(-s, e, +-1)` gets longer and the **norm gap** widens; the "Lattice attack result"
flips to FAILED more often and the meter turns warn/bad. The parameter changed is
the **noise standard deviation**; bigger noise destroys the shortness that the
attack depends on. *(With `n <= 8` the brute-force baseline may still find s — see
Q13.)*

**A13.** **No, you did not break LWE.** The lattice attack (short-vector
extraction) failed; only the **brute-force teaching baseline** succeeded, and that
is exhaustive search over the tiny `{0..4}^n` secret space (`n <= 8`), not
reduction. The app explicitly says "do NOT count baseline recovery as a lattice
break." A genuine break requires reduction to expose the `(-s, e, +-1)` vector.

**A14.** `n = 256`, `q = 3329`, required `beta ~ 400`, cost `~2^117`. The app does
**not** report recovery — it only prints the regime. That is correct: real Kyber
parameters put the required block size, and therefore the cost, far out of reach,
so no recovery should (or does) occur. The button is hard-wired never to claim a
break.

## Part E — Parameter Explorer

**A15.** The status flips around the marked **security threshold near n = 50** (the
chart draws an orange line there). With `cost ~ 2^(0.292*beta)` and the required
`beta` growing with dimension, large real parameters drive `beta` — and thus the
exponent — high enough that the attack, while still well-defined, is
computationally infeasible. *(Remind students this is intuition, not a certified
estimate; real numbers come from an LWE estimator.)*

### Stretch

**S1.** Raising `beta` (2 -> 6-8) makes BKZ enumerate larger blocks and can produce
more block improvements and a shorter head vector, sometimes flipping a borderline
instance to SUCCESS. A larger `beta` yields a stronger reduction that is more
likely to expose the short `(-s, e, +-1)` target — at the cost of more work. Note
the instance is freshly randomized, so this is a tendency, not a guarantee; the
deeper point is that block size is the strength knob, mirroring real BKZ.
