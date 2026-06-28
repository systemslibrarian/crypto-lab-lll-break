# Student Worksheet: LLL / BKZ / LWE

Name: ______________________________   Date: ______________

Work through the five exhibits in the live demo
(`https://systemslibrarian.github.io/crypto-lab-lll-break/` or `npm run dev`).
Every question requires using the app. Write your answers in the space provided.
Questions get harder as you go.

---

## Part A — What Is a Lattice? (Exhibit 1)

**Q1.** Set b1 length and b2 length to any values. Read the `det(Lambda)` label.
Now drag the b2 *angle* slider while keeping both lengths fixed. Does the
determinant change? In one sentence, say what geometric quantity the determinant
equals.

> ______________________________________________________________________
>
> ______________________________________________________________________

**Q2.** Click "Same lattice, different basis". The green arrows change but the
`det(Lambda)` label does not. Why must the determinant stay the same when the
basis changes this way?

> ______________________________________________________________________
>
> ______________________________________________________________________

---

## Part B — Gram-Schmidt (Exhibit 2)

**Q3.** With the default matrix (`3 1 / 1 2`), click "Next" until you reach Step 2
and record `mu21`. In your own words, what does `mu21` measure?

> mu21 = __________   Meaning: ________________________________________
>
> ______________________________________________________________________

**Q4.** Continue to Step 3 (the Lovasz check). Record the two numbers shown
(`0.75 * ||b1~||^2` and `||b2~ + mu21*b1~||^2`) and whether the condition holds.

> Left = __________   Right = __________   Holds? ______________________

**Q5.** Change the matrix to `1 9 / 10 0` and step through again. Does the Lovasz
condition still hold? What does a *failed* condition tell LLL to do next?

> ______________________________________________________________________
>
> ______________________________________________________________________

---

## Part C — LLL Step-by-Step (Exhibit 3)

**Q6.** Load the "Classic example" preset (`19 2 / 7 1`) and click "Step" repeatedly
to the end (or "Auto"). Record the final reduced basis vectors and the final
"Swap count".

> Reduced basis: ______________________________   Swaps: ______________

**Q7.** Watch the "Lattice det |det B|" line throughout the run. Compare it before
the first step and after the last step. What does its behavior prove about what LLL
does to the lattice?

> ______________________________________________________________________
>
> ______________________________________________________________________

**Q8.** Read the `Transform U` line on the final step and its `det(U)` value. What
relationship does `U` express between the reduced and original bases, and why does
`det(U)` have the value it does?

> U = ______________________________   det(U) = ______   Why: __________
>
> ______________________________________________________________________

**Q9.** Reset, then move the **delta** slider to its minimum (~0.5) and run to the
end; note the swap count. Set delta near its maximum (~0.999) and run again. Which
delta produced more swaps, and what does delta control?

> delta low swaps = ______   delta high swaps = ______
>
> Interpretation: ______________________________________________________

---

## Part D — Break a Toy LWE Instance (Exhibit 4)

**Q10.** Set `n = 4`, `q = 71`, `sigma = 2`, `beta = 2`. Click "Generate LWE
Instance" and record the printed secret `s` and the number of samples `m`. Confirm
the relationship between `m` and `n`.

> s = ______________________   m = ______   Relationship: m = n + ______

**Q11.** With the same instance, click "Run LLL/BKZ Attack". Find the
"Structure (-s, e, +-1) read off this vector" block. Does the recovered secret
block match the `s` printed in Q10? What was the embedding coordinate value?

> Recovered s = ______________________   embed coord = ______   Match? ____

**Q12.** Raise `sigma` to a large value (e.g. 9-10), regenerate, and run the attack
several times. Describe what happens to the "Lattice attack result" and the
norm-gap meter. Which parameter did you change, and why does it hurt the attack?

> ______________________________________________________________________
>
> ______________________________________________________________________

**Q13.** Find a run (high `sigma`, small `n <= 8`) where the output says
"Lattice attack result: FAILED" but "Teaching baseline: brute force found the
secret". Did you break LWE? Explain the difference between the two recovery paths.

> ______________________________________________________________________
>
> ______________________________________________________________________

**Q14.** Click "Try Kyber-512 parameters". Record the printed `n`, `q`, required
`beta`, and cost exponent. Does the app report recovering the secret? Why is that
the correct behavior?

> n = ______   q = ______   beta ~ ______   cost ~ 2^______   Recovered? ____
>
> ______________________________________________________________________

---

## Part E — Parameter Explorer (Exhibit 5)

**Q15.** Drag the `n` slider from small to large and watch the SECURE panel's
"Estimated beta" and "Attack cost". Roughly where does the status flip from "still
toy scale" to "practical attacks out of reach"? Using `cost ~ 2^(0.292*beta)`,
explain in one or two sentences why large real parameters resist LLL/BKZ.

> Flip near n = ______
>
> ______________________________________________________________________
>
> ______________________________________________________________________

---

### Stretch (optional)

**S1.** In Exhibit 4, pick one seeded-feeling failure at `beta = 2`, then raise
`beta` to 6-8 and rerun on a freshly generated instance of the same `n/q/sigma`.
Describe any change in tours, block improvements, or success. Why might a larger
`beta` help expose the target vector?

> ______________________________________________________________________
>
> ______________________________________________________________________
