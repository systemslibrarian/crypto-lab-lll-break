# crypto-lab-lll-break

## What It Is

This demo shows LLL lattice reduction and toy BKZ reduction on Learning With Errors (LWE) embeddings in the browser. It focuses on the attacker workflow: reduce a basis, search for short vectors, and test whether secret recovery is possible at a chosen parameter set. The core problem illustrated is how SVP-approximation quality affects practical LWE attacks. This is a post-quantum cryptography educational model, not a production cryptanalytic tool and not evidence that LLL breaks real Kyber parameters.

## When to Use It

- Teaching why post-quantum parameters matter in lattice cryptography. This demo makes the relationship between dimension, modulus, noise, and attack feasibility visible.
- Explaining LLL and BKZ mechanics to students or engineers. Step traces and block-improvement logs show what each reduction phase actually changes.
- Comparing toy insecure settings against Kyber-like settings. The same pipeline can be run at small and large parameters to show where attacks stop being effective.
- Demonstrating how primal LWE embeddings are constructed. The matrix and reduction output are shown directly so learners can inspect the attack surface.
- Not for real-world security assessment. This browser implementation is intentionally simplified and does not replace specialized lattice estimators or high-performance reduction libraries.

## Live Demo

https://systemslibrarian.github.io/crypto-lab-lll-break/

The live app lets you step through LLL, inspect Gram-Schmidt/Lovasz behavior, generate toy LWE instances, and run LLL/BKZ-based recovery attempts. It does not perform encryption/decryption; it demonstrates reduction-and-recovery attack dynamics for educational analysis. Controls include lattice dimension presets, delta, LWE parameters (`n`, `q`, `sigma`), and BKZ block size (`beta`).

## Exhibits and Learning Goals

| Exhibit | Topic | What you should be able to say afterward |
| --- | --- | --- |
| 1 | What Is a Lattice? | Many bases generate the same lattice; the determinant is invariant under unimodular changes. |
| 2 | Gram-Schmidt | Compute `b*` and `mu`, and evaluate the Lovasz condition from them. |
| 3 | LLL Step-by-Step | Predict size-reduce / swap / advance, and see `reducedBasis = U * original` with `det(U) = +-1` proving the lattice is unchanged. |
| 4 | Break a Toy LWE Instance | Derive why `(-s, e, +-1)` is short in the primal embedding, and recover `s` by genuine reduction — never confusing it with the brute-force baseline. |
| 5 | Parameter Explorer | Explain why tiny parameters fall to LLL/BKZ while Kyber-like dimensions push the required `beta` (and cost) out of reach. |

A **Reproducible Labs** panel adds a seed (deterministic runs), shareable lab links, and named canonical labs; a **Challenges** section offers predict-then-reveal exercises.

## Documentation

Classroom-ready materials live in [`docs/`](docs/):

- [Glossary](docs/glossary.md) — every term, tied to where it appears in the app.
- [Model limitations](docs/model-limitations.md) — exactly what is modeled vs. simplified; read before teaching.
- [Worksheet](docs/worksheet.md) and [answer key](docs/answer-key.md).
- Timed labs: [30-min](docs/lab-30-min.md), [60-min](docs/lab-60-min.md), [90-min](docs/lab-90-min.md).

## Accuracy Boundaries

- **Real:** the LLL reduction, the unimodular-transform invariant, the primal/Kannan embedding, exact small-block SVP enumeration, and short-vector secret recovery.
- **Simplified:** toy dimensions, plain floating-point GSO, no enumeration pruning or projected-block handling, and an intuition-only cost model (`cost ~ 2^(0.292*beta)`).
- **Do not infer:** that LLL/BKZ break real Kyber, or that the brute-force "teaching baseline" is a lattice attack. For real analysis use a lattice estimator (Albrecht et al.) and a production reducer (fplll, G6K). See [docs/model-limitations.md](docs/model-limitations.md).

## References and Further Reading

A short, tiered path (each line says why it is worth reading):

1. **LLL intuition** — any "LLL for beginners" lattice-reduction notes: build a mental model of size reduction and swaps before the formal proofs.
2. **LWE background** — Regev, *On Lattices, Learning with Errors...*: the foundational hardness assumption these schemes rest on.
3. **Primal/Kannan embedding** — surveys of primal lattice attacks on LWE: explains why `(-s, e, 1)` is short and how the embedding is built.
4. **BKZ and cost** — Chen-Nguyen *BKZ 2.0* and the root-Hermite-factor literature: how block size sets reduction strength and attack cost.
5. **Estimators and libraries** — the Albrecht et al. lattice estimator, and `fplll` / `G6K`: the tools used for actual parameter security, not a browser toy.

## How to Run Locally

```bash
git clone https://github.com/systemslibrarian/crypto-lab-lll-break
cd crypto-lab-lll-break
npm install
npm run dev
```

No environment variables are required.

## Part of the Crypto-Lab Suite

One of 60+ live browser demos at [systemslibrarian.github.io/crypto-lab](https://systemslibrarian.github.io/crypto-lab/) — spanning Atbash (600 BCE) through NIST FIPS 203/204/205 (2024).

---

*"Whether you eat or drink, or whatever you do, do all to the glory of God." — 1 Corinthians 10:31*
