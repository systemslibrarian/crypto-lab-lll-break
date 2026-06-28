# crypto-lab-lll-break

## What It Is

This demo shows LLL lattice reduction and toy BKZ reduction on Learning With Errors (LWE) embeddings in the browser. It focuses on the attacker workflow: reduce a basis, search for short vectors, and test whether secret recovery is possible at a chosen parameter set. The core problem illustrated is how SVP-approximation quality affects practical LWE attacks. This is a post-quantum cryptography educational model, not a production cryptanalytic tool and not evidence that LLL breaks real Kyber parameters.

## When to Use It

- Teaching why post-quantum parameters matter in lattice cryptography. This demo makes the relationship between dimension, modulus, noise, and attack feasibility visible.
- Explaining LLL and BKZ mechanics to students or engineers. Step traces and block-improvement logs show what each reduction phase actually changes.
- Comparing toy insecure settings against Kyber-like settings. The same pipeline can be run at small and large parameters to show where attacks stop being effective.
- Demonstrating how primal LWE embeddings are constructed. The matrix and reduction output are shown directly so learners can inspect the attack surface.
- Do NOT use it for real-world security assessment. This browser implementation is intentionally simplified and does not replace specialized lattice estimators or high-performance reduction libraries.

## Live Demo

**[systemslibrarian.github.io/crypto-lab-lll-break](https://systemslibrarian.github.io/crypto-lab-lll-break/)**

The live app lets you step through LLL, inspect Gram-Schmidt/Lovász behavior, generate toy LWE instances, and run LLL/BKZ-based recovery attempts. It does not perform encryption/decryption; it demonstrates reduction-and-recovery attack dynamics for educational analysis. Controls include lattice dimension presets, delta, LWE parameters (`n`, `q`, `sigma`), and BKZ block size (`beta`).

## What Can Go Wrong

- **Under-parameterized LWE** — a small dimension `n`, small modulus `q`, or low noise `sigma` lets LLL/BKZ recover the secret quickly. Security comes from the combination of parameters, not any single knob.
- **Mistaking LLL for the real attack** — LLL only guarantees an exponential approximation factor; real cryptanalysis uses BKZ with a large block size `beta`, and cost is measured with a lattice estimator, not one LLL run.
- **Too little noise** — if the error is tiny relative to `q`, the embedded unique-shortest-vector gap is large and even weak reduction finds the secret.
- **Treating a browser toy as an estimator** — production parameter choices need vetted tooling (the lattice estimator, sieving libraries), not a simplified in-browser reducer.
- **Assuming "post-quantum" means "unbreakable"** — lattice schemes are only as strong as their parameters; the lesson here is that bad parameters fall to ordinary classical reduction.

## Real-World Usage

- **PQC parameter selection** — Kyber / ML-KEM, Dilithium / ML-DSA, and Falcon parameters are chosen so the best known BKZ attack costs at least the target security level.
- **The lattice estimator** — the standard tool (Albrecht et al.) models exactly this reduce-and-recover attack to price LWE and NTRU instances during standardization.
- **BKZ and SVP records** — progressive BKZ, BKZ 2.0, and lattice sieving (G6K) drive the Darmstadt SVP/LWE challenges that calibrate real-world attack cost.
- **LLL beyond LWE** — Coppersmith's method (RSA small-root and Håstad attacks), knapsack and early lattice-scheme breaks, and integer-relation detection all rely on LLL reduction.

## How to Run Locally

```bash
git clone https://github.com/systemslibrarian/crypto-lab-lll-break
cd crypto-lab-lll-break
npm install
npm run dev
```

## Related Demos

- [crypto-lab-nonce-lattice](https://systemslibrarian.github.io/crypto-lab-nonce-lattice/) — lattice attack recovering ECDSA keys from biased nonces (the Hidden Number Problem).
- [crypto-lab-lwe-hints](https://systemslibrarian.github.io/crypto-lab-lwe-hints/) — how side-channel "hints" weaken an LWE instance.
- [crypto-lab-kyber-vault](https://systemslibrarian.github.io/crypto-lab-kyber-vault/) — ML-KEM (Kyber), the lattice KEM whose parameters this attack pressures.
- [crypto-lab-frodo-vault](https://systemslibrarian.github.io/crypto-lab-frodo-vault/) — FrodoKEM, plain (unstructured) LWE for comparison.

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

---

*One of 120+ browser demos in the [Crypto Lab](https://crypto-lab.systemslibrarian.dev/) suite.*

*"So whether you eat or drink or whatever you do, do it all for the glory of God." — 1 Corinthians 10:31*
