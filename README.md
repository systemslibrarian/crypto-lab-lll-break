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
