# crypto-lab-lll-break

## What It Is

This demo shows LLL lattice reduction and toy BKZ reduction on Learning With Errors (LWE) embeddings in the browser. It focuses on the attacker workflow: reduce a basis, search for short vectors, and test whether secret recovery is possible at a chosen parameter set. The core problem illustrated is how SVP-approximation quality affects practical LWE attacks. This is a post-quantum cryptography educational model, not a production cryptanalytic tool and not evidence that LLL breaks real Kyber parameters. A sixth exhibit turns the same reducer on a historical target: a toy Merkle-Hellman trapdoor knapsack, where LLL recovers message bits from public data alone, and where the lesson is why NP-completeness never secured the scheme.

## When to Use It

- Teaching why post-quantum parameters matter in lattice cryptography. This demo makes the relationship between dimension, modulus, noise, and attack feasibility visible.
- Explaining LLL and BKZ mechanics to students or engineers. Step traces and block-improvement logs show what each reduction phase actually changes.
- Comparing toy insecure settings against Kyber-like settings. The same pipeline can be run at small and large parameters to show where attacks stop being effective.
- Demonstrating how primal LWE embeddings are constructed. The matrix and reduction output are shown directly so learners can inspect the attack surface.
- Do NOT use it for real-world security assessment. This browser implementation is intentionally simplified and does not replace specialized lattice estimators or high-performance reduction libraries.

## Live Demo

**[systemslibrarian.github.io/crypto-lab-lll-break](https://systemslibrarian.github.io/crypto-lab-lll-break/)**

The live app lets you step through LLL, inspect Gram-Schmidt/Lovász behavior, generate toy LWE instances, and run LLL/BKZ-based recovery attempts. The LWE exhibits do not perform encryption/decryption; they demonstrate reduction-and-recovery attack dynamics for educational analysis. Exhibit 6 is the deliberate exception, because contrasting a legitimate Merkle-Hellman trapdoor decryption with a public-only lattice attack on the same ciphertext is the whole point of it. Controls include lattice dimension presets, delta, LWE parameters (`n`, `q`, `sigma`), and BKZ block size (`beta`).

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
| 1 | What Is a Lattice? | Many bases generate the same lattice; the determinant is invariant under unimodular changes. Opens with a plain-language primer (lattice / basis / "many bases, one lattice") before any jargon, and the determinant label carries a hover/click tooltip: "area of one grid cell — it never changes no matter which basis you pick." |
| 2 | Gram-Schmidt | Compute `b*` and `mu`, and evaluate the Lovasz condition from them. |
| 3 | LLL Step-by-Step | Predict size-reduce / swap / advance, and see `reducedBasis = U * original` with `det(U) = +-1` proving the lattice is unchanged. The green (current) basis is animated so you watch the vectors shrink and re-order per step. A `delta` vs `beta` callout links this exhibit's LLL knob to Exhibit 4's BKZ knob. |
| 4 | Break a Toy LWE Instance | Derive why `(-s, e, +-1)` is short in the primal embedding, and recover `s` by genuine reduction — never confusing it with the brute-force baseline. A color-coded block view of the embedding matrix (secret / error-q-ary / embed columns) animates the bottom-row-minus-`sum(s_i·row_i)` cancellation so the short vector visibly emerges; the norm-gap confidence meter carries a "what am I looking at" caption. |
| 5 | Parameter Explorer | Explain why tiny parameters fall to LLL/BKZ while Kyber-like dimensions push the required `beta` (and cost) out of reach. The `n≈50` line and `beta ≈ n/(2·log2(q/sigma))` formula are labelled on-screen as a rough teaching model (not a security estimate — real numbers come from the lattice estimator), with a one-line derivation tooltip for `beta`. |
| 6 | Knapsack Trapdoors and the Limits of Worst-Case Hardness | Build a toy Merkle-Hellman trapdoor knapsack, encrypt with it, decrypt with the trapdoor, then recover the message bits from **public data alone** by handing a Lagarias-Odlyzko lattice to the same LLL engine. Say why "subset sum is NP-complete" never protected it: NP-completeness is about the *hardest* instances, and a cryptosystem is a *distribution*. |

Every canvas carries a persistent legend (red = your input basis, green = current/reduced basis, dashed cyan = Gram-Schmidt directions, dots = lattice points) so the visuals read without reverse-engineering the code.

A **Reproducible Labs** panel adds a seed (deterministic runs), shareable lab links, and named canonical labs; a **Challenges** section offers predict-then-reveal exercises.

## Exhibit 6 in detail: knapsack trapdoors, and what NP-completeness never bought

**Merkle-Hellman is historical and broken. Do not use any of this for anything.**

The exhibit builds a real toy Merkle-Hellman key, encrypts a bit vector under it,
decrypts with the trapdoor, and then attacks the ciphertext with nothing but the
public weights.

### Exactly what the attack does and does not do

- It recovers **message bits** from one ciphertext, in selected toy low-density
  instances.
- It does **not** reconstruct the original private trapdoor or find an equivalent
  decryption pair. The modulus `m`, multiplier `w` and superincreasing sequence
  `a'` remain unknown to this attack.
- It is **not** Shamir's attack. Shamir's attack seeks a usable decryption pair
  `(W*, M*)`, which need not equal the pair used by the encrypter; that is what
  historically broke the basic scheme, and it is not implemented here.
  Lagarias and Odlyzko draw that same line themselves, describing their method
  as complementary to attacks "based on recovering trapdoor information".
- **Failure to recover is not a security result.** It is one browser-scale LLL
  run on one toy instance. A failure is reachable and reproducible in the app:
  seed `30`, `n = 16`, `k = 1`, message `0110100101100101`. The trapdoor still
  decrypts that instance perfectly, so what failed is the attack and not the
  scheme. In the fixed seeds 1-80 sweep, 7 fail at `k = 1`, 2 fail at `k = 3`,
  and none fail at `k = 5`. Those are measured counts at one dimension, not a
  threshold or a claim of monotonic behavior for each seed.
- A small experiment **does not establish a density threshold**, and this lab
  does not display one. Even an in-budget `n = 8` instance can fail here: seed
  `10`, `k = 1`, message `00001100` makes the trapdoor decrypt exactly while
  the two-pass public attack finds no candidate. That bounded counterexample
  means this experiment does not isolate either dimension or density as the
  binding cause of success or failure.
- NP-completeness of generic subset sum says nothing about this structured
  distribution. Nothing here establishes or refutes worst-case lattice hardness,
  and nothing here bears on ML-KEM/Kyber or any production parameter set.

### The density convention, stated exactly

```
d = n / log2( max_i a_i )        over the PUBLIC weights
```

That is the convention of Lagarias-Odlyzko (1985, Section 1) and, identically,
Coster et al. (1992, equation 1.2). Lagarias-Odlyzko state `0.645` for their
shortest-vector result; Coster et al. later give `0.9408...` for their modified
embedding. Turning either shortest-vector statement into an attack bound is
**conditional on a shortest-vector oracle**, which LLL is not; the rigorous
Lagarias-Odlyzko result for their actual LLL route is of order `1/n`. Those
theorems sample *independently uniform* weights. Merkle-Hellman's public weights
are instead generated jointly as `w * a' mod m` from one superincreasing `a'`
and shared `w,m`. The figures are quoted as history, not borrowed as authority.

### The embedding

The lattice is the **original Lagarias-Odlyzko Algorithm SV basis (2.4), exactly
as printed** - no scaling constant, no centring:

```
b_1     = (1, 0, ..., 0, -a_1)
...
b_n     = (0, 0, ..., 1, -a_n)
b_{n+1} = (0, 0, ..., 0,  M)
```

This is **not** the Coster et al. embedding and is not labelled as one. Coster
et al. restate Lagarias-Odlyzko with a scaling constant `N` in their Section 2,
and their Section 3 modification - the one carrying the `0.9408...` bound -
centres the coordinates at one half. Neither is implemented here.

Candidate extraction is the paper's own step (3): a reduced row is a candidate
only if each of its first `n` coordinates is either `0` or one fixed non-zero `lambda`,
and a candidate is accepted only when `sum(bits_i * a_i) == M` holds exactly in
integer arithmetic. The search is bounded and disclosed - `n+1` rows per pass,
two passes (step (4) repeats on `sum(a_i) - M`), so at most `2(n+1)` candidates
are ever tested. **There is no brute-force fallback.** The exhaustive enumerator
in the module exists only to answer whether an instance is *ambiguous*, is
labelled as verifier-side ground truth, and is never reported as a recovery.

### What the attacker is given

`attackPublicOnly()` receives exactly `{ publicWeights, target, delta }`. It
cannot reach the private sequence, the modulus, the multiplier, its inverse or
the plaintext, and ground-truth comparison happens in the caller after it
returns. Two checks enforce that rather than asserting it: the attack is run
against an instance built from pinned public literals with no private key in
existence, and it is run through a `Proxy` that throws on any property read
outside those three.

### Arithmetic

Key generation, encryption, decryption and every candidate verification use
`BigInt`. Algorithm SV is run only for positive weights and its published
domain `1 <= M < sum(a_i)`. The exhibit rejects an all-zero message before key
generation. It can still encrypt and trapdoor-decrypt the all-one endpoint, but
the public lattice route refuses that endpoint by name rather than presenting it
as an Algorithm SV run.

The lattice goes to this lab's existing `number[][]` reducer. Before either
pass, the input screen requires `(n+1) * E^2 <= Number.MAX_SAFE_INTEGER`, where
`E = max(max a_i, M, sum(a_i)-M)`. Including the complement target is
load-bearing because Algorithm SV step (4) reduces a second basis. This screen
keeps both input bases exact, but it does **not** pretend every intermediate
Gram-Schmidt value is exact: `mu` uses division, and size reduction can grow
working-basis entries.

After each pass, a separate BigInt postcheck reconstructs the returned basis as
`U * inputBasis` and requires `det(U) = +/-1`. A failed postcheck is a refusal,
not a result. Thus floating-point decisions can change which exact basis is
returned or cause no candidate to be found, while no row is accepted unless it
is proved to belong to the original lattice and its decoded bits independently
satisfy the public subset sum in `BigInt`.

One consequence is a limit of the model rather than a fact about knapsacks: low
density needs large weights, and large weights spend that budget quadratically,
so low density and large `n` cannot be had together here.

### Primary sources

- R. C. Merkle and M. E. Hellman, "Hiding Information and Signatures in Trapdoor
  Knapsacks", *IEEE Transactions on Information Theory* **IT-24**(5), September
  1978, pp. 525-530. doi:[10.1109/TIT.1978.1055927](https://doi.org/10.1109/TIT.1978.1055927).
  Author copy: [www-ee.stanford.edu/~hellman/publications/30.pdf](https://www-ee.stanford.edu/~hellman/publications/30.pdf).
  Sections used: I p. 525, equations (1) and (2); II p. 526, equations (3)-(8)
  and the n=5 worked example, which this repository reproduces as a pinned
  known-answer check.
- J. C. Lagarias and A. M. Odlyzko, "Solving Low-Density Subset Sum Problems",
  *Journal of the ACM* **32**(1), January 1985, pp. 229-246.
  doi:[10.1145/2455.2461](https://doi.org/10.1145/2455.2461). Sections used: 1
  p. 230 (density); 2 pp. 232-233 (Algorithm SV, basis (2.4), Propositions
  2.1-2.2); 3 pp. 233-234 (the sample space and Theorems 3.3, 3.5).
- J. C. Lagarias, "Performance Analysis of Shamir's Attack on the Basic
  Merkle-Hellman Knapsack Cryptosystem", ICALP 1984, pp. 312-323.
  [Author copy](https://dept.math.lsa.umich.edu/~lagarias/doc/1218attack.pdf),
  Section 2: the goal is a decryption pair `(W*, M*)` that need not be the
  encrypter's original `(W, M)`. Cited only to state the boundary between
  Shamir's attack and the Algorithm SV lesson accurately.
- A. K. Lenstra, H. W. Lenstra Jr. and L. Lovasz, "Factoring Polynomials with
  Rational Coefficients", *Mathematische Annalen* **261**, 1982, pp. 515-534.
  doi:[10.1007/BF01457454](https://doi.org/10.1007/BF01457454). Sections used:
  1 pp. 516-518, the reduced-basis definition (1.4)-(1.5) and Propositions
  (1.11) and (1.12). (1.12) is the licence for inspecting every reduced row
  rather than only the first.
- M. J. Coster, A. Joux, B. A. LaMacchia, A. M. Odlyzko, C.-P. Schnorr and
  J. Stern, "Improved Low-Density Subset Sum Algorithms", *Computational
  Complexity* **2**, 1992, pp. 111-128.
  doi:[10.1007/BF01201999](https://doi.org/10.1007/BF01201999). Cited for the
  density definition, the oracle-conditional bounds, and the fact that its
  embeddings differ from the one implemented here.

The two load-bearing inequalities were read off **rendered images of the printed
pages**, not off extracted text: the OCR of both papers renders the `>=` in
Merkle-Hellman (2) as `>` and drops the prime on `a'` in (8). Both corruptions
change the mathematics. The paper's own n=5 example independently settles the
first - its step-2 residue is exactly `a'_2 = 196` and its printed answer is
`x_2 = 1`, which a strict `>` cannot produce.

## Documentation

Classroom-ready materials live in [`docs/`](docs/):

- [Glossary](docs/glossary.md) — every term, tied to where it appears in the app.
- [Model limitations](docs/model-limitations.md) — exactly what is modeled vs. simplified; read before teaching.
- [Worksheet](docs/worksheet.md) and [answer key](docs/answer-key.md).
- Timed labs: [30-min](docs/lab-30-min.md), [60-min](docs/lab-60-min.md), [90-min](docs/lab-90-min.md).

---

*Part of the [Crypto Lab](https://crypto-lab.systemslibrarian.dev/) suite.*

*"So whether you eat or drink or whatever you do, do it all for the glory of God." — 1 Corinthians 10:31*
