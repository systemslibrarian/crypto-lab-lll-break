## What It Is

Browser-based demonstration of the LLL lattice reduction algorithm
(Lenstra, Lenstra, Lovasz, 1982) and its successor BKZ - the tools
attackers use against lattice-based cryptography.

Every post-quantum system in the crypto-lab suite (ML-KEM, FrodoKEM,
ML-DSA, FALCON) claims security based on the hardness of lattice problems.
This demo makes that claim tangible: watch LLL reduce a bad basis to a
good one step-by-step, embed a toy LWE instance into a lattice and
recover the secret in milliseconds, then watch the exact same attack fail
catastrophically at Kyber's actual parameters.

Implements the full LLL algorithm with exact Gram-Schmidt orthogonalization,
Lovasz condition checking, and step-by-step animation. Includes a toy LWE
attacker using the primal embedding and BKZ reduction, and a parameter
explorer showing exactly where and why security emerges.

## When to Use It

- Understanding why post-quantum parameter choices (n, q, sigma) are not arbitrary
- Seeing what "lattice problems are hard" actually means operationally
- Teaching the connection between SVP hardness and LWE security
- Exploring the gap between toy insecure parameters and real deployments

## Live Demo

https://systemslibrarian.github.io/crypto-lab-lll-break/

## GitHub Pages Deployment

1. Push this repository to the `main` branch on GitHub.
2. In repository settings, open Pages and set Source to GitHub Actions.
3. Ensure Actions are enabled for the repository.
4. The workflow at `.github/workflows/deploy.yml` will build with Vite and publish `dist`.
5. Site URL will be: `https://systemslibrarian.github.io/crypto-lab-lll-break/`.

This project already sets `base: '/crypto-lab-lll-break/'` in `vite.config.ts`,
which is required for correct asset paths on GitHub Pages project sites.

## What Can Go Wrong

- LLL uses floating-point Gram-Schmidt coefficients. For high-dimensional
  or poorly conditioned lattices, numerical precision can affect results.
  This is a known limitation of floating-point LLL - production implementations
  use multiple-precision arithmetic (MPFR). The demo documents this.
- BKZ with large block sizes (beta > 8) is not browser-feasible. The demo
  simulates the cost rather than running it. Real BKZ-400 would take years.
- The simplified security estimator (beta ~= n / (2*log2(q/sigma))) is an
  approximation. The exact lattice estimator (lattice-estimator.readthedocs.io)
  is far more precise but requires Python.

## Real-World Usage

LLL has broken: Merkle-Hellman knapsack cryptosystem (1983), RSA with
small private exponent (Boneh-Durfee), truncated linear congruential generators,
and DSA with biased nonces. BKZ is the current state-of-the-art for lattice
attacks and is used to set security parameters for all NIST PQC standards.

The required BKZ block size beta for breaking Kyber-512 is approximately 400,
with attack cost ~2^117. NIST's security analysis confirms this. The demo
shows this failure explicitly - not just claims it.
