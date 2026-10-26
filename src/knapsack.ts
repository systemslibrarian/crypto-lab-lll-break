/**
 * Exhibit 6 — Merkle-Hellman trapdoor knapsacks, and the Lagarias-Odlyzko
 * low-density lattice attack that reads the message bits back out of PUBLIC
 * data alone.
 *
 * This file implements two things from their primary sources and nothing else:
 *
 *   The cryptosystem — R. C. Merkle and M. E. Hellman, "Hiding Information and
 *   Signatures in Trapdoor Knapsacks", IEEE Transactions on Information Theory
 *   IT-24(5), September 1978, pp. 525-530. Section II, p. 526, equations
 *   (1)-(8), and the n=5 worked example on the same page.
 *
 *   The attack — J. C. Lagarias and A. M. Odlyzko, "Solving Low-Density Subset
 *   Sum Problems", Journal of the ACM 32(1), January 1985, pp. 229-246.
 *   Algorithm SV, Section 2, pp. 232-233, including the basis printed at (2.4)
 *   and the candidate test printed at step (3).
 *
 * The reducer is this lab's own `lllReduce` — the same engine Exhibits 3 and 4
 * drive. Nothing here reimplements lattice reduction.
 *
 * WHAT THIS IS NOT:
 *   - NOT Shamir's attack. Shamir recovers the TRAPDOOR (w, m, a'). This
 *     recovers MESSAGE BITS from one ciphertext and learns nothing about the
 *     trapdoor. Lagarias-Odlyzko draw that line themselves (Section 1, p. 230),
 *     calling their method complementary to trapdoor-recovery attacks.
 *   - NOT the Coster et al. embedding. The basis below is the ORIGINAL
 *     Lagarias-Odlyzko (2.4): no scaling constant N, no centring at 1/2.
 *     Coster et al. (Computational Complexity 2, 1992, pp. 111-128) restate
 *     Lagarias-Odlyzko WITH a constant N in their Section 2, and their Section 3
 *     modification — the one carrying the 0.9408... bound — centres the
 *     coordinates. Neither is implemented here.
 *   - NOT a security result in either direction. Merkle-Hellman is historical
 *     and broken. A failure of this browser-scale attack says nothing about any
 *     scheme's security, and nothing here bears on ML-KEM/Kyber or on
 *     worst-case lattice hardness.
 */
import { lllReduce } from './lll';
import { nextUint32 } from './rng';

export interface MerkleHellmanPrivateKey {
  /** The superincreasing sequence a'_1 .. a'_n of Merkle-Hellman (1). */
  readonly privateWeights: readonly bigint[];
  /** The modulus m, which (8) requires to exceed the sum of a'. */
  readonly modulus: bigint;
  /** The multiplier w of (3), invertible modulo m. */
  readonly multiplier: bigint;
  /** w^-1 mod m, used by (4). */
  readonly multiplierInverse: bigint;
}

export interface MerkleHellmanKeyPair {
  readonly privateKey: MerkleHellmanPrivateKey;
  /** a_i = w * a'_i mod m, Merkle-Hellman (3). This is the public key. */
  readonly publicWeights: readonly bigint[];
}

/**
 * Everything the lattice attacker is given. There is deliberately nothing else
 * in this object: no modulus, no multiplier, no private sequence, no plaintext.
 */
export interface PublicKnapsackInstance {
  readonly publicWeights: readonly bigint[];
  /** The ciphertext S = sum(x_i * a_i). Lagarias-Odlyzko call it M. */
  readonly target: bigint;
  /** LLL's delta. Lenstra-Lenstra-Lovasz (1.5) allows any fixed 1/4 < d < 1. */
  readonly delta: number;
}

export type CandidateRejectionReason =
  | 'coordinates-outside-bit-set'
  | 'public-sum-mismatch';

export interface CandidateRejection {
  readonly pass: 1 | 2;
  readonly row: number;
  readonly reason: CandidateRejectionReason;
  readonly detail: string;
}

export interface KnapsackAttackResult {
  readonly outcome: 'recovered' | 'no-candidate' | 'refused';
  /** The recovered bits, or null. Never a guess: see `verifyCandidate`. */
  readonly bits: readonly number[] | null;
  /** The fixed non-zero lambda of Algorithm SV step (3). */
  readonly lambda: number | null;
  /** 1 = the target as given; 2 = Algorithm SV step (4), the complement. */
  readonly pass: 1 | 2 | null;
  readonly row: number | null;
  readonly winningRow: readonly number[] | null;
  /** Total reduced-basis rows examined across both passes. Bounded by 2(n+1). */
  readonly rowsInspected: number;
  readonly rejections: readonly CandidateRejection[];
  /**
   * The exact basis handed to LLL on PASS 1 only. Pass 2 (Algorithm SV step
   * (4)) builds a second, different basis on sum(a_i) - M; it is deliberately
   * not carried here, so `winningRow` may come from a basis this field does not
   * describe. `pass` says which, and the UI says so too rather than letting a
   * reader assume the row is a row of the basis on screen.
   */
  readonly basis: readonly (readonly number[])[];
  /** What LLL returned on pass 1, with the same caveat as `basis`. */
  readonly reducedBasis: readonly (readonly number[])[];
  /** Pass 1's swap count. */
  readonly swapCount: number;
  /**
   * True only after independently checking, in BigInt, that every returned
   * reduced row equals U times the input basis and det(U) = +/-1. Floating
   * Gram-Schmidt chooses the row operations; this postcondition proves those
   * operations still returned a basis of the same integer lattice.
   */
  readonly exactLatticeTransformVerified: boolean;
  /** Set only when outcome === 'refused'. */
  readonly refusal: string | null;
}

export interface ArithmeticBudget {
  /** E = max(max public weight, M, sum(a_i) - M): both Algorithm SV passes. */
  readonly maxEntry: number;
  /** (n+1) * E^2 — a conservative bound for either input basis's row norms. */
  readonly required: number;
  readonly limit: number;
  readonly withinBudget: boolean;
  /** limit / required, as a plain ratio. Below 1 the instance is refused. */
  readonly headroom: number;
}

/** Exhaustive enumeration is only honest at toy sizes; past this it refuses. */
export const MAX_ENUMERATION_BITS = 16;

/* ------------------------------------------------------------------------ */
/* Exact integer helpers. Everything a verdict depends on runs through these. */
/* ------------------------------------------------------------------------ */

function extendedGcd(a: bigint, b: bigint): { g: bigint; x: bigint; y: bigint } {
  let [oldR, r] = [a, b];
  let [oldS, s] = [1n, 0n];
  let [oldT, t] = [0n, 1n];
  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
    [oldT, t] = [t, oldT - q * t];
  }
  return { g: oldR, x: oldS, y: oldT };
}

/** w^-1 mod m, or null when w is not invertible (gcd(w, m) != 1). */
export function modInverse(w: bigint, m: bigint): bigint | null {
  if (m <= 1n) return null;
  const a = ((w % m) + m) % m;
  if (a === 0n) return null;
  const { g, x } = extendedGcd(a, m);
  if (g !== 1n) return null;
  return ((x % m) + m) % m;
}

export function sumOf(weights: readonly bigint[]): bigint {
  let total = 0n;
  for (const w of weights) total += w;
  return total;
}

/**
 * Merkle-Hellman (1): a'_i > sum of all earlier a'_j, STRICTLY, for every i.
 * Returns the first index that violates it, or -1.
 */
export function firstNonSuperincreasingIndex(weights: readonly bigint[]): number {
  let running = 0n;
  for (let i = 0; i < weights.length; i += 1) {
    if (weights[i] <= running) return i;
    running += weights[i];
  }
  return -1;
}

export function isSuperincreasing(weights: readonly bigint[]): boolean {
  return weights.length > 0 && firstNonSuperincreasingIndex(weights) === -1;
}

/* ------------------------------------------------------------------------ */
/* Key validation. Each of the three conditions is checked separately so a     */
/* failure names which one broke, and each has its own owning test.            */
/* ------------------------------------------------------------------------ */

export function validatePrivateKey(key: MerkleHellmanPrivateKey): void {
  const { privateWeights, modulus, multiplier, multiplierInverse } = key;

  if (privateWeights.length === 0) {
    throw new Error('A private knapsack needs at least one weight');
  }
  for (let i = 0; i < privateWeights.length; i += 1) {
    if (privateWeights[i] <= 0n) {
      throw new Error(`Private weight a'${i + 1} must be positive`);
    }
  }

  // Merkle-Hellman (1).
  const bad = firstNonSuperincreasingIndex(privateWeights);
  if (bad !== -1) {
    throw new Error(
      `Private sequence is not superincreasing: a'${bad + 1} = ${privateWeights[bad]} ` +
        `does not exceed the sum ${sumOf(privateWeights.slice(0, bad))} of the weights before it ` +
        '(Merkle-Hellman 1978, equation (1))',
    );
  }

  // Merkle-Hellman (8). Without it the modular sum in (7) is not the integer
  // sum, and the greedy pass in (2) decodes a different message.
  const total = sumOf(privateWeights);
  if (modulus <= total) {
    throw new Error(
      `Modulus m = ${modulus} must exceed the private-weight sum ${total} ` +
        '(Merkle-Hellman 1978, equation (8))',
    );
  }

  // Merkle-Hellman (3) needs w invertible mod m.
  const inverse = modInverse(multiplier, modulus);
  if (inverse === null) {
    throw new Error(
      `Multiplier w = ${multiplier} is not invertible modulo m = ${modulus}; ` +
        'Merkle-Hellman 1978 requires gcd(w, m) = 1',
    );
  }
  if (inverse !== multiplierInverse) {
    throw new Error(
      `Stored multiplier inverse ${multiplierInverse} is not w^-1 mod m (expected ${inverse})`,
    );
  }
}

/** Merkle-Hellman (3): a_i = w * a'_i mod m, for every i. */
export function derivePublicWeights(key: MerkleHellmanPrivateKey): bigint[] {
  validatePrivateKey(key);
  return key.privateWeights.map((a) => (key.multiplier * a) % key.modulus);
}

/* ------------------------------------------------------------------------ */
/* Key generation                                                             */
/* ------------------------------------------------------------------------ */

/** Uniform in [lo, hi] by rejection sampling, so no modulo bias. */
function randomBigInt(lo: bigint, hi: bigint, next: () => number): bigint {
  if (hi < lo) throw new Error('Empty range');
  const span = hi - lo + 1n;
  const bits = span.toString(2).length;
  const words = Math.ceil(bits / 32);
  const mask = (1n << BigInt(bits)) - 1n;
  for (let guard = 0; guard < 10_000; guard += 1) {
    let out = 0n;
    for (let i = 0; i < words; i += 1) out = (out << 32n) | BigInt(next() >>> 0);
    out &= mask;
    if (out < span) return lo + out;
  }
  throw new Error('Rejection sampling failed to terminate');
}

/**
 * Generate a toy Merkle-Hellman key pair.
 *
 * The private sequence follows the paper's own prescription (p. 526): a'_i is
 * drawn uniformly from [(2^(i-1) - 1) * 2^k + 1, 2^(i-1) * 2^k]. That range
 * GUARANTEES (1), since the earlier weights sum to at most 2^k(2^(i-1) - 1).
 * The paper uses k = 100; this lab uses a small k so the lattice stays inside
 * the exact-arithmetic budget (see `arithmeticBudget`). `validatePrivateKey`
 * still checks (1) independently — a generator that guarantees an invariant and
 * a checker that confirms it fail in different ways.
 */
export function generateKeyPair(
  n: number,
  scaleBits: number,
  next: () => number = nextUint32,
): MerkleHellmanKeyPair {
  if (!Number.isInteger(n) || n < 2 || n > 24) {
    throw new Error('Knapsack length n must be an integer in [2, 24]');
  }
  if (!Number.isInteger(scaleBits) || scaleBits < 1 || scaleBits > 40) {
    throw new Error('Private-weight scale must be an integer number of bits in [1, 40]');
  }

  const scale = 1n << BigInt(scaleBits);
  const privateWeights: bigint[] = [];
  for (let i = 1; i <= n; i += 1) {
    const step = 1n << BigInt(i - 1);
    const lo = (step - 1n) * scale + 1n;
    const hi = step * scale;
    privateWeights.push(randomBigInt(lo, hi, next));
  }

  const total = sumOf(privateWeights);
  // (8), with room above it so m is not pinned to the sum.
  const modulus = randomBigInt(total + 1n, 2n * total, next);

  let multiplier = 0n;
  let multiplierInverse: bigint | null = null;
  for (let guard = 0; guard < 10_000 && multiplierInverse === null; guard += 1) {
    multiplier = randomBigInt(2n, modulus - 2n, next);
    multiplierInverse = modInverse(multiplier, modulus);
  }
  if (multiplierInverse === null) {
    throw new Error('Could not find a multiplier invertible modulo m');
  }

  const privateKey: MerkleHellmanPrivateKey = {
    privateWeights,
    modulus,
    multiplier,
    multiplierInverse,
  };
  return { privateKey, publicWeights: derivePublicWeights(privateKey) };
}

/* ------------------------------------------------------------------------ */
/* Encrypt / decrypt                                                          */
/* ------------------------------------------------------------------------ */

function assertBits(bits: readonly number[], n: number): void {
  if (bits.length !== n) {
    throw new Error(`Expected ${n} message bits, got ${bits.length}`);
  }
  for (let i = 0; i < bits.length; i += 1) {
    if (bits[i] !== 0 && bits[i] !== 1) {
      throw new Error(`Message bit ${i + 1} must be 0 or 1, got ${bits[i]}`);
    }
  }
}

/**
 * S = sum(x_i * a_i), the paper's S = a*x.
 *
 * BIT ORDER: bits[i] multiplies publicWeights[i], i = 0..n-1, which is index
 * i+1 in both papers. Bits read left to right in ASCENDING weight index. That
 * is not most-significant-first, and it is not the order the greedy loop in
 * `decryptWithTrapdoor` runs (that runs i = n down to 1).
 */
export function encryptBits(publicWeights: readonly bigint[], bits: readonly number[]): bigint {
  assertBits(bits, publicWeights.length);
  let target = 0n;
  for (let i = 0; i < publicWeights.length; i += 1) {
    if (bits[i] === 1) target += publicWeights[i];
  }
  return target;
}

export interface TrapdoorDecryption {
  readonly bits: readonly number[];
  /** S' = w^-1 * S mod m, Merkle-Hellman (4). */
  readonly transformedSum: bigint;
  /** What the greedy pass could not account for. Zero on a well-formed decrypt. */
  readonly remainder: bigint;
  readonly consistent: boolean;
}

/**
 * Merkle-Hellman (4) then (2).
 *
 * The inequality in (2) is `>=`, not `>`. The paper's own n=5 example settles
 * it: with a' = (171, 196, 457, 1191, 2410) and S' = 3797 the residue at step 2
 * is exactly 196 = a'_2, and the printed answer x_2 = 1 is reachable only with
 * `>=`. The rendered page carries the same glyph.
 */
export function decryptWithTrapdoor(
  key: MerkleHellmanPrivateKey,
  target: bigint,
): TrapdoorDecryption {
  validatePrivateKey(key);
  const n = key.privateWeights.length;
  const transformedSum = (key.multiplierInverse * (((target % key.modulus) + key.modulus) % key.modulus)) % key.modulus;

  let residue = transformedSum;
  const bits = new Array<number>(n).fill(0);
  for (let i = n - 1; i >= 0; i -= 1) {
    if (residue >= key.privateWeights[i]) {
      bits[i] = 1;
      residue -= key.privateWeights[i];
    }
  }
  return {
    bits,
    transformedSum,
    remainder: residue,
    consistent: residue === 0n,
  };
}

/* ------------------------------------------------------------------------ */
/* Density, and the exact-arithmetic budget                                   */
/* ------------------------------------------------------------------------ */

/**
 * d = n / log2(max_i a_i), over the PUBLIC weights.
 *
 * Lagarias-Odlyzko Section 1, p. 230, and identically Coster et al. (1.2).
 * The exhibit prints this formula beside the number it produces.
 */
export function density(publicWeights: readonly bigint[]): number {
  if (publicWeights.length === 0) return 0;
  let max = publicWeights[0];
  for (const a of publicWeights) if (a > max) max = a;
  if (max <= 1n) return Number.POSITIVE_INFINITY;
  return publicWeights.length / Math.log2(Number(max));
}

/**
 * The lattice is handed to a `number[][]` reducer with a floating-point
 * Gram-Schmidt, so both Algorithm SV input bases have to fit in exact integers:
 *
 *   1. every input entry is a safe integer;
 *   2. (n+1) * E^2 <= MAX_SAFE_INTEGER, so each input row's squared norm is
 *      exact before the floating-point Gram-Schmidt work begins.
 *
 * E includes M AND sum(a_i)-M because Algorithm SV step (4) reduces a second
 * basis. Measuring only M can approve a second pass that is actually outside
 * the declared boundary.
 *
 * This is an INPUT screen, not a claim that every intermediate dot product is
 * exact. `mu` uses division, and size reduction can grow working-basis entries.
 * After each pass, `attackPublicOnly` therefore verifies in BigInt that the
 * returned basis is exactly U times the input and det(U)=+/-1. A failed
 * postcondition is a refusal. A surviving candidate must separately satisfy the
 * public subset sum in BigInt.
 */
export function arithmeticBudget(
  publicWeights: readonly bigint[],
  target: bigint,
): ArithmeticBudget {
  const n = publicWeights.length;
  const total = sumOf(publicWeights);
  const complementTarget = total - target;
  const absolute = (value: bigint) => (value < 0n ? -value : value);
  let max = absolute(target);
  if (absolute(complementTarget) > max) max = absolute(complementTarget);
  for (const a of publicWeights) if (absolute(a) > max) max = absolute(a);
  const maxEntry = Number(max);
  const required = (n + 1) * maxEntry * maxEntry;
  const limit = Number.MAX_SAFE_INTEGER;
  const entriesExact = Number.isSafeInteger(maxEntry);
  return {
    maxEntry,
    required,
    limit,
    withinBudget: entriesExact && Number.isSafeInteger(required) && required <= limit,
    headroom: required > 0 ? limit / required : Number.POSITIVE_INFINITY,
  };
}

/* ------------------------------------------------------------------------ */
/* The Lagarias-Odlyzko lattice                                               */
/* ------------------------------------------------------------------------ */

/**
 * Algorithm SV basis (2.4), EXACTLY as printed:
 *
 *   b_1     = (1, 0, ..., 0, -a_1)
 *   b_2     = (0, 1, ..., 0, -a_2)
 *   ...
 *   b_n     = (0, 0, ..., 1, -a_n)
 *   b_{n+1} = (0, 0, ..., 0,  M)
 *
 * No scaling constant, no centring. The Coster et al. restatement uses +N*a_i
 * and N*s; that is a different basis and is not what this builds.
 */
export function lagariasOdlyzkoBasis(
  publicWeights: readonly bigint[],
  target: bigint,
): number[][] {
  const n = publicWeights.length;
  const toSafeInteger = (value: bigint, label: string): number => {
    const converted = Number(value);
    if (!Number.isSafeInteger(converted)) {
      throw new Error(`${label} does not fit exactly in a JavaScript safe integer`);
    }
    return converted;
  };
  const basis: number[][] = [];
  for (let i = 0; i < n; i += 1) {
    const row = new Array<number>(n + 1).fill(0);
    row[i] = 1;
    row[n] = -toSafeInteger(publicWeights[i], `public weight a_${i + 1}`);
    basis.push(row);
  }
  const last = new Array<number>(n + 1).fill(0);
  last[n] = toSafeInteger(target, 'subset-sum target M');
  basis.push(last);
  return basis;
}

/* ------------------------------------------------------------------------ */
/* Exact postcondition for the floating-point LLL engine                     */
/* ------------------------------------------------------------------------ */

function safeIntegerMatrix(matrix: readonly (readonly number[])[]): bigint[][] | null {
  if (matrix.length === 0) return null;
  const width = matrix[0].length;
  if (width === 0 || matrix.some((row) => row.length !== width)) return null;
  const out: bigint[][] = [];
  for (const row of matrix) {
    const converted: bigint[] = [];
    for (const value of row) {
      if (!Number.isSafeInteger(value)) return null;
      converted.push(BigInt(value));
    }
    out.push(converted);
  }
  return out;
}

/** Exact fraction-free determinant over integers (Bareiss elimination). */
function determinantBigInt(matrix: readonly (readonly bigint[])[]): bigint {
  const n = matrix.length;
  if (n === 0 || matrix.some((row) => row.length !== n)) return 0n;
  if (n === 1) return matrix[0][0];
  const work = matrix.map((row) => row.slice());
  let sign = 1n;
  let previousPivot = 1n;
  for (let k = 0; k < n - 1; k += 1) {
    let pivot = k;
    while (pivot < n && work[pivot][k] === 0n) pivot += 1;
    if (pivot === n) return 0n;
    if (pivot !== k) {
      [work[k], work[pivot]] = [work[pivot], work[k]];
      sign = -sign;
    }
    const pivotValue = work[k][k];
    for (let i = k + 1; i < n; i += 1) {
      for (let j = k + 1; j < n; j += 1) {
        work[i][j] =
          (work[i][j] * pivotValue - work[i][k] * work[k][j]) / previousPivot;
      }
      work[i][k] = 0n;
    }
    previousPivot = pivotValue;
  }
  return sign * work[n - 1][n - 1];
}

function exactLatticeTransformError(
  inputBasis: readonly (readonly number[])[],
  reducedBasis: readonly (readonly number[])[],
  transform: readonly (readonly number[])[],
): string | null {
  const input = safeIntegerMatrix(inputBasis);
  const reduced = safeIntegerMatrix(reducedBasis);
  const u = safeIntegerMatrix(transform);
  if (input === null || reduced === null || u === null) {
    return 'a basis or transform entry was not an exact safe integer';
  }
  const n = input.length;
  if (
    input.some((row) => row.length !== n) ||
    reduced.length !== n ||
    reduced.some((row) => row.length !== n) ||
    u.length !== n ||
    u.some((row) => row.length !== n)
  ) {
    return 'the basis or transform was not square with the expected dimension';
  }
  const det = determinantBigInt(u);
  if (det !== 1n && det !== -1n) {
    return `the exact transform determinant was ${det}, not +/-1`;
  }
  for (let row = 0; row < n; row += 1) {
    for (let col = 0; col < n; col += 1) {
      let reconstructed = 0n;
      for (let j = 0; j < n; j += 1) reconstructed += u[row][j] * input[j][col];
      if (reconstructed !== reduced[row][col]) {
        return `reducedBasis != U * inputBasis at row ${row + 1}, column ${col + 1}`;
      }
    }
  }
  return null;
}

/**
 * Independent exact check for the integer-lattice invariant claimed by the UI.
 * This deliberately shares no floating-point dot product with `lllReduce`.
 */
export function verifyExactLatticeTransform(
  inputBasis: readonly (readonly number[])[],
  reducedBasis: readonly (readonly number[])[],
  transform: readonly (readonly number[])[],
): boolean {
  return exactLatticeTransformError(inputBasis, reducedBasis, transform) === null;
}

/**
 * Lagarias-Odlyzko (3.4): e = sum(e_i * b_i) + b_{n+1} is in the lattice, and
 * its last coordinate is M - sum(e_i a_i) = 0.
 *
 * This is asserted on every instance the exhibit builds. An embedding that does
 * not contain the plaintext vector is teaching a broken attack no matter what
 * the reduction happens to return.
 */
export function plaintextVectorInLattice(
  publicWeights: readonly bigint[],
  target: bigint,
  bits: readonly number[],
): number[] {
  assertBits(bits, publicWeights.length);
  const basis = lagariasOdlyzkoBasis(publicWeights, target);
  const n = publicWeights.length;
  const combo = new Array<number>(n + 1).fill(0);
  for (let i = 0; i < n; i += 1) {
    if (bits[i] !== 1) continue;
    for (let c = 0; c <= n; c += 1) combo[c] += basis[i][c];
  }
  for (let c = 0; c <= n; c += 1) combo[c] += basis[n][c];
  return combo;
}

/* ------------------------------------------------------------------------ */
/* Candidate extraction and verification                                      */
/* ------------------------------------------------------------------------ */

/**
 * The only thing that may turn a lattice vector into a claimed plaintext.
 * Exact BigInt equality against the PUBLIC subset sum: sum(bits_i * a_i) == M.
 */
export function verifyCandidate(
  publicWeights: readonly bigint[],
  target: bigint,
  bits: readonly number[],
): boolean {
  if (bits.length !== publicWeights.length) return false;
  let sum = 0n;
  for (let i = 0; i < bits.length; i += 1) {
    if (bits[i] === 1) sum += publicWeights[i];
    else if (bits[i] !== 0) return false;
  }
  return sum === target;
}

/**
 * Algorithm SV step (3): a row is a candidate only if its first n coordinates
 * are all 0 or all one fixed non-zero lambda. Returns the bits, or null with
 * the reason the row was rejected.
 *
 * Exported so the rule the paper states can be tested at the level it is
 * stated. Relaxing it does not usually produce a false recovery -- the exact
 * public-sum check downstream catches nearly everything -- which is precisely
 * why it needs its own test: a row like (1, 0, 2, 0, ...) does not encode the
 * bit vector (1, 0, 1, 0, ...), and reading one as the other is fabricating a
 * plaintext the lattice vector never carried. It would take only a coincidence
 * in the sum for that fabrication to be reported as a break.
 */
export function readCandidateBits(
  row: readonly number[],
  n: number,
): { bits: number[]; lambda: number } | { bits: null; detail: string } {
  const head = row.slice(0, n);
  let lambda = 0;
  for (const v of head) {
    if (v === 0) continue;
    if (lambda === 0) lambda = v;
    else if (v !== lambda) {
      return {
        bits: null,
        detail: `coordinates take at least two distinct non-zero values (${lambda} and ${v})`,
      };
    }
  }
  if (lambda === 0) {
    return { bits: null, detail: 'all leading coordinates are zero' };
  }
  return { bits: head.map((v) => (v === 0 ? 0 : 1)), lambda };
}

/* ------------------------------------------------------------------------ */
/* The attack. PUBLIC DATA ONLY.                                              */
/* ------------------------------------------------------------------------ */

/**
 * Lagarias-Odlyzko Algorithm SV, steps (1)-(4).
 *
 * The parameter is destructured immediately and exactly three properties are
 * read. There is no module-level state holding a private key, a modulus, a
 * multiplier or a plaintext for this function to reach. Ground-truth
 * comparison — "is this the message that was actually sent?" — happens in the
 * caller, AFTER this returns.
 *
 * The search is bounded and disclosed: n+1 rows per pass, two passes, so at
 * most 2(n+1) candidates are ever tested. There is no brute-force fallback of
 * any kind. The licence for looking at every row rather than only the first is
 * Lenstra-Lenstra-Lovasz Proposition (1.12), p. 518, which bounds |b_j| for
 * every j; Proposition (1.11) bounds only |b_1|.
 */
export function attackPublicOnly(instance: PublicKnapsackInstance): KnapsackAttackResult {
  const { publicWeights, target, delta } = instance;

  const empty = {
    bits: null,
    lambda: null,
    pass: null,
    row: null,
    winningRow: null,
    rejections: [] as CandidateRejection[],
    basis: [] as number[][],
    reducedBasis: [] as number[][],
    swapCount: 0,
    exactLatticeTransformVerified: false,
  };

  if (publicWeights.length < 2) {
    return { ...empty, outcome: 'refused', rowsInspected: 0, refusal: 'Need at least two public weights' };
  }
  for (let i = 0; i < publicWeights.length; i += 1) {
    if (publicWeights[i] <= 0n) {
      return {
        ...empty,
        outcome: 'refused',
        rowsInspected: 0,
        refusal: `Algorithm SV requires positive public weights; a_${i + 1} = ${publicWeights[i]}`,
      };
    }
  }
  const total = sumOf(publicWeights);
  // Lagarias-Odlyzko Section 2, immediately before Algorithm SV: the two
  // trivial endpoints are outside its stated domain. The UI rejects the
  // all-zero message before key generation and still permits trapdoor
  // encryption/decryption of the all-one endpoint, but neither endpoint may be
  // dressed up as an execution of Algorithm SV.
  if (target < 1n || target >= total) {
    return {
      ...empty,
      outcome: 'refused',
      rowsInspected: 0,
      refusal: `Algorithm SV requires 1 <= M < sum(a_i); got M = ${target}, sum(a_i) = ${total}`,
    };
  }
  // Lenstra-Lenstra-Lovasz (1.5): the constant may be any fixed 1/4 < d < 1.
  if (!Number.isFinite(delta) || delta <= 0.25 || delta >= 1) {
    return {
      ...empty,
      outcome: 'refused',
      rowsInspected: 0,
      refusal: `delta must satisfy 1/4 < delta < 1 (Lenstra-Lenstra-Lovasz 1982, (1.5)); got ${delta}`,
    };
  }
  const budget = arithmeticBudget(publicWeights, target);
  if (!budget.withinBudget) {
    return {
      ...empty,
      outcome: 'refused',
      rowsInspected: 0,
      refusal:
        `Instance is outside the exact-arithmetic budget: (n+1)*E^2 = ${budget.required.toExponential(3)} ` +
        `exceeds MAX_SAFE_INTEGER = ${budget.limit}. Reducing it could make this Number-based ` +
        'Gram-Schmidt arithmetic inexact before the independent BigInt postcheck, so it is refused rather than attacked.',
    };
  }

  const n = publicWeights.length;
  const rejections: CandidateRejection[] = [];
  let rowsInspected = 0;
  let firstBasis: number[][] = [];
  let firstReduced: number[][] = [];
  let firstSwaps = 0;

  // Step (1)-(3) on M, then step (4): the same on sum(a_i) - M, complemented.
  const passes: { pass: 1 | 2; passTarget: bigint; complement: boolean }[] = [
    { pass: 1, passTarget: target, complement: false },
    { pass: 2, passTarget: total - target, complement: true },
  ];

  for (const { pass, passTarget, complement } of passes) {
    if (passTarget < 0n) continue;
    const basis = lagariasOdlyzkoBasis(publicWeights, passTarget);
    let reduction: ReturnType<typeof lllReduce>;
    try {
      reduction = lllReduce(basis.map((r) => r.slice()), delta);
    } catch (error) {
      const message = error instanceof Error ? error.message : String(error);
      if (!/LLL did not converge within operation bound/.test(message)) throw error;
      return {
        ...empty,
        outcome: 'refused',
        rowsInspected,
        refusal: `Lattice reduction did not complete on pass ${pass}: ${message}`,
      };
    }
    const postconditionError = exactLatticeTransformError(
      basis,
      reduction.reducedBasis,
      reduction.transform,
    );
    if (postconditionError !== null) {
      return {
        ...empty,
        outcome: 'refused',
        rowsInspected,
        refusal:
          `Exact lattice-transform postcheck failed on pass ${pass}: ${postconditionError}. ` +
          'No candidate from this floating-point reduction is accepted.',
      };
    }
    if (pass === 1) {
      firstBasis = basis;
      firstReduced = reduction.reducedBasis.map((r) => r.slice());
      firstSwaps = reduction.swapCount;
    }

    for (let r = 0; r < reduction.reducedBasis.length; r += 1) {
      rowsInspected += 1;
      const row = reduction.reducedBasis[r];
      const read = readCandidateBits(row, n);
      if (read.bits === null) {
        rejections.push({
          pass,
          row: r,
          reason: 'coordinates-outside-bit-set',
          detail: read.detail,
        });
        continue;
      }
      const bits = complement ? read.bits.map((b) => 1 - b) : read.bits;
      // The ONLY thing that may turn a lattice vector into a plaintext claim.
      if (!verifyCandidate(publicWeights, target, bits)) {
        rejections.push({
          pass,
          row: r,
          reason: 'public-sum-mismatch',
          detail: `sum(bits * a) != M for the candidate read off row ${r + 1}`,
        });
        continue;
      }
      return {
        outcome: 'recovered',
        bits,
        lambda: read.lambda,
        pass,
        row: r,
        winningRow: row.slice(),
        rowsInspected,
        rejections,
        basis: firstBasis,
        reducedBasis: firstReduced,
        swapCount: firstSwaps,
        exactLatticeTransformVerified: true,
        refusal: null,
      };
    }
  }

  return {
    outcome: 'no-candidate',
    bits: null,
    lambda: null,
    pass: null,
    row: null,
    winningRow: null,
    rowsInspected,
    rejections,
    basis: firstBasis,
    reducedBasis: firstReduced,
    swapCount: firstSwaps,
    exactLatticeTransformVerified: true,
    refusal: null,
  };
}

/* ------------------------------------------------------------------------ */
/* Verifier-side ground truth. NEVER part of the attack.                      */
/* ------------------------------------------------------------------------ */

/**
 * Every bit vector whose public sum equals the target, by exhaustive search.
 *
 * This exists to answer one question the attack cannot and must not answer:
 * is the instance AMBIGUOUS — does more than one message produce this
 * ciphertext? It is only affordable because n is toy, it is labelled as ground
 * truth wherever it is displayed, and its result is NEVER reported as a
 * recovery. Calling it a break would be exactly the brute-force-disguised-as-
 * an-attack mistake Exhibit 4 already exists to guard against.
 */
export function enumerateSolutions(
  publicWeights: readonly bigint[],
  target: bigint,
  limit = 8,
): { solutions: number[][]; count: number; truncated: boolean } {
  const n = publicWeights.length;
  if (n > MAX_ENUMERATION_BITS) {
    throw new Error(
      `Exhaustive enumeration is refused above ${MAX_ENUMERATION_BITS} weights (asked for ${n})`,
    );
  }
  const solutions: number[][] = [];
  let count = 0;
  const totalMasks = 1 << n;
  for (let mask = 0; mask < totalMasks; mask += 1) {
    let sum = 0n;
    for (let i = 0; i < n; i += 1) if ((mask >> i) & 1) sum += publicWeights[i];
    if (sum !== target) continue;
    count += 1;
    if (solutions.length < limit) {
      solutions.push(Array.from({ length: n }, (_, i) => (mask >> i) & 1));
    }
  }
  return { solutions, count, truncated: count > solutions.length };
}
