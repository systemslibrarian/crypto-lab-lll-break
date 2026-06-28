import { lllReduce } from './lll';
import { dot, gramSchmidt, norm, type Basis, type Vec } from './lattice';
import { randomIntInclusive, randomUnit } from './rng';

interface EnumerationState {
  bestNormSq: number;
  bestCoeffs: number[] | null;
}

export interface BKZBlockImprovement {
  tour: number;
  blockStart: number;
  blockEnd: number;
  headNormBefore: number;
  insertedNorm: number;
}

export interface BKZTourLog {
  tour: number;
  improvementsInTour: number;
  converged: boolean;
}

function centeredMod(x: number, q: number): number {
  const r = ((x % q) + q) % q;
  return r > q / 2 ? r - q : r;
}

function gaussianDiscrete(sigma: number): number {
  const u1 = randomUnit();
  const u2 = randomUnit();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return Math.round(z * sigma);
}

function combine(basis: Basis, coeffs: number[]): Vec {
  const dim = basis[0].length;
  const out = new Array<number>(dim).fill(0);
  for (let i = 0; i < basis.length; i += 1) {
    for (let j = 0; j < dim; j += 1) {
      out[j] += coeffs[i] * basis[i][j];
    }
  }
  return out;
}

export function enumerateSVP(basis: Basis): Vec | null {
  const n = basis.length;
  if (n === 0 || n > 8) {
    return null;
  }

  const gs = gramSchmidt(basis);
  const B = gs.gso.map((v) => dot(v, v));
  const coeffs = new Array<number>(n).fill(0);
  const state: EnumerationState = {
    bestNormSq: Number.POSITIVE_INFINITY,
    bestCoeffs: null,
  };

  for (const v of basis) {
    const ns = dot(v, v);
    if (ns > 0 && ns < state.bestNormSq) {
      state.bestNormSq = ns;
    }
  }

  const recurse = (k: number, partialNormSq: number): void => {
    if (partialNormSq >= state.bestNormSq) {
      return;
    }
    if (k < 0) {
      let nonZero = false;
      for (const c of coeffs) {
        if (c !== 0) {
          nonZero = true;
          break;
        }
      }
      if (!nonZero) {
        return;
      }
      const vec = combine(basis, coeffs);
      const ns = dot(vec, vec);
      if (ns > 0 && ns < state.bestNormSq) {
        state.bestNormSq = ns;
        state.bestCoeffs = coeffs.slice();
      }
      return;
    }

    if (B[k] < 1e-12) {
      return;
    }

    let center = 0;
    for (let j = k + 1; j < n; j += 1) {
      center -= coeffs[j] * gs.mu[j][k];
    }

    const radiusSq = (state.bestNormSq - partialNormSq) / B[k];
    if (radiusSq < 0) {
      return;
    }
    const radius = Math.sqrt(radiusSq);
    const lo = Math.ceil(center - radius);
    const hi = Math.floor(center + radius);

    for (let z = lo; z <= hi; z += 1) {
      coeffs[k] = z;
      const diff = z - center;
      recurse(k - 1, partialNormSq + diff * diff * B[k]);
    }
    coeffs[k] = 0;
  };

  recurse(n - 1, 0);

  if (!state.bestCoeffs) {
    return null;
  }
  return combine(basis, state.bestCoeffs);
}

export function bkzReduce(
  inputBasis: Basis,
  beta: number,
): {
  reducedBasis: Basis;
  tours: number;
  improvements: number;
  blockImprovements: BKZBlockImprovement[];
  tourLogs: BKZTourLog[];
} {
  if (beta < 2) {
    throw new Error('BKZ beta must be >= 2');
  }
  let basis = lllReduce(inputBasis).reducedBasis;
  const n = basis.length;
  let improvements = 0;
  let tours = 0;
  const blockImprovements: BKZBlockImprovement[] = [];
  const tourLogs: BKZTourLog[] = [];

  const maxTours = 20;
  while (tours < maxTours) {
    tours += 1;
    let improvedThisTour = false;
    let improvementsInTour = 0;

    for (let i = 0; i < n; i += 1) {
      const end = Math.min(i + beta, n);
      const block = basis.slice(i, end).map((v) => v.slice());
      if (block.length < 2) {
        continue;
      }
      const short = enumerateSVP(block);
      if (!short) {
        continue;
      }
      const shortNorm = norm(short);
      const headNorm = norm(block[0]);
      if (shortNorm + 1e-9 < headNorm) {
        blockImprovements.push({
          tour: tours,
          blockStart: i,
          blockEnd: end - 1,
          headNormBefore: headNorm,
          insertedNorm: shortNorm,
        });
        // Lattice-preserving insertion (what real BKZ does): insert the
        // enumerated short vector as an extra generator and LLL-reduce the
        // (block.length + 1) set. `short` is an integer combination of the block,
        // so the set is rank-deficient and LLL yields exactly one zero vector;
        // dropping it gives a true basis of the SAME block lattice with `short`
        // pulled to the front. Overwriting block[0] instead (the old shortcut)
        // only preserves the lattice when short's coefficient on block[0] is +-1.
        const generating = [short.slice(), ...block];
        const blockReduced = lllReduce(generating).reducedBasis.filter((v) => dot(v, v) > 1e-9);
        for (let k = 0; k < block.length && k < blockReduced.length; k += 1) {
          basis[i + k] = blockReduced[k];
        }
        improvedThisTour = true;
        improvements += 1;
        improvementsInTour += 1;
      }
    }

    if (!improvedThisTour) {
      tourLogs.push({ tour: tours, improvementsInTour: 0, converged: true });
      break;
    }
    tourLogs.push({ tour: tours, improvementsInTour, converged: false });
    basis = lllReduce(basis).reducedBasis;
  }

  return { reducedBasis: basis, tours, improvements, blockImprovements, tourLogs };
}

/**
 * Primal (Kannan) embedding of an LWE instance b = A*s + e (mod q).
 *
 * The lattice is spanned by the rows of the (n + m + 1) x (n + m + 1) matrix
 *
 *     [ I_n   A^T   0 ]   <- n rows: column i of A^T is column i of A
 *     [ 0     q*I_m 0 ]   <- m rows: the q-ary relations on the error block
 *     [ 0     b^T   1 ]   <- 1 row:  the embedding / target shift
 *
 * Columns are grouped as [ n secret coords | m error coords | 1 embed coord ].
 * Take the last row once and subtract s_i times row i for each i:
 *
 *     (0, b, 1) - sum_i s_i (e_i, A[:,i], 0) = (-s, b - A*s, 1) = (-s, e, 1)
 *
 * after reducing the middle block modulo q with the q*I_m rows (since
 * b - A*s == e (mod q)). That vector has norm ~ sqrt(|s|^2 + |e|^2 + 1), so for
 * toy parameters it is by far the shortest non-zero lattice vector and LLL/BKZ
 * surfaces it. Reading its first n coordinates (negated, since the embed coord
 * comes out as +1) hands back the secret -- a genuine reduction attack, not a
 * search over secrets. This is the construction every other piece of the demo
 * depends on, so it must be the textbook one.
 */
export function buildLWELattice(A: number[][], b: number[], q: number): Basis {
  const m = A.length;
  const n = A[0]?.length ?? 0;
  if (m === 0 || n === 0 || b.length !== m) {
    throw new Error('Invalid LWE dimensions');
  }
  const dim = n + m + 1;
  const basis: Basis = [];

  // n rows: [ I_n | A^T | 0 ]
  for (let i = 0; i < n; i += 1) {
    const row = new Array<number>(dim).fill(0);
    row[i] = 1;
    for (let r = 0; r < m; r += 1) {
      row[n + r] = ((A[r][i] % q) + q) % q;
    }
    basis.push(row);
  }

  // m rows: [ 0 | q*I_m | 0 ]
  for (let r = 0; r < m; r += 1) {
    const row = new Array<number>(dim).fill(0);
    row[n + r] = q;
    basis.push(row);
  }

  // 1 row: [ 0 | b^T | 1 ]
  const last = new Array<number>(dim).fill(0);
  for (let r = 0; r < m; r += 1) {
    last[n + r] = ((b[r] % q) + q) % q;
  }
  last[dim - 1] = 1;
  basis.push(last);

  return basis;
}

export async function generateLWE(
  n: number,
  q: number,
  sigma: number,
): Promise<{ A: number[][]; b: number[]; secret: number[]; errors: number[] }> {
  const m = n + 4;
  const A: number[][] = Array.from({ length: m }, () =>
    Array.from({ length: n }, () => randomIntInclusive(0, q - 1)),
  );
  const secret = Array.from({ length: n }, () => randomIntInclusive(0, 4));
  const errors = Array.from({ length: m }, () => gaussianDiscrete(sigma));

  const outB = new Array<number>(m).fill(0);
  for (let r = 0; r < m; r += 1) {
    let acc = 0;
    for (let c = 0; c < n; c += 1) {
      acc += A[r][c] * secret[c];
    }
    acc += errors[r];
    outB[r] = ((acc % q) + q) % q;
  }

  return { A, b: outB, secret, errors };
}

function checkCandidate(
  A: number[][],
  b: number[],
  q: number,
  candidate: number[],
  sigma = 1,
): { ok: boolean; residualNormSq: number } {
  let residualNormSq = 0;
  for (let r = 0; r < A.length; r += 1) {
    let acc = 0;
    for (let c = 0; c < candidate.length; c += 1) {
      acc += A[r][c] * candidate[c];
    }
    const e = centeredMod(b[r] - acc, q);
    residualNormSq += e * e;
  }
  const avgSq = residualNormSq / Math.max(A.length, 1);
  // Accept when the per-sample residual energy looks like genuine LWE noise.
  // The true secret yields avgSq ~ sigma^2; a wrong secret yields ~ q^2/12 (a
  // huge gap), so a 3-sigma RMS bound (9*sigma^2) cleanly separates them. The
  // old fixed bound of 9 assumed sigma=1 and wrongly rejected the real secret
  // whenever sigma was raised past ~3, reporting success as failure.
  const s = Math.max(sigma, 1);
  return { ok: avgSq <= 9 * s * s, residualNormSq };
}

function bruteForceRecover(
  A: number[][],
  b: number[],
  q: number,
  n: number,
  sigma = 1,
): number[] | null {
  if (n > 8) {
    return null;
  }
  let best: number[] | null = null;
  let bestScore = Number.POSITIVE_INFINITY;
  const cur = new Array<number>(n).fill(0);

  const rec = (idx: number): void => {
    if (idx === n) {
      const check = checkCandidate(A, b, q, cur, sigma);
      if (check.residualNormSq < bestScore) {
        bestScore = check.residualNormSq;
        best = cur.slice();
      }
      return;
    }
    for (let x = 0; x <= 4; x += 1) {
      cur[idx] = x;
      rec(idx + 1);
    }
  };

  rec(0);
  if (!best) {
    return null;
  }
  const valid = checkCandidate(A, b, q, best, sigma);
  return valid.ok ? best : null;
}

function toModQ(x: number, q: number): number {
  return ((x % q) + q) % q;
}

function recoverFromShortVector(
  shortVec: Vec,
  A: number[][],
  b: number[],
  q: number,
  n: number,
  sigma = 1,
): number[] | null {
  const tail = shortVec[shortVec.length - 1];
  if (Math.abs(tail) !== 1) {
    return null;
  }

  const sign = -tail;
  const raw = shortVec.slice(0, n).map((x) => sign * x);
  const candidates: number[][] = [];
  candidates.push(raw.map((x) => toModQ(x, q)));
  candidates.push(raw.map((x) => toModQ(centeredMod(x, q), q)));

  let best: number[] | null = null;
  let bestScore = Number.POSITIVE_INFINITY;
  for (const cand of candidates) {
    const check = checkCandidate(A, b, q, cand, sigma);
    if (check.residualNormSq < bestScore) {
      best = cand;
      bestScore = check.residualNormSq;
    }
    if (check.ok) {
      return cand;
    }
  }

  return best;
}

export interface LWEAttackResult {
  success: boolean;
  recovered: number[] | null;
  shortVec: Vec | null;
  recoverMethod: 'short-vector' | 'bruteforce' | 'none';
}

/**
 * Recover the LWE secret from an ALREADY-REDUCED embedding basis. Kept separate
 * from the reduction step so the UI can run its chosen reduction (BKZ-beta) once
 * and recover from that exact basis -- otherwise the beta slider would be cosmetic,
 * since a hidden re-run of plain LLL would decide success.
 */
export function recoverFromReducedBasis(
  reduced: Basis,
  A: number[][],
  b: number[],
  q: number,
  n: number,
  sigma = 1,
): LWEAttackResult {
  // The genuine lattice attack: among the reduced vectors, the LWE secret rides
  // on the embedding-shaped one (last coordinate +-1, encoding (-s, e, 1)). Scan
  // every candidate shortest-first so a near-shortest match is not missed when
  // the globally shortest vector happens to lie outside the embedding coset.
  const candidates = reduced
    .map((v) => ({ v: v.slice(), ns: dot(v, v) }))
    .filter((c) => c.ns > 0)
    .sort((a, b2) => a.ns - b2.ns);

  // Report the overall shortest vector for the norm-gap meter regardless of
  // whether recovery succeeds, so the UI can show how close reduction got.
  const shortVec: Vec | null = candidates.length > 0 ? candidates[0]!.v.slice() : null;

  for (const c of candidates) {
    if (Math.abs(c.v[c.v.length - 1]) !== 1) {
      continue;
    }
    const recoveredFromShort = recoverFromShortVector(c.v, A, b, q, n, sigma);
    if (recoveredFromShort) {
      const check = checkCandidate(A, b, q, recoveredFromShort, sigma);
      if (check.ok) {
        return { success: true, recovered: recoveredFromShort, shortVec: c.v, recoverMethod: 'short-vector' };
      }
    }
  }

  const recovered = bruteForceRecover(A, b, q, n, sigma);
  if (!recovered) {
    return { success: false, recovered: null, shortVec, recoverMethod: 'none' };
  }

  return { success: true, recovered, shortVec, recoverMethod: 'bruteforce' };
}

/**
 * Convenience end-to-end attack: build the embedding, reduce it with `beta`
 * (default 2 == LLL), then recover. The UI passes the user's chosen beta so a
 * larger block size produces a genuinely stronger reduction.
 */
export function attackLWE(
  A: number[][],
  b: number[],
  q: number,
  n: number,
  sigma = 1,
  beta = 2,
): LWEAttackResult {
  if (n >= 20) {
    return { success: false, recovered: null, shortVec: null, recoverMethod: 'none' };
  }

  const lattice = buildLWELattice(A, b, q);
  const reduced = beta > 2 ? bkzReduce(lattice, beta).reducedBasis : lllReduce(lattice).reducedBasis;
  return recoverFromReducedBasis(reduced, A, b, q, n, sigma);
}
