import { lllReduce } from './lll';
import { dot, gramSchmidt, norm, type Basis, type Vec } from './lattice';

interface EnumerationState {
  bestNormSq: number;
  bestCoeffs: number[] | null;
}

function randomUint32(): number {
  const a = new Uint32Array(1);
  crypto.getRandomValues(a);
  return a[0] ?? 0;
}

function randomIntInclusive(min: number, max: number): number {
  const span = max - min + 1;
  return min + (randomUint32() % span);
}

function randomUnit(): number {
  return (randomUint32() + 1) / 4294967297;
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
): { reducedBasis: Basis; tours: number; improvements: number } {
  if (beta < 2) {
    throw new Error('BKZ beta must be >= 2');
  }
  let basis = lllReduce(inputBasis).reducedBasis;
  const n = basis.length;
  let improvements = 0;
  let tours = 0;

  const maxTours = 20;
  while (tours < maxTours) {
    tours += 1;
    let improvedThisTour = false;

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
        block[0] = short;
        const blockReduced = lllReduce(block).reducedBasis;
        for (let k = 0; k < blockReduced.length; k += 1) {
          basis[i + k] = blockReduced[k];
        }
        improvedThisTour = true;
        improvements += 1;
      }
    }

    if (!improvedThisTour) {
      break;
    }
    basis = lllReduce(basis).reducedBasis;
  }

  return { reducedBasis: basis, tours, improvements };
}

export function buildLWELattice(A: number[][], b: number[], q: number): Basis {
  const m = A.length;
  const n = A[0]?.length ?? 0;
  if (m === 0 || n === 0 || b.length !== m) {
    throw new Error('Invalid LWE dimensions');
  }
  const dim = n + m + 1;
  const basis: Basis = [];

  for (let i = 0; i < n; i += 1) {
    const row = new Array<number>(dim).fill(0);
    row[i] = q;
    basis.push(row);
  }

  for (let r = 0; r < m; r += 1) {
    const row = new Array<number>(dim).fill(0);
    for (let c = 0; c < n; c += 1) {
      row[c] = ((A[r][c] % q) + q) % q;
    }
    row[n + r] = 1;
    basis.push(row);
  }

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
  return { ok: avgSq <= 9, residualNormSq };
}

function bruteForceRecover(A: number[][], b: number[], q: number, n: number): number[] | null {
  if (n > 8) {
    return null;
  }
  let best: number[] | null = null;
  let bestScore = Number.POSITIVE_INFINITY;
  const cur = new Array<number>(n).fill(0);

  const rec = (idx: number): void => {
    if (idx === n) {
      const check = checkCandidate(A, b, q, cur);
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
  const valid = checkCandidate(A, b, q, best);
  return valid.ok ? best : null;
}

export function attackLWE(
  A: number[][],
  b: number[],
  q: number,
  n: number,
): { success: boolean; recovered: number[] | null; shortVec: Vec | null } {
  if (n >= 16) {
    return { success: false, recovered: null, shortVec: null };
  }

  const lattice = buildLWELattice(A, b, q);
  const reduced = lllReduce(lattice).reducedBasis;

  let shortVec: Vec | null = null;
  let shortNormSq = Number.POSITIVE_INFINITY;
  for (const v of reduced) {
    const ns = dot(v, v);
    if (ns > 0 && ns < shortNormSq) {
      shortNormSq = ns;
      shortVec = v.slice();
    }
  }

  const recovered = bruteForceRecover(A, b, q, n);
  if (!recovered) {
    return { success: false, recovered: null, shortVec };
  }

  return { success: true, recovered, shortVec };
}
