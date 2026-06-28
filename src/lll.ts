import {
  cloneBasis,
  dot,
  gramSchmidt,
  lovaszCondition,
  norm,
  sizeReduce,
  swapVectors,
  type Basis,
  type Vec,
} from './lattice';

export interface LLLStep {
  type: 'size-reduce' | 'swap' | 'advance' | 'done';
  i: number;
  j?: number;
  before: Basis;
  after: Basis;
  gsoBefore: Vec[];
  gsoAfter: Vec[];
  muBefore: number[][];
  muAfter: number[][];
  lovaszSatisfied?: boolean;
  reductionCoeff?: number;
  // Running unimodular transform U after this step: reducedBasis = U * inputBasis.
  // Every entry is an integer and |det(U)| = 1, which is the proof that the
  // lattice is unchanged while the basis is rewritten.
  transformAfter: number[][];
  description: string;
}

function cloneMu(mu: number[][]): number[][] {
  return mu.map((row) => row.slice());
}

function identity(n: number): number[][] {
  return Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)),
  );
}

// Mirror a basis row operation b_k <- b_k - c*b_j onto the transform matrix, so
// U keeps tracking reducedBasis = U * inputBasis exactly.
function transformRowReduce(u: number[][], k: number, j: number, c: number): void {
  if (c === 0) return;
  for (let col = 0; col < u[k].length; col += 1) {
    u[k][col] -= c * u[j][col];
  }
}

function transformSwap(u: number[][], k: number): void {
  const tmp = u[k];
  u[k] = u[k + 1];
  u[k + 1] = tmp;
}

/** Integer matrix product M * B over rows (each row of M selects rows of B). */
export function applyTransform(u: number[][], basis: Basis): Basis {
  return u.map((urow) => {
    const out = new Array<number>(basis[0].length).fill(0);
    for (let j = 0; j < urow.length; j += 1) {
      const coeff = urow[j];
      if (coeff === 0) continue;
      for (let c = 0; c < out.length; c += 1) {
        out[c] += coeff * basis[j][c];
      }
    }
    return out;
  });
}

export function determinant(matrix: number[][]): number {
  const n = matrix.length;
  const m = matrix.map((row) => row.slice());
  let det = 1;
  for (let i = 0; i < n; i += 1) {
    let pivot = i;
    for (let r = i + 1; r < n; r += 1) {
      if (Math.abs(m[r][i]) > Math.abs(m[pivot][i])) {
        pivot = r;
      }
    }
    if (Math.abs(m[pivot][i]) < 1e-12) {
      return 0;
    }
    if (pivot !== i) {
      const tmp = m[i];
      m[i] = m[pivot];
      m[pivot] = tmp;
      det *= -1;
    }
    const piv = m[i][i];
    det *= piv;
    for (let r = i + 1; r < n; r += 1) {
      const f = m[r][i] / piv;
      for (let c = i; c < n; c += 1) {
        m[r][c] -= f * m[i][c];
      }
    }
  }
  return det;
}

export function lllReduce(
  inputBasis: Basis,
  delta = 0.75,
): { reducedBasis: Basis; steps: LLLStep[]; swapCount: number; transform: number[][] } {
  const basis = cloneBasis(inputBasis);
  const n = basis.length;
  const steps: LLLStep[] = [];
  // U starts as the identity and accumulates every row operation, so that at all
  // times basis === U * inputBasis. See LLLStep.transformAfter.
  const u = identity(n);

  if (n <= 1) {
    const { gso, mu } = gramSchmidt(basis);
    steps.push({
      type: 'done',
      i: 0,
      before: cloneBasis(basis),
      after: cloneBasis(basis),
      gsoBefore: cloneBasis(gso),
      gsoAfter: cloneBasis(gso),
      muBefore: cloneMu(mu),
      muAfter: cloneMu(mu),
      transformAfter: cloneMu(u),
      description: 'LLL complete (trivial basis).',
    });
    return { reducedBasis: basis, steps, swapCount: 0, transform: u };
  }

  let k = 1;
  let swapCount = 0;
  let watchdog = 0;
  const maxOps = 20000;

  while (k < n) {
    watchdog += 1;
    if (watchdog > maxOps) {
      throw new Error('LLL did not converge within operation bound');
    }

    for (let j = k - 1; j >= 0; j -= 1) {
      const before = cloneBasis(basis);
      const beforeGS = gramSchmidt(basis);
      const coeff = sizeReduce(basis, k, j);
      transformRowReduce(u, k, j, coeff);
      const after = cloneBasis(basis);
      const afterGS = gramSchmidt(basis);
      steps.push({
        type: 'size-reduce',
        i: k,
        j,
        before,
        after,
        gsoBefore: cloneBasis(beforeGS.gso),
        gsoAfter: cloneBasis(afterGS.gso),
        muBefore: cloneMu(beforeGS.mu),
        muAfter: cloneMu(afterGS.mu),
        reductionCoeff: coeff,
        transformAfter: cloneMu(u),
        description:
          coeff === 0
            ? `k=${k}, j=${j}: size reduction skipped (round(mu)=0).`
            : `k=${k}, j=${j}: size reduce with coefficient ${coeff}.`,
      });
    }

    const before = cloneBasis(basis);
    const beforeGS = gramSchmidt(basis);
    const ok = lovaszCondition(beforeGS.gso, beforeGS.mu, k - 1, delta);

    if (ok) {
      const after = cloneBasis(basis);
      const afterGS = gramSchmidt(basis);
      steps.push({
        type: 'advance',
        i: k,
        before,
        after,
        gsoBefore: cloneBasis(beforeGS.gso),
        gsoAfter: cloneBasis(afterGS.gso),
        muBefore: cloneMu(beforeGS.mu),
        muAfter: cloneMu(afterGS.mu),
        lovaszSatisfied: true,
        transformAfter: cloneMu(u),
        description: `k=${k}: Lovasz satisfied, advance to k=${k + 1}.`,
      });
      k += 1;
    } else {
      swapVectors(basis, k - 1);
      transformSwap(u, k - 1);
      swapCount += 1;
      const after = cloneBasis(basis);
      const afterGS = gramSchmidt(basis);
      steps.push({
        type: 'swap',
        i: k,
        before,
        after,
        gsoBefore: cloneBasis(beforeGS.gso),
        gsoAfter: cloneBasis(afterGS.gso),
        muBefore: cloneMu(beforeGS.mu),
        muAfter: cloneMu(afterGS.mu),
        lovaszSatisfied: false,
        transformAfter: cloneMu(u),
        description: `k=${k}: Lovasz violated, swap b[${k - 1}] and b[${k}].`,
      });
      k = Math.max(k - 1, 1);
    }
  }

  const finalGS = gramSchmidt(basis);
  steps.push({
    type: 'done',
    i: n - 1,
    before: cloneBasis(basis),
    after: cloneBasis(basis),
    gsoBefore: cloneBasis(finalGS.gso),
    gsoAfter: cloneBasis(finalGS.gso),
    muBefore: cloneMu(finalGS.mu),
    muAfter: cloneMu(finalGS.mu),
    transformAfter: cloneMu(u),
    description: `LLL complete: ${swapCount} swaps, ${steps.length} trace steps.`,
  });

  return { reducedBasis: basis, steps, swapCount, transform: u };
}

export function basisProfile(basis: Basis): number[] {
  const { gso } = gramSchmidt(basis);
  return gso.map((v) => {
    const n = norm(v);
    return n > 0 ? Math.log2(n) : Number.NEGATIVE_INFINITY;
  });
}

export function orthogonalityDefect(basis: Basis): number {
  if (basis.length === 0) {
    return 1;
  }
  const dim = basis[0].length;
  if (basis.length !== dim) {
    throw new Error('Orthogonality defect requires square basis matrix');
  }
  let prod = 1;
  for (const v of basis) {
    prod *= norm(v);
  }
  const detAbs = Math.abs(determinant(basis));
  if (detAbs < 1e-12) {
    return Number.POSITIVE_INFINITY;
  }
  return prod / detAbs;
}

export function shortestVectorInBasis(basis: Basis): Vec {
  let best = basis[0];
  let bestN = dot(best, best);
  for (let i = 1; i < basis.length; i += 1) {
    const candN = dot(basis[i], basis[i]);
    if (candN < bestN) {
      best = basis[i];
      bestN = candN;
    }
  }
  return best.slice();
}
