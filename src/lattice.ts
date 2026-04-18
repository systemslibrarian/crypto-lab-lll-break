export type Vec = number[];
export type Basis = Vec[];

function assertSameLength(a: Vec, b: Vec): void {
  if (a.length !== b.length) {
    throw new Error('Vector length mismatch');
  }
}

function assertBasis(basis: Basis): void {
  if (basis.length === 0) {
    throw new Error('Basis must be non-empty');
  }
  const dim = basis[0].length;
  if (dim === 0) {
    throw new Error('Vectors must be non-empty');
  }
  for (const v of basis) {
    if (v.length !== dim) {
      throw new Error('All basis vectors must have equal dimension');
    }
  }
}

export function cloneBasis(basis: Basis): Basis {
  return basis.map((v) => v.slice());
}

export function dot(a: Vec, b: Vec): number {
  assertSameLength(a, b);
  let s = 0;
  for (let i = 0; i < a.length; i += 1) {
    s += a[i] * b[i];
  }
  return s;
}

export function norm(v: Vec): number {
  return Math.sqrt(dot(v, v));
}

export function add(a: Vec, b: Vec): Vec {
  assertSameLength(a, b);
  const out = new Array<number>(a.length);
  for (let i = 0; i < a.length; i += 1) {
    out[i] = a[i] + b[i];
  }
  return out;
}

export function sub(a: Vec, b: Vec): Vec {
  assertSameLength(a, b);
  const out = new Array<number>(a.length);
  for (let i = 0; i < a.length; i += 1) {
    out[i] = a[i] - b[i];
  }
  return out;
}

export function scale(v: Vec, c: number): Vec {
  const out = new Array<number>(v.length);
  for (let i = 0; i < v.length; i += 1) {
    out[i] = v[i] * c;
  }
  return out;
}

export function gramSchmidt(basis: Basis): { gso: Vec[]; mu: number[][] } {
  assertBasis(basis);
  const n = basis.length;
  const dim = basis[0].length;

  const gso: Vec[] = Array.from({ length: n }, () => new Array<number>(dim).fill(0));
  const mu: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));

  for (let i = 0; i < n; i += 1) {
    let v = basis[i].slice();
    for (let j = 0; j < i; j += 1) {
      const denom = dot(gso[j], gso[j]);
      if (denom === 0) {
        mu[i][j] = 0;
        continue;
      }
      const coeff = dot(basis[i], gso[j]) / denom;
      mu[i][j] = coeff;
      const proj = scale(gso[j], coeff);
      v = sub(v, proj);
    }
    gso[i] = v;
  }

  return { gso, mu };
}

export function isSizeReduced(mu: number[][], i: number, j: number): boolean {
  return Math.abs(mu[i][j]) <= 0.5;
}

export function lovaszCondition(
  gso: Vec[],
  mu: number[][],
  i: number,
  delta = 0.75,
): boolean {
  if (i < 0 || i + 1 >= gso.length) {
    throw new Error('Invalid Lovasz index');
  }
  const left = delta * dot(gso[i], gso[i]);
  const rhsVec = add(gso[i + 1], scale(gso[i], mu[i + 1][i]));
  const right = dot(rhsVec, rhsVec);
  return left <= right + 1e-12;
}

export function sizeReduce(basis: Basis, i: number, j: number): number {
  const { mu } = gramSchmidt(basis);
  const c = Math.round(mu[i][j]);
  if (c === 0) {
    return 0;
  }
  const bi = basis[i];
  const bj = basis[j];
  for (let k = 0; k < bi.length; k += 1) {
    bi[k] -= c * bj[k];
  }
  return c;
}

export function swapVectors(basis: Basis, k: number): void {
  if (k < 0 || k + 1 >= basis.length) {
    throw new Error('Invalid swap index');
  }
  const tmp = basis[k];
  basis[k] = basis[k + 1];
  basis[k + 1] = tmp;
}
