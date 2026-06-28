// Shared randomness source. Two modes:
//   - crypto (default): cryptographically strong, non-reproducible. Best for free
//     exploration so learners never see a fake/predictable "random" basis.
//   - seeded: a deterministic mulberry32 PRNG, so a teacher can hand every student
//     the same instance via a seed in a shareable lab link.
// Neither path uses Math.random (the repo guardrail forbids it), and the seeded
// PRNG is pure arithmetic so it is portable and reproducible across machines.

export interface Rng {
  nextUint32(): number;
}

class CryptoRng implements Rng {
  nextUint32(): number {
    const a = new Uint32Array(1);
    crypto.getRandomValues(a);
    return a[0] ?? 0;
  }
}

class SeededRng implements Rng {
  private state: number;

  constructor(seed: number) {
    // Spread small/zero seeds across the state space so seed=0 still mixes well.
    this.state = (seed ^ 0x9e3779b9) >>> 0;
  }

  nextUint32(): number {
    // mulberry32
    this.state = (this.state + 0x6d2b79f5) | 0;
    let t = this.state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return (t ^ (t >>> 14)) >>> 0;
  }
}

let current: Rng = new CryptoRng();
let activeSeed: number | null = null;

/** Switch RNG mode. Pass null for crypto randomness, or an integer seed for a
 *  reproducible deterministic stream. Resets the seeded stream from the start. */
export function setSeed(seed: number | null): void {
  if (seed === null || !Number.isFinite(seed)) {
    current = new CryptoRng();
    activeSeed = null;
  } else {
    const s = Math.trunc(seed) >>> 0;
    current = new SeededRng(s);
    activeSeed = s;
  }
}

/** The active seed, or null when in crypto mode. */
export function getSeed(): number | null {
  return activeSeed;
}

export function nextUint32(): number {
  return current.nextUint32();
}

export function randomIntInclusive(min: number, max: number): number {
  const span = max - min + 1;
  return min + (nextUint32() % span);
}

/** Uniform in the open interval (0, 1) -- never 0, so log(u) is safe (Box-Muller). */
export function randomUnit(): number {
  return (nextUint32() + 1) / 4294967297;
}
