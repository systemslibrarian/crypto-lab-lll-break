/**
 * E-K-09 algorithmic verification for the Merkle-Hellman exhibit.
 *
 * Every check has a NAME, and `--only "<name>"` runs exactly one of them. That
 * is what lets `scripts/knapsack-mutations.mjs` point a single named check at a
 * single mutation and show that this check, and not the suite in general, is
 * what catches it.
 *
 * Failures are raised through `node:assert/strict`, so a caught mutation fails
 * with ERR_ASSERTION rather than a crash. The mutation harness rejects anything
 * else as an infrastructure failure rather than counting it as a kill.
 *
 * Known-answer values are PINNED LITERALS. The 1978 paper's own worked example
 * is the anchor: its five public weights, its transformed sum and its recovered
 * bits are copied from page 526 and are not computed by this repository.
 */
import assert from 'node:assert/strict';
import {
  arithmeticBudget,
  attackPublicOnly,
  decryptWithTrapdoor,
  density,
  derivePublicWeights,
  encryptBits,
  enumerateSolutions,
  firstNonSuperincreasingIndex,
  generateKeyPair,
  isSuperincreasing,
  lagariasOdlyzkoBasis,
  modInverse,
  plaintextVectorInLattice,
  readCandidateBits,
  sumOf,
  verifyCandidate,
  verifyExactLatticeTransform,
  type MerkleHellmanPrivateKey,
  type PublicKnapsackInstance,
} from '../src/knapsack';
import { setSeed } from '../src/rng';

/**
 * Capture a thrown error rather than using `assert.throws(fn, /re/)`.
 *
 * `assert.throws` with a pattern RE-THROWS the actual error when it does not
 * match, so a mutation caught that way surfaces as the code's own Error and not
 * as an ERR_ASSERTION. The mutation harness insists on ERR_ASSERTION precisely
 * to tell a caught mutation from a crashed one, so every negative check below
 * goes through here and asserts on the captured message.
 */
function expectThrow(fn: () => unknown): Error {
  try {
    fn();
  } catch (error) {
    return error as Error;
  }
  assert.fail('expected a throw, but the call returned normally');
}

const checks: { name: string; run: () => void }[] = [];
function check(name: string, run: () => void): void {
  checks.push({ name, run });
}

/* ------------------------------------------------------------------ */
/* Pinned fixtures                                                      */
/* ------------------------------------------------------------------ */

/**
 * Merkle and Hellman 1978, p. 526, the paper's own n = 5 example.
 * Every number below is copied from the printed page.
 */
const PAPER = {
  privateWeights: [171n, 196n, 457n, 1191n, 2410n],
  modulus: 8443n,
  multiplier: 2550n,
  multiplierInverse: 3950n,
  publicWeights: [5457n, 1663n, 216n, 6013n, 7439n],
  ciphertext: 15115n,
  transformedSum: 3797n,
  // The paper states x5=1, x4=1, x3=0, x2=1, x1=0.
  bits: [0, 1, 0, 1, 1],
};

const paperKey: MerkleHellmanPrivateKey = {
  privateWeights: PAPER.privateWeights,
  modulus: PAPER.modulus,
  multiplier: PAPER.multiplier,
  multiplierInverse: PAPER.multiplierInverse,
};

/** Seeded instances, pinned as literals so a drift in key generation shows. */
const PINNED = [
  {
    label: 'n=8, k=14, seed 1234',
    n: 8,
    scaleBits: 14,
    seed: 1234,
    bits: '01101001',
    publicWeights: [672190n, 3094370n, 1937532n, 5012433n, 5607925n, 4539592n, 6072687n, 195388n],
    target: 10835215n,
    density: 0.355021,
  },
  {
    label: 'n=10, k=9, seed 1234',
    n: 10,
    scaleBits: 9,
    seed: 1234,
    bits: '0110100101',
    publicWeights: [510405n, 515269n, 222428n, 219216n, 430876n, 155128n, 420384n, 245075n, 258843n, 788142n],
    target: 2201790n,
    density: 0.510514,
  },
  {
    label: 'n=12, k=6, seed 1234',
    n: 12,
    scaleBits: 6,
    seed: 1234,
    bits: '011010010110',
    publicWeights: [332133n, 321978n, 336601n, 355692n, 342084n, 211288n, 126459n, 167346n, 56583n, 36598n, 396325n, 100221n],
    target: 1600932n,
    density: 0.645289,
  },
];

const DELTA = 0.99;
const toBits = (s: string): number[] => s.split('').map((c) => (c === '1' ? 1 : 0));

function seededKey(n: number, scaleBits: number, seed: number) {
  setSeed(seed);
  try {
    return generateKeyPair(n, scaleBits);
  } finally {
    setSeed(null);
  }
}

/* ------------------------------------------------------------------ */
/* The cryptosystem                                                     */
/* ------------------------------------------------------------------ */

check("the 1978 paper's own worked example reproduces exactly", () => {
  const derived = derivePublicWeights(paperKey);
  assert.deepEqual(
    derived.map(String),
    PAPER.publicWeights.map(String),
    'public weights must match the five printed on Merkle-Hellman p. 526',
  );
  assert.equal(
    encryptBits(PAPER.publicWeights, PAPER.bits).toString(),
    PAPER.ciphertext.toString(),
    'the paper sums a_2 + a_4 + a_5 = 15115',
  );
  const decrypted = decryptWithTrapdoor(paperKey, PAPER.ciphertext);
  assert.equal(
    decrypted.transformedSum.toString(),
    PAPER.transformedSum.toString(),
    "S' = w^-1 * S mod m must be the printed 3797",
  );
  assert.deepEqual(decrypted.bits, PAPER.bits, 'the paper prints x = (0,1,0,1,1)');
  assert.equal(decrypted.remainder.toString(), '0', 'the greedy pass must account for all of S-prime');
});

check("the greedy pass uses >= and not >, as the paper's own example requires", () => {
  // At step i=2 of the paper's example the residue is exactly a'_2 = 196, and
  // the printed answer is x_2 = 1. A strict `>` yields x_2 = 0, a different
  // message, and a non-zero remainder. This check is what pins the glyph.
  const decrypted = decryptWithTrapdoor(paperKey, PAPER.ciphertext);
  assert.equal(decrypted.bits[1], 1, 'x_2 must be 1 when the residue equals a_2 exactly');
  assert.equal(decrypted.remainder.toString(), '0');

  // And directly: a one-weight knapsack whose target equals the weight.
  const exact: MerkleHellmanPrivateKey = {
    privateWeights: [5n, 11n],
    modulus: 17n,
    multiplier: 3n,
    multiplierInverse: modInverse(3n, 17n)!,
  };
  const pub = derivePublicWeights(exact);
  const got = decryptWithTrapdoor(exact, encryptBits(pub, [0, 1]));
  assert.deepEqual(got.bits, [0, 1], 'a residue equal to a weight must set that bit');
});

check('a private sequence that is not superincreasing is refused', () => {
  assert.equal(isSuperincreasing([1n, 2n, 4n, 8n]), true);
  // 3 does not exceed 1 + 2.
  assert.equal(isSuperincreasing([1n, 2n, 3n]), false);
  assert.equal(firstNonSuperincreasingIndex([1n, 2n, 3n]), 2);
  // Equality is a violation too: (1) is strict.
  assert.equal(isSuperincreasing([1n, 1n]), false);

  assert.match(
    expectThrow(() =>
      validate({
        privateWeights: [1n, 2n, 3n],
        modulus: 101n,
        multiplier: 7n,
        multiplierInverse: modInverse(7n, 101n)!,
      }),
    ).message,
    /not superincreasing/,
    'a non-superincreasing private sequence must be refused by name',
  );
});

check('a modulus that does not exceed the private-weight sum is refused', () => {
  const weights = [2n, 5n, 11n];
  const total = sumOf(weights);
  assert.equal(total.toString(), '18');
  // m == sum is a violation: (8) is strict.
  assert.match(
    expectThrow(() =>
      validate({ privateWeights: weights, modulus: total, multiplier: 7n, multiplierInverse: 0n }),
    ).message,
    /must exceed the private-weight sum/,
    'm equal to the sum is a violation: (8) is strict',
  );
  assert.match(
    expectThrow(() =>
      validate({ privateWeights: weights, modulus: 17n, multiplier: 7n, multiplierInverse: 0n }),
    ).message,
    /must exceed the private-weight sum/,
  );
  // One above the sum is legal.
  const ok: MerkleHellmanPrivateKey = {
    privateWeights: weights,
    modulus: 19n,
    multiplier: 7n,
    multiplierInverse: modInverse(7n, 19n)!,
  };
  assert.doesNotThrow(() => validate(ok));
});

check('a multiplier that is not invertible modulo m is refused', () => {
  assert.equal(modInverse(6n, 9n), null, 'gcd(6,9) = 3, so no inverse exists');
  assert.equal(modInverse(7n, 19n)?.toString(), '11');
  assert.equal(((7n * 11n) % 19n).toString(), '1');
  assert.match(
    expectThrow(() =>
      validate({
        privateWeights: [2n, 5n, 11n],
        modulus: 18n + 27n,
        multiplier: 15n,
        multiplierInverse: 1n,
      }),
    ).message,
    /not invertible modulo/,
    'gcd(15, 45) = 15, so (3) has no inverse and the key must be refused',
  );
});

check('public weights are w*a_prime mod m at every index', () => {
  for (const fixture of PINNED) {
    const pair = seededKey(fixture.n, fixture.scaleBits, fixture.seed);
    const priv = pair.privateKey;
    for (let i = 0; i < fixture.n; i += 1) {
      assert.equal(
        pair.publicWeights[i].toString(),
        ((priv.multiplier * priv.privateWeights[i]) % priv.modulus).toString(),
        `${fixture.label}: public weight ${i + 1} must be w*a'_${i + 1} mod m`,
      );
    }
    assert.deepEqual(
      pair.publicWeights.map(String),
      fixture.publicWeights.map(String),
      `${fixture.label}: seeded key generation drifted from its pinned public weights`,
    );
  }
});

check('trapdoor decryption inverts encryption on many generated keys', () => {
  let tested = 0;
  for (let seed = 1; seed <= 40; seed += 1) {
    const pair = seededKey(8, 12, seed);
    assert.equal(
      isSuperincreasing(pair.privateKey.privateWeights),
      true,
      `seed ${seed}: generated private sequence must satisfy (1)`,
    );
    assert.equal(
      pair.privateKey.modulus > sumOf(pair.privateKey.privateWeights),
      true,
      `seed ${seed}: generated modulus must satisfy (8)`,
    );
    for (const pattern of ['01101001', '11111111', '10000000', '00000001']) {
      const bits = toBits(pattern);
      const target = encryptBits(pair.publicWeights, bits);
      const got = decryptWithTrapdoor(pair.privateKey, target);
      assert.deepEqual(got.bits, bits, `seed ${seed}, message ${pattern}: round trip failed`);
      assert.equal(got.remainder.toString(), '0');
      tested += 1;
    }
  }
  assert.equal(tested, 160, 'the sweep must actually have run');
});

check('reading the message bits in the reverse order is a different message', () => {
  const pair = seededKey(8, 12, 5);
  const forward = toBits('01101001');
  const reversed = [...forward].reverse();
  assert.notDeepEqual(forward, reversed, 'the fixture must be order-sensitive to mean anything');
  const target = encryptBits(pair.publicWeights, forward);
  assert.notEqual(
    encryptBits(pair.publicWeights, reversed).toString(),
    target.toString(),
    'bit order is load-bearing: b_i multiplies a_i, ascending index',
  );
  assert.deepEqual(
    decryptWithTrapdoor(pair.privateKey, target).bits,
    forward,
    'the trapdoor returns the ascending-index reading, not the reversed one',
  );
});

/* ------------------------------------------------------------------ */
/* The lattice                                                          */
/* ------------------------------------------------------------------ */

check('the Lagarias-Odlyzko basis is the printed (2.4) and contains the plaintext vector', () => {
  const pub = [5457n, 1663n, 216n, 6013n, 7439n];
  const basis = lagariasOdlyzkoBasis(pub, 15115n);
  assert.equal(basis.length, 6, 'the basis is (n+1)-dimensional');
  assert.deepEqual(basis[0], [1, 0, 0, 0, 0, -5457], 'b_1 = (1,0,...,0,-a_1)');
  assert.deepEqual(basis[4], [0, 0, 0, 0, 1, -7439], 'b_n = (0,...,1,-a_n)');
  assert.deepEqual(basis[5], [0, 0, 0, 0, 0, 15115], 'b_{n+1} = (0,...,0,M)');

  // Lagarias-Odlyzko (3.4). An embedding that does not contain the plaintext
  // vector is teaching a broken attack whatever the reduction returns.
  for (const fixture of [...PINNED]) {
    const pair = seededKey(fixture.n, fixture.scaleBits, fixture.seed);
    const bits = toBits(fixture.bits);
    const target = encryptBits(pair.publicWeights, bits);
    const embedded = plaintextVectorInLattice(pair.publicWeights, target, bits);
    assert.deepEqual(
      embedded.slice(0, fixture.n),
      bits,
      `${fixture.label}: the first n coordinates must be the plaintext bits`,
    );
    assert.equal(
      embedded[fixture.n],
      0,
      `${fixture.label}: the last coordinate must be M - sum(e_i a_i) = 0`,
    );
  }
});

check('the reduced basis is postchecked as an exact unimodular transform', () => {
  const identity = [
    [1, 0],
    [0, 1],
  ];
  const unimodular = [
    [1, 2],
    [0, 1],
  ];
  assert.equal(
    verifyExactLatticeTransform(identity, unimodular, unimodular),
    true,
    'U * I = U with det(U)=1 must pass',
  );
  assert.equal(
    verifyExactLatticeTransform(identity, [[2, 0], [0, 1]], [[2, 0], [0, 1]]),
    false,
    'a transform with determinant 2 changes the lattice and must fail',
  );
  assert.equal(
    verifyExactLatticeTransform(identity, [[1, 1], [0, 1]], identity),
    false,
    'a reduced basis that is not U times the input must fail',
  );
  assert.equal(
    verifyExactLatticeTransform(identity, [[Number.MAX_SAFE_INTEGER + 1, 0], [0, 1]], identity),
    false,
    'an inexact Number entry must fail before BigInt conversion',
  );

  for (const fixture of PINNED) {
    const result = attackPublicOnly({
      publicWeights: fixture.publicWeights,
      target: fixture.target,
      delta: DELTA,
    });
    assert.notEqual(result.outcome, 'refused', `${fixture.label}: live reduction unexpectedly refused`);
    assert.equal(
      result.exactLatticeTransformVerified,
      true,
      `${fixture.label}: every accepted reduction must carry the exact postcondition`,
    );
  }
});

check('density follows n over log2 of the largest public weight', () => {
  // Pinned by hand: max weight 1024 = 2^10, so d = 8/10 exactly.
  assert.equal(density([3n, 1024n, 77n, 5n, 9n, 11n, 13n, 17n]), 0.8);
  assert.equal(density([2n, 4n, 65536n, 8n]), 4 / 16);
  for (const fixture of PINNED) {
    const d = density(fixture.publicWeights);
    assert.equal(
      Math.abs(d - fixture.density) < 1e-6,
      true,
      `${fixture.label}: density ${d} drifted from its pinned ${fixture.density}`,
    );
  }
});

check('an instance outside the exact-arithmetic budget is refused, not attacked', () => {
  // This is inside the page's own controls. Measuring pass 1 alone approves it,
  // but Algorithm SV step (4) has a much larger target and must make the whole
  // two-pass run refuse.
  const pair = seededKey(12, 10, 1);
  const target = encryptBits(pair.publicWeights, toBits('100000000000'));
  const complementTarget = sumOf(pair.publicWeights) - target;
  let passOneMax = target;
  for (const a of pair.publicWeights) if (a > passOneMax) passOneMax = a;
  const passOneOnlyRequired = 13 * Number(passOneMax) ** 2;
  assert.equal(
    passOneOnlyRequired < Number.MAX_SAFE_INTEGER,
    true,
    'the regression fixture must expose the old pass-1-only budget bug',
  );
  assert.equal(complementTarget.toString(), '47730761');
  const budget = arithmeticBudget(pair.publicWeights, target);
  assert.equal(budget.withinBudget, false, 'this fixture exists to be over budget');
  assert.equal(budget.maxEntry, Number(complementTarget), 'E must include the pass-2 target');
  assert.equal(
    budget.required > Number.MAX_SAFE_INTEGER,
    true,
    'the complete two-pass input screen must exceed MAX_SAFE_INTEGER',
  );
  const result = attackPublicOnly({ publicWeights: pair.publicWeights, target, delta: DELTA });
  assert.equal(result.outcome, 'refused', 'an over-budget instance must be refused');
  assert.match(String(result.refusal), /exact-arithmetic budget/);
  assert.equal(result.bits, null, 'a refusal must never carry bits');

  // And the in-budget fixtures really are in budget, so the guard is not just
  // refusing everything.
  for (const fixture of PINNED) {
    assert.equal(
      arithmeticBudget(fixture.publicWeights, fixture.target).withinBudget,
      true,
      `${fixture.label} must be inside the budget`,
    );
  }
});

check('Algorithm SV refuses weights and targets outside its published domain', () => {
  const weights = [3n, 7n, 12n];
  for (const target of [-1n, 0n, 22n, 23n]) {
    const result = attackPublicOnly({ publicWeights: weights, target, delta: DELTA });
    assert.equal(result.outcome, 'refused', `target ${target} must be refused`);
    assert.match(String(result.refusal), /1 <= M < sum\(a_i\)/);
    assert.equal(result.rowsInspected, 0, 'an out-of-domain target must never reach LLL');
  }
  for (const bad of [[0n, 7n, 12n], [-3n, 7n, 12n]]) {
    const result = attackPublicOnly({ publicWeights: bad, target: 7n, delta: DELTA });
    assert.equal(result.outcome, 'refused');
    assert.match(String(result.refusal), /positive public weights/);
    assert.equal(result.rowsInspected, 0);
  }
});

check('delta outside the Lenstra-Lenstra-Lovasz range is refused', () => {
  const f = PINNED[0];
  for (const delta of [0.25, 0.1, 1, 1.5, Number.NaN]) {
    const result = attackPublicOnly({
      publicWeights: f.publicWeights,
      target: f.target,
      delta,
    });
    assert.equal(result.outcome, 'refused', `delta ${delta} must be refused`);
    assert.match(String(result.refusal), /1\/4 < delta < 1/);
  }
  assert.equal(
    attackPublicOnly({ publicWeights: f.publicWeights, target: f.target, delta: 0.75 }).outcome,
    'recovered',
    'a legal delta must still work',
  );
});

/* ------------------------------------------------------------------ */
/* The attack                                                           */
/* ------------------------------------------------------------------ */

check('the public-only attack recovers the message on pinned low-density instances', () => {
  for (const fixture of PINNED) {
    const bits = toBits(fixture.bits);
    const result = attackPublicOnly({
      publicWeights: fixture.publicWeights,
      target: fixture.target,
      delta: DELTA,
    });
    assert.equal(result.outcome, 'recovered', `${fixture.label}: expected a recovery`);
    assert.deepEqual(result.bits, bits, `${fixture.label}: recovered the wrong bits`);
    const independentSum = result.bits!.reduce(
      (sum, bit, i) => sum + (bit === 1 ? fixture.publicWeights[i] : 0n),
      0n,
    );
    assert.equal(
      independentSum.toString(),
      fixture.target.toString(),
      `${fixture.label}: an independent sum must verify the recovered bits`,
    );
    assert.equal(result.exactLatticeTransformVerified, true);
    // The search is bounded and that bound is asserted, not just documented.
    assert.equal(
      result.rowsInspected <= 2 * (fixture.n + 1),
      true,
      `${fixture.label}: inspected ${result.rowsInspected} rows, over the 2(n+1) bound`,
    );
  }
});

check('Algorithm SV step (4) recovers a fixture that needs the complement pass', () => {
  const bits = toBits('0110100101100101');
  const pair = seededKey(16, 1, 5);
  const target = encryptBits(pair.publicWeights, bits);
  const result = attackPublicOnly({ publicWeights: pair.publicWeights, target, delta: DELTA });
  assert.equal(result.outcome, 'recovered');
  assert.equal(result.pass, 2, 'this fixture must cross the complement-pass boundary');
  assert.equal(result.rowsInspected > 17, true, 'pass 1 must have been exhausted first');
  assert.deepEqual(result.bits, bits, 'complementing the pass-2 candidate must recover the sent bits');
  const independentSum = result.bits!.reduce(
    (sum, bit, i) => sum + (bit === 1 ? pair.publicWeights[i] : 0n),
    0n,
  );
  assert.equal(independentSum.toString(), target.toString());
  assert.equal(result.exactLatticeTransformVerified, true);
});

check('the attack recovers with no private key in existence at all', () => {
  // Public literals only. No key pair object is constructed anywhere in this
  // check, so there is nothing private for the attack to have read.
  const publicWeights = [672190n, 3094370n, 1937532n, 5012433n, 5607925n, 4539592n, 6072687n, 195388n];
  const target = 10835215n;
  const result = attackPublicOnly({ publicWeights, target, delta: DELTA });
  assert.equal(result.outcome, 'recovered');
  assert.deepEqual(result.bits, toBits('01101001'));
});

check('the attack reads exactly the three declared properties of its input', () => {
  const f = PINNED[0];
  const allowed = new Set(['publicWeights', 'target', 'delta']);
  const touched: string[] = [];
  const backing = { publicWeights: f.publicWeights, target: f.target, delta: DELTA };
  const spy = new Proxy(backing as unknown as PublicKnapsackInstance, {
    get(obj, prop) {
      const key = String(prop);
      touched.push(key);
      if (!allowed.has(key)) {
        throw new Error(`attackPublicOnly read an undeclared property: ${key}`);
      }
      return (obj as Record<string, unknown>)[key];
    },
  });
  const result = attackPublicOnly(spy);
  assert.equal(result.outcome, 'recovered', 'the attack must still work through the proxy');
  assert.deepEqual(
    [...new Set(touched)].sort(),
    ['delta', 'publicWeights', 'target'],
    'the attack must read exactly the three declared properties and nothing else',
  );
});

check('a candidate whose public sum differs from the target is rejected', () => {
  const f = PINNED[0];
  const right = toBits(f.bits);
  assert.equal(verifyCandidate(f.publicWeights, f.target, right), true);

  // Flip one bit: still a legal bit vector, no longer the right sum.
  for (let i = 0; i < right.length; i += 1) {
    const wrong = [...right];
    wrong[i] = wrong[i] === 1 ? 0 : 1;
    assert.equal(
      verifyCandidate(f.publicWeights, f.target, wrong),
      false,
      `flipping bit ${i + 1} must fail the public-sum check`,
    );
  }
  // Wrong length is not a candidate either.
  assert.equal(verifyCandidate(f.publicWeights, f.target, right.slice(1)), false);

  // And a real attack on a target nobody encrypted must not report a recovery.
  const bogus = attackPublicOnly({
    publicWeights: f.publicWeights,
    target: f.target + 1n,
    delta: DELTA,
  });
  assert.notEqual(bogus.outcome, 'refused', 'the bogus target is still inside the budget');
  if (bogus.outcome === 'recovered') {
    assert.equal(
      verifyCandidate(f.publicWeights, f.target + 1n, bogus.bits!),
      true,
      'any reported recovery must satisfy the sum it was asked about',
    );
  }
});

check('candidate coordinates outside the permitted bit set are rejected', () => {
  const f = PINNED[0];
  // Values other than 0 and 1 are not bits, whatever a lattice row holds.
  assert.equal(verifyCandidate(f.publicWeights, f.target, [0, 1, 1, 2, 1, 0, 0, 1]), false);
  assert.equal(verifyCandidate(f.publicWeights, f.target, [0, 1, 1, -1, 1, 0, 0, 1]), false);

  // Algorithm SV step (3) at the level the paper states it. A row is a
  // candidate only when its leading coordinates are all 0 or all ONE fixed
  // non-zero lambda. These four cases are what that sentence means.
  const one = readCandidateBits([1, 0, 1, 0, -5], 4);
  assert.deepEqual(one.bits, [1, 0, 1, 0], 'lambda = 1 is the ordinary case');
  assert.equal((one as { lambda: number }).lambda, 1);

  const scaled = readCandidateBits([2, 0, 2, 0, -5], 4);
  assert.deepEqual(scaled.bits, [1, 0, 1, 0], 'a fixed lambda = 2 is still legal');
  assert.equal((scaled as { lambda: number }).lambda, 2);

  const negated = readCandidateBits([-1, 0, -1, 0, 5], 4);
  assert.deepEqual(negated.bits, [1, 0, 1, 0], 'lambda = -1 covers the negated vector');

  // The one that matters: TWO distinct non-zero values is not a bit vector,
  // and must be refused rather than flattened into one.
  const mixed = readCandidateBits([1, 0, 2, 0, -5], 4);
  assert.equal(mixed.bits, null, 'a row with two distinct non-zero values is not a candidate');
  assert.match(
    String((mixed as { detail: string }).detail),
    /two distinct non-zero values/,
    'and the refusal must say why',
  );

  const zero = readCandidateBits([0, 0, 0, 0, 7], 4);
  assert.equal(zero.bits, null, 'an all-zero head carries no message');

  // And the rule is live in the pipeline, not merely present in a helper: a
  // real reduction produces rows refused for exactly this reason.
  const sweepRejects = attackPublicOnly({
    publicWeights: PINNED[2].publicWeights,
    target: PINNED[2].target + 7n,
    delta: DELTA,
  }).rejections.filter((r) => r.reason === 'coordinates-outside-bit-set');
  assert.equal(
    sweepRejects.length > 0,
    true,
    'a real reduction must produce rows refused for coordinates outside the bit set',
  );

  // Every row the attack actually accepted had first-n coordinates drawn from
  // {0, lambda} for ONE lambda. Assert that on the winning row itself.
  const result = attackPublicOnly({
    publicWeights: f.publicWeights,
    target: f.target,
    delta: DELTA,
  });
  assert.equal(result.outcome, 'recovered');
  const head = result.winningRow!.slice(0, f.n);
  const lambda = result.lambda!;
  assert.notEqual(lambda, 0, 'lambda must be non-zero');
  for (const v of head) {
    assert.equal(
      v === 0 || v === lambda,
      true,
      `winning row coordinate ${v} is neither 0 nor lambda=${lambda}`,
    );
  }
  // And at least one row was rejected somewhere in the sweep, so the rejection
  // path is exercised rather than merely present.
  const sweep = attackPublicOnly({
    publicWeights: PINNED[2].publicWeights,
    target: PINNED[2].target + 7n,
    delta: DELTA,
  });
  assert.equal(
    sweep.rejections.length > 0,
    true,
    'a target nobody encrypted must produce rejected rows',
  );
});

check('a bounded higher-density sweep shows honest failures and claims no threshold', () => {
  // Seeds 1..80 at ONE dimension, with one fixed message. Every figure printed
  // in the learner-visible hint is owned here. The first draft tested only the
  // k=1 half of its sentence; its untested "none at k=3" half was false.
  const bits = toBits('0110100101100101');
  const expected = new Map<number, number[]>([
    [1, [30, 37, 46, 64, 71, 73, 78]],
    [3, [24, 27]],
    [5, []],
  ]);
  for (const [scaleBits, expectedFailures] of expected) {
    const failures: number[] = [];
    let attempted = 0;
    for (let seed = 1; seed <= 80; seed += 1) {
      const pair = seededKey(16, scaleBits, seed);
      const target = encryptBits(pair.publicWeights, bits);
      if (!arithmeticBudget(pair.publicWeights, target).withinBudget) continue;
      attempted += 1;
      const result = attackPublicOnly({ publicWeights: pair.publicWeights, target, delta: DELTA });
      if (result.outcome !== 'recovered' || result.bits!.join('') !== bits.join('')) {
        failures.push(seed);
      }
    }
    assert.equal(attempted, 80, `k=${scaleBits}: every sweep instance must be inside the budget`);
    assert.deepEqual(
      failures,
      expectedFailures,
      `k=${scaleBits}: the deterministic failure set drifted`,
    );
  }
});

check('an in-budget n=8 failure prevents dimension or density generalisations', () => {
  const bits = toBits('00001100');
  const pair = seededKey(8, 1, 10);
  const target = encryptBits(pair.publicWeights, bits);
  assert.equal(
    arithmeticBudget(pair.publicWeights, target).withinBudget,
    true,
    'the failure must not be an arithmetic-screen refusal',
  );
  assert.deepEqual(
    decryptWithTrapdoor(pair.privateKey, target).bits,
    bits,
    'the trapdoor must still decrypt the instance exactly',
  );
  const result = attackPublicOnly({ publicWeights: pair.publicWeights, target, delta: DELTA });
  assert.equal(result.outcome, 'no-candidate');
  assert.equal(result.bits, null);
  assert.equal(result.exactLatticeTransformVerified, true);
  assert.equal(result.rowsInspected, 2 * 9, 'both nine-row passes must be exhausted');
});

check('the deterministic failure fixture fails cleanly and claims nothing', () => {
  const bits = toBits('0110100101100101');
  const pair = seededKey(16, 1, 30);
  const target = encryptBits(pair.publicWeights, bits);
  // The trapdoor still works perfectly on it: the scheme is fine, the ATTACK
  // is what failed, and conflating those is the whole point of this check.
  assert.deepEqual(
    decryptWithTrapdoor(pair.privateKey, target).bits,
    bits,
    'the trapdoor must still decrypt the instance the attack cannot break',
  );
  const result = attackPublicOnly({ publicWeights: pair.publicWeights, target, delta: DELTA });
  assert.equal(result.outcome, 'no-candidate', 'seed 30 is the pinned failure');
  assert.equal(result.bits, null, 'a failure must never carry bits');
  assert.equal(result.rejections.length > 0, true, 'and must say which rows it rejected');
  assert.equal(
    result.rowsInspected,
    2 * 17,
    'a failure must have inspected every row of both Algorithm SV passes',
  );
});

check('exhaustive enumeration is verifier-side ground truth and refuses large n', () => {
  const f = PINNED[0];
  const truth = enumerateSolutions(f.publicWeights, f.target);
  assert.equal(truth.count, 1, 'the pinned instance has exactly one solution');
  assert.deepEqual(truth.solutions[0], toBits(f.bits));

  // A deliberately ambiguous instance, so the ambiguity branch is exercised on
  // something real rather than asserted to be unreachable. These are not
  // Merkle-Hellman weights and are not presented as any: 1+4 = 2+3.
  const ambiguous = enumerateSolutions([1n, 2n, 3n, 4n], 5n);
  assert.equal(ambiguous.count, 2, '1+4 and 2+3 both make 5');

  assert.match(
    expectThrow(() => enumerateSolutions(new Array(17).fill(1n), 3n)).message,
    /refused above 16 weights/,
    'enumeration must refuse rather than hang',
  );
});

/* ------------------------------------------------------------------ */

function validate(key: MerkleHellmanPrivateKey): void {
  // Routed through the public entry point so the checks exercise what the app
  // calls, not a private helper.
  derivePublicWeights(key);
}

const onlyIndex = process.argv.indexOf('--only');
const only = onlyIndex === -1 ? null : process.argv[onlyIndex + 1];
const selected = only === null ? checks : checks.filter((c) => c.name === only);

if (only !== null && selected.length !== 1) {
  console.error(`# checks ${selected.length}`);
  console.error(`No unique check named ${JSON.stringify(only)} (matched ${selected.length})`);
  process.exit(2);
}

let pass = 0;
let fail = 0;
for (const c of selected) {
  try {
    c.run();
    pass += 1;
    console.log(`ok - ${c.name}`);
  } catch (error) {
    fail += 1;
    const err = error as NodeJS.ErrnoException;
    console.log(`not ok - ${c.name}`);
    console.log(`# error code ${err.code ?? 'NONE'}`);
    console.log(String(err.stack ?? error));
  }
}
setSeed(null);
console.log(`# checks ${selected.length}`);
console.log(`# pass ${pass}`);
console.log(`# fail ${fail}`);
if (fail > 0) process.exit(1);
console.log('verify:knapsack passed');
