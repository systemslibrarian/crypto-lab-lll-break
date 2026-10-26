/**
 * E-K-09 mutation gate for the Merkle-Hellman / Lagarias-Odlyzko exhibit.
 *
 * Fifteen load-bearing edits to `src/knapsack.ts` and four to `src/main.ts`. For
 * every one the run proves six things, and restores the source in a
 * `finally`-style path whether or not the mutation was caught:
 *
 *   1. the UNMODIFIED baseline passes the named owning check,
 *   2. the patch produces a NON-EMPTY source diff,
 *   3. the named owning check turns red,
 *   4. pure checks turn red on ERR_ASSERTION; the browser check reports exactly
 *      one named Playwright assertion failure — never a crash, missing module
 *      or empty run,
 *   5. the source is restored BYTE FOR BYTE,
 *   6. the clean control passes again after restoration.
 *
 * A mutation that cannot be applied, or that is killed by something other than
 * its own named assertion, is a failure of this harness rather than evidence.
 *
 * The mutations against `src/main.ts` cover DOM wiring and have no pure check
 * to own them. They are owned by Playwright tests instead, and because a
 * browser test runs against the BUILT bundle this
 * harness rebuilds and asserts the bundle hash actually moved before believing
 * the run -- a mutation applied to a file the build did not pick up is
 * indistinguishable from a mutation nothing caught.
 *
 * Run from the repository root:  node scripts/knapsack-mutations.mjs [name]
 */
import { createHash } from 'node:crypto';
import { spawnSync } from 'node:child_process';
import { readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

const KNAPSACK = 'src/knapsack.ts';
const MAIN = 'src/main.ts';
const CHECKS = 'scripts/knapsack-verify.ts';
const UI_SPEC = 'tests/ui.spec.ts';

/** @type {{name:string, file:string, anchor:string, replacement:string, check:string, kind:'pure'|'ui'}[]} */
const mutations = [
  {
    name: 'superincreasing validation bypassed',
    file: KNAPSACK,
    anchor: '  const bad = firstNonSuperincreasingIndex(privateWeights);\n  if (bad !== -1) {',
    replacement: '  const bad = firstNonSuperincreasingIndex(privateWeights);\n  if (false && bad !== -1) {',
    check: 'a private sequence that is not superincreasing is refused',
    kind: 'pure',
  },
  {
    name: 'the modulus condition m > sum(a_prime) is not enforced',
    file: KNAPSACK,
    anchor: '  if (modulus <= total) {',
    replacement: '  if (false && modulus <= total) {',
    check: 'a modulus that does not exceed the private-weight sum is refused',
    kind: 'pure',
  },
  {
    name: 'a non-invertible multiplier is accepted',
    file: KNAPSACK,
    anchor: '  if (inverse === null) {',
    replacement: '  if (false && inverse === null) {',
    check: 'a multiplier that is not invertible modulo m is refused',
    kind: 'pure',
  },
  {
    name: 'the greedy pass uses a strict inequality',
    file: KNAPSACK,
    anchor: '    if (residue >= key.privateWeights[i]) {',
    replacement: '    if (residue > key.privateWeights[i]) {',
    check: "the greedy pass uses >= and not >, as the paper's own example requires",
    kind: 'pure',
  },
  {
    name: 'public weights are derived from a single fixed index',
    file: KNAPSACK,
    anchor: '  return key.privateWeights.map((a) => (key.multiplier * a) % key.modulus);',
    replacement: '  return key.privateWeights.map(() => (key.multiplier * key.privateWeights[0]) % key.modulus);',
    check: 'public weights are w*a_prime mod m at every index',
    kind: 'pure',
  },
  {
    name: 'the public-sum equality check always passes',
    file: KNAPSACK,
    anchor: '    if (bits[i] === 1) sum += publicWeights[i];',
    replacement: '    if (bits[i] === 1) sum += publicWeights[i];\n    if (bits.length > 0) return true;',
    check: 'a candidate whose public sum differs from the target is rejected',
    kind: 'pure',
  },
  {
    name: 'bit-coordinate validation accepts mixed non-zero values',
    file: KNAPSACK,
    anchor: '    else if (v !== lambda) {',
    replacement: '    else if (false && v !== lambda) {',
    check: 'candidate coordinates outside the permitted bit set are rejected',
    kind: 'pure',
  },
  {
    name: 'the exact-arithmetic budget is not enforced',
    file: KNAPSACK,
    anchor: '  if (!budget.withinBudget) {',
    replacement: '  if (false && !budget.withinBudget) {',
    check: 'an instance outside the exact-arithmetic budget is refused, not attacked',
    kind: 'pure',
  },
  {
    name: 'the arithmetic screen forgets the Algorithm SV complement target',
    file: KNAPSACK,
    anchor: '  if (absolute(complementTarget) > max) max = absolute(complementTarget);',
    replacement: '  if (false && absolute(complementTarget) > max) max = absolute(complementTarget);',
    check: 'an instance outside the exact-arithmetic budget is refused, not attacked',
    kind: 'pure',
  },
  {
    name: 'the published Algorithm SV target domain is bypassed',
    file: KNAPSACK,
    anchor: '  if (target < 1n || target >= total) {',
    replacement: '  if (false && (target < 1n || target >= total)) {',
    check: 'Algorithm SV refuses weights and targets outside its published domain',
    kind: 'pure',
  },
  {
    name: 'the delta range of Lenstra-Lenstra-Lovasz (1.5) is not enforced',
    file: KNAPSACK,
    anchor: '  if (!Number.isFinite(delta) || delta <= 0.25 || delta >= 1) {',
    replacement: '  if (false) {',
    check: 'delta outside the Lenstra-Lenstra-Lovasz range is refused',
    kind: 'pure',
  },
  {
    name: 'the Algorithm SV basis silently becomes a scaled Coster-style one',
    file: KNAPSACK,
    anchor: '    row[n] = -toSafeInteger(publicWeights[i], `public weight a_${i + 1}`);',
    replacement: '    row[n] = 2 * toSafeInteger(publicWeights[i], `public weight a_${i + 1}`);',
    check: 'the Lagarias-Odlyzko basis is the printed (2.4) and contains the plaintext vector',
    kind: 'pure',
  },
  {
    name: 'the exact postcheck accepts a non-unimodular transform',
    file: KNAPSACK,
    anchor: '  if (det !== 1n && det !== -1n) {',
    replacement: '  if (false && det !== 1n && det !== -1n) {',
    check: 'the reduced basis is postchecked as an exact unimodular transform',
    kind: 'pure',
  },
  {
    name: 'the attack is handed a non-unimodular transform at its exact postcheck',
    file: KNAPSACK,
    anchor: '      reduction.transform,',
    replacement:
      '      reduction.transform.map((row, index) => index === 0 ? row.map((value) => 2 * value) : row),',
    check: 'the public-only attack recovers the message on pinned low-density instances',
    kind: 'pure',
  },
  {
    name: 'Algorithm SV step 4 repeats pass 1 instead of using the complement target',
    file: KNAPSACK,
    anchor: '    { pass: 2, passTarget: total - target, complement: true },',
    replacement: '    { pass: 2, passTarget: target, complement: true },',
    check: 'Algorithm SV step (4) recovers a fixture that needs the complement pass',
    kind: 'pure',
  },
  {
    name: 'changing the message no longer retires the displayed result',
    file: MAIN,
    anchor:
      "mhBits.addEventListener('input', () => {\n  mhRetire('Message changed. Previous key, ciphertext and verdict retired. Press Generate.');\n});",
    replacement: "mhBits.addEventListener('input', () => {\n  mhStatus.textContent = 'Message changed.';\n});",
    check: 'changing the message retires the key, the lattice and the verdict',
    kind: 'ui',
  },
  {
    name: 'seeded knapsack generation no longer rewinds its random stream',
    file: MAIN,
    anchor: '  if (activeSeed !== null) setSeed(activeSeed);',
    replacement: '  if (false && activeSeed !== null) setSeed(activeSeed);',
    check: 'seeded knapsack generation rewinds to the same instance',
    kind: 'ui',
  },
  {
    name: 'shared lab links omit the knapsack dimension',
    file: MAIN,
    anchor: "  p.set('mhn', mhN.value);",
    replacement: "  // p.set('mhn', mhN.value);",
    check: 'a copied lab link preserves every knapsack control',
    kind: 'ui',
  },
  {
    name: 'an Algorithm SV endpoint is labelled eligible for LLL',
    file: MAIN,
    anchor: '  if (budget.withinBudget && inAlgorithmSvDomain) {',
    replacement: '  if (budget.withinBudget) {',
    check: 'the all-one endpoint decrypts but is refused as outside Algorithm SV',
    kind: 'ui',
  },
];

function run(command, args, timeout = 900_000, env = process.env) {
  const result = spawnSync(command, args, {
    encoding: 'utf8',
    maxBuffer: 64 * 1024 * 1024,
    timeout,
    env,
  });
  if (result.error || result.signal || !Number.isInteger(result.status)) {
    throw new Error(
      `${command}: process did not exit normally (${result.error?.message ?? result.signal ?? 'no exit status'})`,
    );
  }
  const stdout = result.stdout ?? '';
  const stderr = result.stderr ?? '';
  return { code: result.status, stdout, stderr, output: `${stdout}${stderr}` };
}

const INFRASTRUCTURE_MARKERS = [
  'Cannot find module',
  'ERR_MODULE_NOT_FOUND',
  'SyntaxError',
  'ReferenceError',
];

function rejectInfrastructureFailure(name, output) {
  for (const marker of INFRASTRUCTURE_MARKERS) {
    if (output.includes(marker)) {
      throw new Error(
        `${name}: "${marker}" is an infrastructure failure, not a mutation kill\n${output.slice(0, 2000)}`,
      );
    }
  }
}

/** Run exactly the one named pure check and insist the result is unambiguous. */
function namedPureCheck(mutation, shouldPass) {
  const result = run('npx', ['tsx', CHECKS, '--only', mutation.check]);
  const total = Number(/^# checks (\d+)$/m.exec(result.output)?.[1] ?? -1);
  const pass = Number(/^# pass (\d+)$/m.exec(result.output)?.[1] ?? -1);
  const fail = Number(/^# fail (\d+)$/m.exec(result.output)?.[1] ?? -1);
  rejectInfrastructureFailure(mutation.name, result.output);
  if (total !== 1) {
    throw new Error(
      `${mutation.name}: expected exactly one check named ${JSON.stringify(mutation.check)}, the runner reported ${total}`,
    );
  }
  if (shouldPass) {
    if (result.code !== 0 || pass !== 1 || fail !== 0) {
      throw new Error(
        `${mutation.name}: expected a green baseline, got exit ${result.code} pass ${pass} fail ${fail}\n${result.output}`,
      );
    }
    return;
  }
  if (result.code === 0 || fail !== 1 || pass !== 0) {
    throw new Error(
      `${mutation.name}: the mutation SURVIVED (exit ${result.code}, pass ${pass}, fail ${fail})\n${result.output.slice(0, 2000)}`,
    );
  }
  if (!/^# error code ERR_ASSERTION$/m.test(result.output)) {
    throw new Error(
      `${mutation.name}: the owning check went red without an assertion failure, which cannot count as a kill\n${result.output.slice(0, 2000)}`,
    );
  }
}

/** Hash of the built bundle, so a mutation the build never picked up is visible. */
function bundleDigest() {
  const dir = 'dist/assets';
  let names;
  try {
    names = readdirSync(dir).sort();
  } catch {
    return 'NO-DIST';
  }
  const hash = createHash('sha256');
  for (const name of names) {
    hash.update(name);
    hash.update(readFileSync(join(dir, name)));
  }
  return hash.digest('hex');
}

/**
 * Run exactly the one named browser test against a freshly built bundle.
 *
 * `vite build` is run HERE rather than left to Playwright's webServer, because
 * locally `reuseExistingServer` means Playwright may attach to a preview server
 * that is already up and never rebuild -- which would test the unmutated bundle
 * and read as "the mutation was not caught".
 */
function namedBrowserCheck(mutation, shouldPass) {
  const build = run('npx', ['tsc', '--noEmit']);
  if (build.code !== 0) {
    throw new Error(`${mutation.name}: type check failed\n${build.output.slice(0, 2000)}`);
  }
  const bundle = run('npx', ['vite', 'build']);
  if (bundle.code !== 0) {
    throw new Error(`${mutation.name}: vite build failed\n${bundle.output.slice(0, 2000)}`);
  }
  const digest = bundleDigest();

  const result = run('npx', [
    'playwright',
    'test',
    '--config',
    'playwright.config.ts',
    '-g',
    mutation.check,
    '--reporter',
    'line',
  ]);
  rejectInfrastructureFailure(mutation.name, result.output);
  if (/0 tests? (?:ran|found)|No tests found/i.test(result.output)) {
    throw new Error(`${mutation.name}: no browser test matched ${JSON.stringify(mutation.check)}`);
  }
  if (shouldPass) {
    if (result.code !== 0 || !/\b1 passed\b/.test(result.output) || /\b1 failed\b/.test(result.output)) {
      throw new Error(`${mutation.name}: expected a green browser baseline\n${result.output.slice(0, 3000)}`);
    }
  } else {
    if (result.code === 0) {
      throw new Error(`${mutation.name}: the mutation SURVIVED its browser test\n${result.output.slice(0, 3000)}`);
    }
    if (!/\b1 failed\b/.test(result.output) || !result.output.includes(mutation.check)) {
      throw new Error(
        `${mutation.name}: the named browser test did not fail exactly once\n${result.output.slice(0, 3000)}`,
      );
    }
    if (!/Error: expect\(/.test(result.output)) {
      throw new Error(
        `${mutation.name}: browser run went red without a Playwright assertion failure\n${result.output.slice(0, 3000)}`,
      );
    }
  }
  return digest;
}

const digestOf = (text) => createHash('sha256').update(text).digest('hex');

const requested = process.argv[2];
const selected = requested ? mutations.filter((m) => m.name === requested) : mutations;
if (selected.length === 0) {
  console.error(`Unknown mutation: ${requested}`);
  process.exit(2);
}

for (const file of [KNAPSACK, MAIN, CHECKS, UI_SPEC]) {
  if (run('git', ['status', '--porcelain', '--', file]).output.trim()) {
    throw new Error(
      `${file} has an existing working-tree change. Preserve or deliberately restore it before ` +
        'replaying mutations; this harness will not overwrite unrelated work.',
    );
  }
}

let killed = 0;
for (const mutation of selected) {
  const original = readFileSync(mutation.file, 'utf8');
  const originalDigest = digestOf(original);
  if (original.split(mutation.anchor).length !== 2) {
    throw new Error(`${mutation.name}: its anchor must occur exactly once in ${mutation.file}`);
  }

  if (mutation.kind === 'pure') namedPureCheck(mutation, true);
  else namedBrowserCheck(mutation, true);
  console.log(`[${mutation.name}]`);
  console.log(`  1/6 baseline GREEN: ${mutation.check}`);

  let bundleBefore = null;
  let bundleAfter = null;
  try {
    const changed = original.replace(mutation.anchor, mutation.replacement);
    if (changed === original) throw new Error(`${mutation.name}: the patch changed nothing`);
    if (mutation.kind === 'ui') bundleBefore = bundleDigest();
    writeFileSync(mutation.file, changed);

    const diff = run('git', ['diff', '--no-ext-diff', '--', mutation.file], 120_000);
    const changedLines = diff.output.split('\n').filter((l) => /^[+-][^+-]/.test(l)).length;
    if (diff.code !== 0 || changedLines === 0) {
      throw new Error(`${mutation.name}: the source diff was empty`);
    }
    console.log(`  2/6 source diff NON-EMPTY (${changedLines} changed lines in ${mutation.file})`);

    if (mutation.kind === 'pure') {
      namedPureCheck(mutation, false);
    } else {
      bundleAfter = namedBrowserCheck(mutation, false);
      if (bundleBefore !== null && bundleAfter === bundleBefore) {
        throw new Error(
          `${mutation.name}: the built bundle hash did not move, so the browser test ran against ` +
            'unmutated code and proves nothing',
        );
      }
      console.log('  2b/6 built bundle hash MOVED, so the browser really ran the mutated code');
    }
    console.log('  3/6 owning check RED');
    console.log(
      mutation.kind === 'pure'
        ? '  4/6 red on ERR_ASSERTION, with exactly one check matched'
        : '  4/6 exactly one named browser test failed on a Playwright assertion',
    );
  } finally {
    writeFileSync(mutation.file, original);
  }

  if (digestOf(readFileSync(mutation.file, 'utf8')) !== originalDigest) {
    throw new Error(`${mutation.name}: the restored source differs from the original`);
  }
  console.log('  5/6 source restored BYTE FOR BYTE');

  if (mutation.kind === 'pure') namedPureCheck(mutation, true);
  else namedBrowserCheck(mutation, true);
  console.log('  6/6 clean control GREEN again after restoration');
  killed += 1;
}

console.log('');
console.log(`BASELINE PASSED for all ${selected.length} mutations.`);
console.log(`ALL MUTATIONS CAUGHT: ${killed} of ${selected.length} killed by their own named checks.`);
if (killed !== selected.length) process.exit(1);
