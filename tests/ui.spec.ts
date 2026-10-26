import { test, expect } from '@playwright/test';
import AxeBuilder from '@axe-core/playwright';

test.beforeEach(async ({ page }) => {
  await page.goto('.');
  // Settle the `rise` fade-in before anything scans: axe otherwise samples the
  // exhibits mid-animation at partial opacity and reports phantom
  // color-contrast failures. (Same killMotion tactic as the e2e a11y gate;
  // the config-level reducedMotion option is not honored by this runner.)
  await page.addStyleTag({
    content: `*,*::before,*::after{animation:none!important;transition:none!important}`,
  });
  // The app renders into #app via main.ts; wait for the first exhibit.
  await expect(page.getByRole('heading', { name: /What Is a Lattice/i })).toBeVisible();
});

test('renders all exhibits and the labs/challenges sections', async ({ page }) => {
  for (const name of [
    /Reproducible Labs/i,
    /What Is a Lattice/i,
    /Gram-Schmidt/i,
    /LLL Step-by-Step/i,
    /Break a Toy LWE Instance/i,
    /Parameter Explorer/i,
    /Knapsack Trapdoors/i,
    /Challenges/i,
  ]) {
    await expect(page.getByRole('heading', { name })).toBeVisible();
  }
});

test('exactly one banner landmark and a working skip link', async ({ page }) => {
  // The shared site header is the only banner; the page title bar must not be one.
  await expect(page.getByRole('banner')).toHaveCount(1);
  const skip = page.getByRole('link', { name: /skip to main content/i });
  await expect(skip).toHaveCount(1);
  await expect(page.locator('#main-content')).toHaveCount(1);
});

test('LLL step shows the unimodular-transform invariant', async ({ page }) => {
  await page.getByRole('button', { name: 'Step', exact: true }).click();
  const metrics = page.locator('#lll-metrics');
  await expect(metrics).toContainText(/det\(U\)=/);
  await expect(metrics).toContainText(/unchanged: same lattice/);
});

test('toy LWE attack recovers by genuine reduction, not the baseline', async ({ page }) => {
  await page.locator('#lab-seed').fill('1234');
  await page.getByRole('button', { name: 'Apply seed' }).click();
  await page.getByRole('button', { name: 'Generate LWE Instance' }).click();
  await expect(page.locator('#lwe-output')).toContainText(/Secret s =/);
  await page.getByRole('button', { name: 'Run LLL/BKZ Attack' }).click();
  const out = page.locator('#lwe-output');
  await expect(out).toContainText(/Lattice attack result: SUCCESS/);
  await expect(out).toContainText(/EXACT MATCH/);
});

test('moving sliders after Generate does not desync n/q/sigma (bug-1 guard)', async ({
  page,
}) => {
  await page.locator('#lwe-n').fill('4');
  await page.getByRole('button', { name: 'Generate LWE Instance' }).click();
  await expect(page.locator('#lwe-output')).toContainText(/Secret s =/);
  // Change n on the slider AFTER generating; the attack must use the instance's n.
  await page.locator('#lwe-n').fill('8');
  await page.getByRole('button', { name: 'Run LLL/BKZ Attack' }).click();
  const out = page.locator('#lwe-output');
  await expect(out).toContainText(/Lattice attack result:/);
  // The recovered secret must have length 4 (the generated n), not 8.
  const text = (await out.textContent()) ?? '';
  const match = text.match(/Recovered secret: \[([^\]]*)\]/);
  if (match) {
    const count = match[1].split(',').filter((s) => s.trim().length > 0).length;
    expect(count).toBe(4);
  }
});

test('Kyber-512 button never reports a lattice recovery', async ({ page }) => {
  await page.getByRole('button', { name: /Try Kyber-512 parameters/i }).click();
  const out = page.locator('#lwe-output');
  await expect(out).toContainText(/Kyber-like test/i);
  await expect(out).toContainText(/does not recover/i);
  await expect(out).not.toContainText(/Lattice attack result: SUCCESS/);
});

test('seeded generation is reproducible from the UI', async ({ page }) => {
  const readSecret = async () => {
    const t = (await page.locator('#lwe-output').textContent()) ?? '';
    return (t.match(/Secret s = \[([^\]]*)\]/) ?? [])[1] ?? '';
  };
  await page.locator('#lab-seed').fill('99');
  await page.getByRole('button', { name: 'Apply seed' }).click();
  await page.getByRole('button', { name: 'Generate LWE Instance' }).click();
  const first = await readSecret();
  await page.getByRole('button', { name: 'Generate LWE Instance' }).click();
  const second = await readSecret();
  expect(first.length).toBeGreaterThan(0);
  expect(second).toBe(first);
});

test('canonical lab preset loads a reproducible setup', async ({ page }) => {
  await page.locator('#lab-preset').selectOption('toy-lwe-success');
  await expect(page.locator('#lwe-output')).toContainText(/Secret s =/);
  await expect(page.locator('#lab-status')).toContainText(/seeded/i);
});

test('no critical or serious accessibility violations', async ({ page }) => {
  const results = await new AxeBuilder({ page })
    .withTags(['wcag2a', 'wcag2aa'])
    .analyze();
  const blocking = results.violations.filter(
    (v) => v.impact === 'critical' || v.impact === 'serious',
  );
  expect(
    blocking,
    blocking.map((v) => `${v.id}: ${v.help}`).join('\n'),
  ).toEqual([]);
});

test('SUCCESS is never printed over a vector that is not the secret', async ({ page }) => {
  // The defect this pins: the verdict was `ok && recoverMethod === 'short-vector'`
  // and never consulted the byte comparison computed a few lines later. `ok` is a
  // sigma-scaled residual check, so at higher noise an arbitrary short vector
  // passes it and the page printed "SUCCESS - secret read off a reduced short
  // vector" above a recovered vector unrelated to the secret. Challenge 5 walks
  // learners straight into this regime and promises a FAILED verdict.
  const out = page.locator('#lwe-output');

  for (const sigma of ['5', '7', '9', '10']) {
    await page.locator('#lab-seed').fill('1234');
    await page.getByRole('button', { name: 'Apply seed' }).click();
    await page.locator('#lwe-sigma').fill(sigma);
    await page.getByRole('button', { name: 'Generate LWE Instance' }).click();
    await expect(out).toContainText(/Secret s =/);
    await page.getByRole('button', { name: 'Run LLL\/BKZ Attack' }).click();
    await expect(out).toContainText(/Lattice attack result:/);

    const text = (await out.textContent()) ?? '';
    if (/Lattice attack result: SUCCESS/.test(text)) {
      // A SUCCESS headline must be backed by the exact-match line.
      expect(text, `sigma ${sigma}: SUCCESS without EXACT MATCH`).toContain('EXACT MATCH');
      expect(text, `sigma ${sigma}: SUCCESS alongside a non-secret vector`).not.toContain(
        'NOT the secret',
      );
    }
    if (/NOT the secret/.test(text)) {
      // ...and a non-secret vector must never carry a SUCCESS headline.
      expect(text, `sigma ${sigma}: non-secret vector reported as SUCCESS`).not.toContain(
        'Lattice attack result: SUCCESS',
      );
    }
  }
});


/* ---------------------------------------------------------------------- */
/* Exhibit 6 - Merkle-Hellman trapdoor knapsacks (E-K-09)                   */
/* ---------------------------------------------------------------------- */

/** Seed, then generate the pinned n=8 / k=14 instance the assertions below use. */
async function knapsackSetup(page: import('@playwright/test').Page): Promise<void> {
  await page.locator('#lab-seed').fill('1234');
  await page.getByRole('button', { name: 'Apply seed' }).click();
  await page.locator('#mh-n').selectOption('8');
  await page.locator('#mh-scale').fill('14');
  await page.locator('#mh-bits').fill('01101001');
  await page.getByRole('button', { name: 'Generate Knapsack Key' }).click();
  await expect(page.locator('#mh-keys')).toContainText('PUBLIC (everything an attacker is given)');
}

/** Pull the public weights and ciphertext straight out of what the page rendered. */
function readPublicInstance(text: string): { weights: bigint[]; target: bigint } {
  const weights = (
    text.match(/a_i = w \* a'_i mod m, Merkle-Hellman \(3\)\n\s*\[([^\]]*)\]/) ?? []
  )[1];
  const target = (text.match(/ciphertext S = sum\(x_i \* a_i\) = (\d+)/) ?? [])[1];
  return {
    weights: (weights ?? '').split(',').map((w) => BigInt(w.trim())),
    target: BigInt(target ?? '0'),
  };
}

test('the knapsack exhibit ships empty, with both downstream controls locked', async ({ page }) => {
  await expect(page.locator('#mh-keys')).toBeEmpty();
  await expect(page.locator('#mh-basis')).toBeEmpty();
  await expect(page.locator('#mh-output')).toBeEmpty();
  await expect(page.locator('#mh-decrypt')).toBeDisabled();
  await expect(page.locator('#mh-attack')).toBeDisabled();
});

test('the generated key satisfies the two Merkle-Hellman conditions on screen', async ({ page }) => {
  await knapsackSetup(page);
  const text = (await page.locator('#mh-keys').textContent()) ?? '';

  // (8) is displayed as a computed boolean, not asserted in prose.
  expect(text).toContain("(m > sum(a'): true, Merkle-Hellman (8))");

  // Recompute the density from the weights the page itself printed, rather
  // than trusting the number beside them.
  const { weights, target } = readPublicInstance(text);
  expect(weights).toHaveLength(8);
  const max = weights.reduce((a, b) => (a > b ? a : b));
  const recomputed = 8 / Math.log2(Number(max));
  const shown = Number((text.match(/= (\d\.\d{4})\n/) ?? [])[1] ?? NaN);
  expect(Math.abs(recomputed - shown)).toBeLessThan(5e-5);

  // The ciphertext must be the subset sum of the bits that were encrypted.
  const bits = '01101001'.split('').map(Number);
  let sum = 0n;
  for (let i = 0; i < 8; i += 1) if (bits[i] === 1) sum += weights[i];
  expect(sum).toBe(target);

  // And the embedding really contains the plaintext vector.
  await expect(page.locator('#mh-basis')).toContainText(
    'last coordinate = M - sum(e_i a_i) = 0  (must be 0)',
  );
});

test('trapdoor decryption reproduces the message that was encrypted', async ({ page }) => {
  await knapsackSetup(page);
  await page.getByRole('button', { name: 'Decrypt with Trapdoor' }).click();
  const out = page.locator('#mh-output');
  await expect(out).toContainText('TRAPDOOR ROUTE');
  await expect(out).toContainText('recovered bits = 01101001');
  await expect(out).toContainText('remainder after the pass = 0 (must be 0: true)');
  await expect(out).toContainText('EXACT MATCH with the message that was sent');
});

test('the knapsack lattice attack recovers message bits from public data alone', async ({ page }) => {
  await knapsackSetup(page);
  const keys = (await page.locator('#mh-keys').textContent()) ?? '';
  const { weights, target } = readPublicInstance(keys);

  await page.getByRole('button', { name: 'Run LLL Attack (public data only)' }).click();
  const out = page.locator('#mh-output');
  await expect(out).toContainText('Lattice attack result: SUCCESS');

  const text = (await out.textContent()) ?? '';
  // What the attacker was withheld is stated, and the claim is checked here by
  // recomputing the subset sum of the bits the page says it recovered.
  expect(text).toContain('withheld: the private sequence, the modulus, the multiplier');
  const recovered = (text.match(/bits = ([01]+)/) ?? [])[1] ?? '';
  expect(recovered).toHaveLength(8);
  let sum = 0n;
  for (let i = 0; i < 8; i += 1) if (recovered[i] === '1') sum += weights[i];
  expect(sum, 'the recovered bits must reproduce the ciphertext exactly').toBe(target);
  expect(text).toContain('EXACT MATCH with the sent message');

  // The two things this exhibit must never be read as claiming.
  expect(text).toContain('The trapdoor is untouched.');
  expect(text).toContain("That is Shamir's attack, which is a different attack");
});

test('a SUCCESS headline is never printed without the exact-sum verification', async ({ page }) => {
  await knapsackSetup(page);
  await page.getByRole('button', { name: 'Run LLL Attack (public data only)' }).click();
  const text = (await page.locator('#mh-output').textContent()) ?? '';
  if (/Lattice attack result: SUCCESS/.test(text)) {
    expect(text).toContain('Verified in exact integer arithmetic: sum(bits_i * a_i) == S  ->  true');
  }
  if (/Lattice attack result: FAILED/.test(text)) {
    expect(text).toContain('This is NOT a security result');
    expect(text).not.toContain('Lattice attack result: SUCCESS');
  }
});

test('changing the message retires the key, the lattice and the verdict', async ({ page }) => {
  await knapsackSetup(page);
  await page.getByRole('button', { name: 'Run LLL Attack (public data only)' }).click();
  await expect(page.locator('#mh-output')).toContainText('Lattice attack result:');

  // One keystroke in the message box invalidates the key, the ciphertext, the
  // lattice and the verdict together. A stale verdict beside changed inputs is
  // a number computed from something the reader is no longer looking at.
  await page.locator('#mh-bits').fill('11101001');

  await expect(page.locator('#mh-output')).toBeEmpty();
  await expect(page.locator('#mh-keys')).toBeEmpty();
  await expect(page.locator('#mh-basis')).toBeEmpty();
  await expect(page.locator('#mh-decrypt')).toBeDisabled();
  await expect(page.locator('#mh-attack')).toBeDisabled();
  await expect(page.locator('#mh-status')).toContainText('retired');
});

test('changing n or the scale also retires the displayed result', async ({ page }) => {
  for (const change of [
    async () => page.locator('#mh-n').selectOption('10'),
    async () => page.locator('#mh-scale').fill('9'),
  ]) {
    await knapsackSetup(page);
    await page.getByRole('button', { name: 'Decrypt with Trapdoor' }).click();
    await expect(page.locator('#mh-output')).toContainText('TRAPDOOR ROUTE');
    await change();
    await expect(page.locator('#mh-output')).toBeEmpty();
    await expect(page.locator('#mh-attack')).toBeDisabled();
  }
});

test('a malformed message is refused by name and generates nothing', async ({ page }) => {
  await page.locator('#mh-bits').fill('0110100');
  await page.getByRole('button', { name: 'Generate Knapsack Key' }).click();
  await expect(page.locator('#mh-status')).toContainText('Input error');
  await expect(page.locator('#mh-keys')).toBeEmpty();
  await expect(page.locator('#mh-attack')).toBeDisabled();

  await page.locator('#mh-bits').fill('0110100x');
  await page.getByRole('button', { name: 'Generate Knapsack Key' }).click();
  await expect(page.locator('#mh-status')).toContainText('Input error');
  await expect(page.locator('#mh-keys')).toBeEmpty();
});

test('the knapsack exhibit is operable from the keyboard alone', async ({ page }) => {
  await page.locator('#lab-seed').fill('1234');
  await page.getByRole('button', { name: 'Apply seed' }).click();
  const generate = page.locator('#mh-generate');
  await generate.focus();
  await expect(generate).toBeFocused();
  await page.keyboard.press('Enter');
  await expect(page.locator('#mh-keys')).toContainText('PUBLIC');
  const attack = page.locator('#mh-attack');
  await attack.focus();
  await expect(attack).toBeFocused();
  await page.keyboard.press('Enter');
  await expect(page.locator('#mh-output')).toContainText('PUBLIC LATTICE ROUTE');
});

test('the knapsack exhibit states its limits without being opened', async ({ page }) => {
  // These are the claims a reader must not have to go looking for, so they are
  // in always-visible prose rather than behind a disclosure.
  const notes = page.locator('#exhibit-6 .heuristic-note');
  await expect(notes.first()).toContainText('Historical and broken');
  await expect(notes.first()).toContainText('not Shamir');
  await expect(notes.first()).toContainText('not a security result');
  await expect(notes.first()).toContainText('ML-KEM/Kyber');
  await expect(page.locator('#mh-limits')).toContainText('conditional on a shortest-vector oracle');
  await expect(page.locator('#mh-limits')).toContainText('does not display a density threshold');
});

test('the knapsack panel reflows at 380px without horizontal overflow', async ({ page }) => {
  await page.setViewportSize({ width: 380, height: 900 });
  await knapsackSetup(page);
  await page.getByRole('button', { name: 'Run LLL Attack (public data only)' }).click();
  await expect(page.locator('#mh-output')).toContainText('PUBLIC LATTICE ROUTE');
  const overflow = await page.evaluate(
    () => document.documentElement.scrollWidth - window.innerWidth,
  );
  expect(overflow, 'the page must not scroll sideways at 380px').toBeLessThanOrEqual(1);
});
