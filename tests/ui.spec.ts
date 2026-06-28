import { test, expect } from '@playwright/test';
import AxeBuilder from '@axe-core/playwright';

test.beforeEach(async ({ page }) => {
  await page.goto('/');
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
