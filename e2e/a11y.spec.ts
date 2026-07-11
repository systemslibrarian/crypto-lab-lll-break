import AxeBuilder from '@axe-core/playwright';
import { expect, test, type Page } from '@playwright/test';

/**
 * WCAG regression gate. Deploys are already gated on the algorithm + guardrail
 * verifiers; this gates them on accessibility the same way. Scans the full page
 * with every collapsible expanded and the interactive demos driven so their
 * generated output is in the DOM, in both themes. Zero tolerance: any WCAG
 * A/AA violation fails the build.
 */

const TAGS = ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'];

/** Neutralise CSS transitions/animations so scans see settled, stable state. */
async function killMotion(page: Page): Promise<void> {
  await page.addStyleTag({
    content: `*,*::before,*::after{
      transition:none!important;
      animation:none!important;
      scroll-behavior:auto!important;
    }`,
  });
}

/** Expand every <details> so collapsed content is scanned too. */
async function openAllDetails(page: Page): Promise<void> {
  await page.evaluate(() => {
    for (const details of document.querySelectorAll('details')) {
      details.open = true;
    }
  });
}

/**
 * Drive the interactive demos so axe scans their generated output, not just the
 * static shell. Each step is best-effort (guarded) so a UI tweak can't turn the
 * gate red for the wrong reason.
 */
async function driveDemos(page: Page): Promise<void> {
  const click = async (name: RegExp | string) => {
    const btn = page.getByRole('button', { name });
    if (await btn.count()) await btn.first().click();
  };

  // Exhibit 2/3: Gram-Schmidt + LLL step-through.
  await click(/^Step$/);
  await click(/^Step$/);

  // Exhibit 4: seed, generate a toy LWE instance, run the reduction attack.
  const seed = page.locator('#lab-seed');
  if (await seed.count()) {
    await seed.fill('1234');
    await click(/Apply seed/i);
  }
  await click(/Generate LWE Instance/i);
  await expect(page.locator('#lwe-output')).toContainText(/Secret s =/i);
  await click(/Run LLL\/BKZ Attack/i);
  await expect(page.locator('#lwe-output')).toContainText(/Lattice attack result:/i);
}

async function scan(page: Page): Promise<void> {
  const results = await new AxeBuilder({ page }).withTags(TAGS).analyze();
  const summary = results.violations.map((v) => ({
    id: v.id,
    impact: v.impact,
    help: v.help,
    nodes: v.nodes.map((n) => n.target.join(' ')).slice(0, 5),
  }));
  expect(summary).toEqual([]);
}

test('no WCAG A/AA violations in dark theme', async ({ page }) => {
  await page.goto('.');
  await killMotion(page);
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'dark');
  await driveDemos(page);
  await openAllDetails(page);
  await scan(page);
});

test('no WCAG A/AA violations in light theme', async ({ page }) => {
  await page.goto('.');
  await killMotion(page);
  await page.locator('#cl-theme-toggle').click();
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
  await driveDemos(page);
  await openAllDetails(page);
  await scan(page);
});
