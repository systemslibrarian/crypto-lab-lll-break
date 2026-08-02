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

/**
 * WCAG 1.4.11 (non-text contrast) regression for text-entry control boundaries.
 * Axe does not flag low-contrast control borders, so we measure them directly:
 * every visible input/select/textarea's rendered border color must reach 3:1
 * against both the control's own fill and the first opaque ancestor surface
 * behind it. Translucent colors are composited against those surfaces first.
 */
async function minimumControlBoundaryRatio(page: Page): Promise<number> {
  return page.locator('input:visible, select:visible, textarea:visible').evaluateAll((elements) => {
    const parse = (value: string): { c: number[]; a: number } => {
      const n = (value.match(/[\d.]+/g) ?? ['0', '0', '0']).map(Number);
      return { c: n.slice(0, 3), a: n[3] ?? 1 };
    };
    const luminance = (parts: number[]): number => {
      const c = parts.map((part) => {
        const v = part / 255;
        return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
      });
      return 0.2126 * (c[0] ?? 0) + 0.7152 * (c[1] ?? 0) + 0.0722 * (c[2] ?? 0);
    };
    const ratio = (a: number[], b: number[]): number => {
      const [la, lb] = [luminance(a), luminance(b)];
      return (Math.max(la, lb) + 0.05) / (Math.min(la, lb) + 0.05);
    };
    const composite = (fg: number[], alpha: number, bg: number[]): number[] =>
      fg.map((v, i) => v * alpha + (bg[i] ?? 0) * (1 - alpha));
    const surfaceBehind = (el: Element): number[] => {
      for (let node = el.parentElement; node; node = node.parentElement) {
        const bg = parse(getComputedStyle(node).backgroundColor);
        if (bg.a >= 1) return bg.c;
      }
      return [255, 255, 255];
    };
    return Math.min(
      ...elements.map((el) => {
        const style = getComputedStyle(el);
        const exterior = surfaceBehind(el);
        const bg = parse(style.backgroundColor);
        const fill = bg.a >= 1 ? bg.c : composite(bg.c, bg.a, exterior);
        const b = parse(style.borderTopColor);
        const border = b.a >= 1 ? b.c : composite(b.c, b.a, fill);
        return Math.min(ratio(border, fill), ratio(border, exterior));
      }),
    );
  });
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
  expect(await minimumControlBoundaryRatio(page)).toBeGreaterThanOrEqual(3);
  await scan(page);
});

test('no WCAG A/AA violations in light theme', async ({ page }) => {
  await page.goto('.');
  await killMotion(page);
  await page.locator('#cl-theme-toggle').click();
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
  await driveDemos(page);
  await openAllDetails(page);
  expect(await minimumControlBoundaryRatio(page)).toBeGreaterThanOrEqual(3);
  await scan(page);
});
