import { expect, test } from '@playwright/test';
import { boot, driveAllStates, NARROW, reportCollected, watchPageErrors } from './gate';

/**
 * WCAG A/AA regression gate.
 *
 * The lab is driven along everything it teaches, and every state is scanned:
 * the arrival state, where Exhibit 4 is entirely empty and its collapse control
 * ships `disabled`; both skip links focused through real Tab presses; the two
 * info tooltips opened by their buttons; Exhibit 1's four sliders taken to the
 * extreme that collapses the determinant to zero, then the unimodular re-basis;
 * the Gram-Schmidt step-through cycled all the way round and then re-run on an
 * already-orthogonal basis; all five LLL presets, both dimensions, delta at
 * 0.999 auto-run to convergence on a Fibonacci basis (which is the only thing
 * that makes the step log overflow at all) and at 0.5, plus the input-error
 * state a non-square matrix produces; Exhibit 4's prerequisite message, then a
 * seeded instance, the embedding-collapse replay run to its end state, a
 * successful attack, a FAILED attack at n=10/sigma=10, and the Kyber-512
 * figures; Exhibit 5 at n=256 and at its minimum; all five canonical labs; a
 * challenge set-up button; and the finished page with all twelve disclosures
 * open. Each of those is scanned in {dark, light} × {1280px, 380px}.
 *
 * Clipboard permission is granted because "Copy lab link" calls
 * `navigator.clipboard.writeText` inside a try/catch that falls back to a
 * DIFFERENT status message on rejection: without the grant the drive would be
 * asserting against the failure branch while claiming to scan the success one.
 *
 * See `gate.ts` for why nothing is injected into the page (this lab reads
 * `prefers-reduced-motion` in its TypeScript as well as its CSS, and a style tag
 * reaches neither), why no disclosure is force-opened, why the lab's defaults
 * are asserted rather than assumed, and why `violations` is not the whole
 * oracle on a page where every section is a gradient.
 */

for (const theme of ['dark'] as const) {
  test(`no WCAG A/AA violations in ${theme} theme`, async ({ page, context }) => {
    test.setTimeout(1_500_000);
    await context.grantPermissions(['clipboard-read', 'clipboard-write']);
    const errors = watchPageErrors(page);
    await boot(page, theme);
    await driveAllStates(page, theme);
    expect(errors, errors.join('\n')).toEqual([]);
    reportCollected();
  });

  test(`no WCAG A/AA violations in ${theme} theme at 380px`, async ({ page, context }) => {
    test.setTimeout(1_500_000);
    await context.grantPermissions(['clipboard-read', 'clipboard-write']);
    const errors = watchPageErrors(page);
    await page.setViewportSize(NARROW);
    await boot(page, theme);
    await driveAllStates(page, `${theme} @380px`);
    expect(errors, errors.join('\n')).toEqual([]);
    reportCollected();
  });
}
