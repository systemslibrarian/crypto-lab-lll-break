import AxeBuilder from '@axe-core/playwright';
import { expect, type Page } from '@playwright/test';
import { auditContrast, formatContrastFailures } from './contrast';
import { auditNonText, formatNonTextFailures } from './nontext';

export const TAGS = ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'];

/** A phone-width viewport, for the WCAG 1.4.10 reflow half of the gate. */
export const NARROW = { width: 380, height: 800 };

/**
 * Shared machinery for the WCAG gate.
 *
 * Five rules govern everything here, each one a correction of the gate this
 * replaces:
 *
 *  1. NOTHING IS INJECTED INTO THE PAGE BEFORE A SCAN. `killMotion()` pushed
 *     `transition:none!important; animation:none!important` through
 *     `addStyleTag`. That BYPASSED this stylesheet's own
 *     `@media (prefers-reduced-motion: reduce)` block instead of exercising it,
 *     so the block was never once measured — and the block is the only thing
 *     that stops `.exhibit`'s `rise` keyframes from running. `rise` animates
 *     `opacity: 0 -> 1`, which is exactly the shape that strands an element
 *     invisible when a reduced-motion rule cancels an animation without
 *     restoring its end state. It does not strand it here: `animation: none`
 *     reverts `.exhibit` to its declared `opacity`, which is the initial `1`.
 *     `expectNotBlank` measures that in every state rather than trusting the
 *     reading, in a page where all six sections are `.exhibit`.
 *
 *     The preference also reaches this lab's TypeScript, not only its CSS.
 *     `main.ts` reads `matchMedia('(prefers-reduced-motion: reduce)')` once at
 *     module scope and branches on it in two places: `animateBasisTransition`
 *     draws the settled basis in a single paint instead of tweening it over
 *     320ms of `requestAnimationFrame`, and the embedding-collapse replay drops
 *     every one of its `wait()`s from 520ms to 0. A style tag cannot reproduce
 *     either. Emulating the preference properly is therefore also what makes
 *     the collapse animation's END STATE reachable synchronously, which is the
 *     state the drive actually needs to scan.
 *
 *  2. IT FORCE-REVEALED EVERY DISCLOSURE. `openAllDetails()` set `details.open
 *     = true` from script on all twelve — the five "Can you answer this now?"
 *     panels, the six challenge answers and the `.model-note`. This gate never
 *     touches `.open`; every disclosure is opened by clicking its own
 *     `<summary>`, which is the route a reader has, and the SHUT state is
 *     scanned first because that is the state the page ships in.
 *
 *  3. IT SCANNED ONCE, AT ONE VIEWPORT, AFTER A BEST-EFFORT DRIVE. Every step
 *     of `driveDemos()` was guarded with `if (await btn.count())`, so a renamed
 *     button silently skipped its exhibit and the gate stayed green having
 *     scanned first paint. Exhibits 1 and 5 were never touched at all, no
 *     slider was ever moved, the LLL presets and the 3D branch were never
 *     loaded, the embedding-collapse replay was never run, and the whole 380px
 *     column had never been scanned once. This drive asserts a state change
 *     after every step and scans after every step, in {dark, light} × {1280,
 *     380}.
 *
 *  4. `violations` IS NOT THE WHOLE ORACLE. See `scan`. On this page in
 *     particular a violations-only assertion cannot see the contrast of ANY
 *     prose, because every `.exhibit` is a `linear-gradient` and axe declines
 *     to resolve a gradient — it files the whole page under `incomplete`.
 *
 *  5. IT HAD NO REFLOW ORACLE AND ITS 1.4.11 CHECK WAS AIMED WHERE THE RULE WAS
 *     ALREADY KEPT. `minimumControlBoundaryRatio()` queried
 *     `input, select, textarea` — precisely the three selectors this
 *     stylesheet's `--line-control` token was written for and correctly applied
 *     to (measured 3.63:1 dark, 3.61:1 light). Every BUTTON on the page drew
 *     its border from `--line`, a surface divider, and was never measured
 *     against anything. See `e2e/nontext.ts`, which judges every control.
 */

/**
 * Wait for every running animation and transition to drain.
 *
 * Transitions drain in waves, not in one batch, so a poll for "nothing running
 * right now" can exit through a gap between waves. Require quiescence to hold
 * for several consecutive frames instead.
 */
export async function settle(page: Page): Promise<void> {
  await page.waitForFunction(
    () => {
      const w = window as unknown as { __quietFrames?: number };
      const running = document.getAnimations().filter((a) => a.playState === 'running');
      w.__quietFrames = running.length === 0 ? (w.__quietFrames ?? 0) + 1 : 0;
      return w.__quietFrames >= 6;
    },
    undefined,
    { timeout: 20_000, polling: 'raf' }
  );
}

/**
 * Assert that reduced motion left the page visible, not merely un-animated.
 *
 * The failure mode this guards against is an element whose only route to its
 * visible state is an animation, in a stylesheet whose reduced-motion block
 * cancels that animation without restoring its end state — the element then
 * renders at `opacity: 0` for every reader with the preference set.
 *
 * This page is one `@keyframes` away from being in exactly that shape.
 * `.exhibit` — the box around all six sections, i.e. every word the lab
 * teaches — runs `rise 0.45s`, whose `from` is `opacity: 0`, and the
 * reduced-motion block cancels it with a blanket `animation: none !important`.
 * It survives only because `.exhibit` never declares an `opacity`, so
 * cancelling the animation reverts it to the initial `1`. That is a property of
 * a declaration nobody would think to protect, so it is measured in every state
 * rather than reasoned about once.
 *
 * `aria-hidden` subtrees are excluded, since text removed from the
 * accessibility tree is not what this check is for. On this page that exclusion
 * now costs nothing: the only `aria-hidden` element left that owns text is
 * none — `.lwe-meter-track` is an empty bar, and the three `.canvas-legend`
 * paragraphs, which DID own text, were exposed by this sweep.
 */
async function expectNotBlank(page: Page, label: string): Promise<void> {
  const invisible = await page.evaluate(() => {
    const out: string[] = [];
    for (const el of Array.from(document.querySelectorAll('body *'))) {
      const own = Array.from(el.childNodes)
        .filter((n) => n.nodeType === Node.TEXT_NODE)
        .map((n) => n.textContent ?? '')
        .join('')
        .trim();
      if (!own) continue;
      // Deliberately hidden subtrees are not "blank", they are closed.
      if (!(el as HTMLElement).checkVisibility?.({ checkVisibilityCSS: true })) continue;
      if (el.closest('[aria-hidden="true"]')) continue;
      let effective = 1;
      let node: Element | null = el;
      while (node) {
        effective *= parseFloat(getComputedStyle(node).opacity);
        node = node.parentElement;
      }
      if (effective === 0) {
        out.push(`${el.tagName.toLowerCase()}.${(el.getAttribute('class') ?? '').trim()}`);
      }
    }
    return Array.from(new Set(out));
  });
  expect(invisible, `no visible text may render at opacity 0 in state: ${label}`).toEqual([]);
}

/**
 * Uncaught page errors and console errors, collected from the moment the page
 * is created. A renderer that throws halfway through leaves an earlier state on
 * screen, and a gate that scans that state reports green for a page that is
 * broken. Attach before `boot`, assert after the drive.
 *
 * This matters more than usual here. `main.ts` builds the entire page from one
 * `innerHTML` template and then resolves ~40 elements through `byId()`, which
 * THROWS on a miss — so a single renamed id leaves a half-wired page that still
 * looks plausible. The old gate would have scanned it happily.
 */
export function watchPageErrors(page: Page): string[] {
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(`pageerror: ${e.message}`));
  page.on('console', (m) => {
    if (m.type() === 'error') errors.push(`console.error: ${m.text()}`);
  });
  return errors;
}

/**
 * Exactly one banner landmark: the shared bar.
 *
 * The OUTCOME is asserted rather than either mechanism, so a change to the
 * nesting is caught as well as a change to the script. Two mechanisms are at
 * work on this page and they disagree about which one is load-bearing:
 * `index.html`'s `dedupeBanner()` demotes any second banner to `role="group"`,
 * while this lab's own `<header class="cl-hero">` sits inside `<main
 * id="main-content">`, which scopes it out of the banner role on its own — and
 * `dedupeBanner` explicitly returns early for exactly that case
 * (`el.closest('main, article, aside, nav, section')`). So neither actually
 * fires here, and the single banner is a property of the markup.
 */
export async function assertSingleBanner(page: Page): Promise<void> {
  const banners = await page.evaluate(() => {
    const scoped = new Set(['MAIN', 'ARTICLE', 'ASIDE', 'NAV', 'SECTION']);
    const isBanner = (el: Element): boolean => {
      if (el.getAttribute('role') === 'banner') return true;
      if (el.tagName !== 'HEADER') return false;
      if (el.getAttribute('role')) return false; // explicit non-banner role wins
      for (let p = el.parentElement; p; p = p.parentElement) if (scoped.has(p.tagName)) return false;
      return true;
    };
    return [...document.querySelectorAll('header,[role="banner"]')].filter(isBanner).length;
  });
  expect(banners, 'exactly one banner landmark').toBe(1);
}

/**
 * Load the page in a known theme with reduced motion actually in effect, and
 * assert the content every scan relies on is really on the page — including the
 * lab's DEFAULTS, which are never assumed.
 *
 * `test.use({ reducedMotion })` silently does nothing on Playwright 1.61.1, so
 * the emulation is applied imperatively BEFORE the navigation and then
 * *asserted* from inside the page. It has to be before `goto` for a second
 * reason specific to this lab: `main.ts` evaluates
 * `matchMedia('(prefers-reduced-motion: reduce)').matches` ONCE at module scope
 * and caches it in `prefersReducedMotion`. A preference applied after load
 * would change the CSS and not the TypeScript, leaving the basis tween and the
 * 520ms collapse waits live in a page whose stylesheet says motion is off.
 *
 * The theme is seeded through `localStorage` rather than by clicking the
 * toggle, which also pins down a real failure mode: `index.html`'s anti-flash
 * script reads `localStorage.getItem('theme')` and both toggles — the shared
 * bar's `#cl-theme-toggle` and the lab's own `#theme-toggle` — write
 * `localStorage.setItem('theme', …)`. If those keys ever drift apart the theme
 * silently stops persisting, and this boot fails on `data-theme` rather than
 * quietly scanning dark twice.
 *
 * The defaults are asserted at length because two of this lab's five exhibits
 * ship EMPTY (the LWE transcript and the embedding grid have no content until
 * Generate is pressed) while the other three ship POPULATED from a specific set
 * of slider positions and a specific 2x2 matrix — and every ratio the gate
 * measures in Exhibits 1-3 is a ratio of that particular rendering. A drive
 * that assumed the defaults would be measuring whichever half of the lab the
 * defaults happened to select.
 */
export async function boot(page: Page, theme: 'dark' | 'light'): Promise<void> {
  // A click on a control that never becomes actionable otherwise burns the whole
  // test timeout and reports nothing useful. 20s turns that silent hang into a
  // named failure naming the locator.
  page.setDefaultTimeout(20_000);
  await page.emulateMedia({ reducedMotion: 'reduce' });
  await page.addInitScript((t) => localStorage.setItem('theme', t), theme);
  await page.goto('.');
  expect(
    await page.evaluate(() => matchMedia('(prefers-reduced-motion: reduce)').matches),
    'reduced-motion emulation must actually be in effect'
  ).toBe(true);
  await expect(page.locator('html')).toHaveAttribute('data-theme', theme);
  await assertSingleBanner(page);

  // The whole page is one `innerHTML` assignment followed by ~40 `byId()`
  // lookups that throw on a miss, so a navigation that resolves proves nothing.
  await expect(page.locator('main#main-content')).toBeVisible();
  await expect(page.locator('section.exhibit')).toHaveCount(7);

  // The reduced-motion block's one job here, asserted rather than assumed: the
  // six `.exhibit` sections must not be stranded at `rise`'s `opacity: 0`.
  expect(
    await page.evaluate(
      () => getComputedStyle(document.querySelector('#exhibit-1') as Element).opacity
    ),
    'reduced motion must leave .exhibit opaque, not stranded at rise{from{opacity:0}}'
  ).toBe('1');

  // ── Exhibit 1: four sliders and a determinant computed from them ─────────
  await expect(page.locator('#b1-angle')).toHaveValue('25');
  await expect(page.locator('#b1-len')).toHaveValue('5');
  await expect(page.locator('#b2-angle')).toHaveValue('95');
  await expect(page.locator('#b2-len')).toHaveValue('4');
  await expect(page.locator('#det-label')).toHaveText('det(Lambda) = 18.794');
  await expect(page.locator('#same-lattice-label')).toBeEmpty();
  await expect(page.locator('#det-tip')).toBeHidden();
  await expect(page.locator('#det-help')).toHaveAttribute('aria-expanded', 'false');

  // ── Exhibit 2: the 2x2 GSO matrix ────────────────────────────────────────
  for (const [id, v] of [
    ['gso-a11', '3'],
    ['gso-a12', '1'],
    ['gso-a21', '1'],
    ['gso-a22', '2'],
  ] as const) {
    await expect(page.locator(`#${id}`)).toHaveValue(v);
  }
  await expect(page.locator('#gso-text')).not.toBeEmpty();

  // ── Exhibit 3: shipped 2D basis, delta and a reduction already computed ──
  await expect(page.locator('#lll-dim')).toHaveValue('2');
  await expect(page.locator('#lll-delta')).toHaveValue('0.75');
  await expect(page.locator('#lll-matrix')).toHaveValue('19 2\n7 1');
  await expect(page.locator('#lll-log')).toContainText('1.');
  await expect(page.locator('#lll-metrics')).toContainText('Step count: 1 /');

  // ── Exhibit 4: ships entirely empty, and its one control ships DISABLED ──
  await expect(page.locator('#lwe-n')).toHaveValue('4');
  await expect(page.locator('#lwe-q')).toHaveValue('71');
  await expect(page.locator('#lwe-sigma')).toHaveValue('2');
  await expect(page.locator('#lwe-beta')).toHaveValue('2');
  await expect(page.locator('#lwe-output')).toBeEmpty();
  await expect(page.locator('#embed-grid tbody')).toHaveCount(0);
  await expect(page.locator('#embed-collapse')).toBeDisabled();
  await expect(page.locator('#lwe-progress')).toHaveJSProperty('value', 0);
  await expect(page.locator('#lwe-meter-text')).toHaveText(
    'Norm-gap confidence: waiting for attack run.'
  );

  // ── Exhibit 5 ────────────────────────────────────────────────────────────
  await expect(page.locator('#explore-n')).toHaveValue('8');
  await expect(page.locator('#explore-left')).not.toBeEmpty();
  await expect(page.locator('#explore-right')).not.toBeEmpty();
  await expect(page.locator('#beta-tip')).toBeHidden();

  // ── The reproducible-lab controls, and the mode the page ships in ────────
  await expect(page.locator('#lab-seed')).toHaveValue('');
  await expect(page.locator('#lab-preset')).toHaveValue('');
  await expect(page.locator('#lab-status')).toHaveText(
    'Mode: crypto-random (not reproducible). Enter a seed for repeatable runs.'
  );

  // Twelve disclosures, all shut.
  await expect(page.locator('details')).toHaveCount(12);
  await expect(page.locator('details[open]')).toHaveCount(0);

  await settle(page);
  await expectNotBlank(page, `${theme} first paint`);
}

/**
 * Assert that `[hidden]` actually hides.
 *
 * `[hidden]` is a UA rule of specificity (0,0,0) in the author-facing sense —
 * any author `display` declaration beats it, including a single class — so
 * `hidden` can be set on an element and do nothing at all. Seven labs in this
 * fleet shipped that. This lab has two elements it toggles by attribute,
 * `#det-tip` and `#beta-tip`, and this measures the computed `display` rather
 * than trusting that `.det-tip` never grows a `display` declaration.
 */
export async function expectHiddenActuallyHides(page: Page): Promise<void> {
  const leaks = await page.evaluate(() => {
    const out: string[] = [];
    for (const id of ['det-tip', 'beta-tip']) {
      const el = document.getElementById(id);
      if (!el) {
        out.push(`${id}: missing`);
        continue;
      }
      const had = el.hasAttribute('hidden');
      el.setAttribute('hidden', '');
      const display = getComputedStyle(el).display;
      if (!had) el.removeAttribute('hidden');
      if (display !== 'none') out.push(`#${id} computes display:${display} while [hidden]`);
    }
    return out;
  });
  expect(leaks, '[hidden] must actually hide — a class-level display beats it').toEqual([]);
}

/**
 * Assert the page does not require horizontal scrolling.
 *
 * WCAG 1.4.10 (Reflow, AA). axe has no rule for this at all, and this page is
 * the shape that breaks it: a 13-column integer matrix, five fixed-size
 * `<canvas>` elements (the widest declares `width="900"`), monospaced dumps of
 * n-dimensional vectors, and a `.controls-row` grid that carries up to eight
 * tracks. Each wide thing is meant to scroll inside its own box; the assertion
 * here is that none of them scrolls the DOCUMENT.
 */
export async function expectNoHorizontalOverflow(page: Page, label: string): Promise<void> {
  const overflow = await page.evaluate(() => {
    const doc = document.documentElement;
    if (doc.scrollWidth <= doc.clientWidth) return null;

    // Only elements that actually push the DOCUMENT sideways are culprits. A
    // wide box inside an `overflow-x: auto` wrapper has a huge bounding rect but
    // is clipped by its scroller and contributes nothing to the document's
    // scroll width — naming it sends you off fixing the wrong element. That cost
    // a run elsewhere in this fleet, and this page has a decoy behind every
    // `pre.panel` and behind `.embed-grid-wrap`.
    const clipped = (el: Element): boolean => {
      let n = el.parentElement;
      while (n && n !== doc) {
        const ox = getComputedStyle(n).overflowX;
        if (ox === 'auto' || ox === 'scroll' || ox === 'hidden' || ox === 'clip') return true;
        n = n.parentElement;
      }
      return false;
    };

    const over = Array.from(document.querySelectorAll('body *'))
      .map((el) => ({ el, r: el.getBoundingClientRect() }))
      .filter((x) => x.r.width > 0 && x.r.right > doc.clientWidth + 1)
      .sort((a, b) => b.r.right - a.r.right);
    const widest = over.filter((x) => !clipped(x.el))[0] ?? over[0];
    return {
      scrollWidth: doc.scrollWidth,
      clientWidth: doc.clientWidth,
      widest: widest
        ? `${clipped(widest.el) ? '[clipped] ' : ''}${widest.el.tagName.toLowerCase()}${widest.el.id ? '#' + widest.el.id : ''}` +
          `${widest.el.getAttribute('class') ? '.' + widest.el.getAttribute('class')!.trim().split(/\s+/).join('.') : ''}` +
          ` @${Math.round(widest.r.width)}px right=${Math.round(widest.r.right)}`
        : '(none identified)',
    };
  });
  expect(overflow, `page must not scroll horizontally in state: ${label}`).toBeNull();
}

/**
 * Every scrolling container must be operable from the keyboard (WCAG 2.1.1). If
 * it holds no focusable content it needs `tabindex="0"`, so it becomes a focus
 * target arrow keys can then scroll.
 *
 * This lab has seven such boxes and they hold most of its evidence: `#lll-log`
 * (the step-by-step reduction trace), `#lll-metrics` (the determinant
 * invariant), `#gso-text`, `#lwe-output` (the whole attack transcript, the
 * recovered secret and the SUCCESS/FAILED verdict), the two parameter-explorer
 * columns, and `.embed-grid-wrap` around the embedding matrix. Six of the seven
 * only overflow once something has been run into them — `#lll-log` is capped at
 * 420px and grows a line per LLL step, `#lwe-output` is empty at first paint —
 * so this is a check that only has an answer in states the drive has to go and
 * build.
 */
export async function expectScrollersReachable(page: Page, label: string): Promise<void> {
  const unreachable = await page.evaluate(() => {
    const FOCUSABLE = 'a[href],button,input,select,textarea,[tabindex]:not([tabindex="-1"])';
    return Array.from(document.querySelectorAll<HTMLElement>('body *'))
      .filter((el) => el.scrollWidth > el.clientWidth + 1 || el.scrollHeight > el.clientHeight + 1)
      .filter((el) => {
        const cs = getComputedStyle(el);
        return (
          ['auto', 'scroll'].includes(cs.overflowX) || ['auto', 'scroll'].includes(cs.overflowY)
        );
      })
      .filter((el) => el.tabIndex < 0 && !el.querySelector(FOCUSABLE))
      .map(
        (el) =>
          `${el.tagName.toLowerCase()}${el.id ? '#' + el.id : ''}.${(el.getAttribute('class') ?? '').trim()}` +
          ` (${el.scrollWidth}x${el.scrollHeight} in ${el.clientWidth}x${el.clientHeight})`
      );
  });
  expect(
    Array.from(new Set(unreachable)),
    `scrolling regions with no keyboard route in state: ${label}`
  ).toEqual([]);
}

/**
 * Every tab stop must show WHERE the focus is (WCAG 2.4.7).
 *
 * This is here because it is the defect a reflow or 2.1.1 fix CREATES: making a
 * wide panel focusable so it can be scrolled from the keyboard adds a tab stop,
 * and a tab stop with no visible indicator is a new 2.4.7 failure. This lab has
 * seven `tabindex="0"` regions added for exactly that reason
 * (`#gso-text`, `#lll-log`, `#lll-metrics`, `#lwe-output`, `#explore-left`,
 * `#explore-right`, `.embed-grid-wrap`), and the stylesheet's focus rule names
 * `a, button, input, select, textarea` — a `<pre>` or a `<div>` matches none of
 * those five.
 *
 * It walks the REAL tab order with real Tab presses rather than calling
 * `focus()` in a loop, because `:focus-visible` is modality-dependent:
 * programmatic focus on a `<div>` does not match it, so a `focus()`-based check
 * would report a failure for every correctly-styled region and a pass for
 * nothing. `outline-style: auto` counts as an indicator — that is the UA focus
 * ring, which is a real one.
 */
export async function expectFocusVisibleThroughTabOrder(
  page: Page,
  label: string
): Promise<void> {
  // Identity is tracked by ELEMENT, in a page-side array, not by a describe()
  // string: five controls on this page share a description (three
  // `button.btn.preset`, two `a.cl-btn` in the shared bar, six
  // `button.btn.challenge-setup`), and a string-keyed set declares the walk
  // "wrapped" at the first repeat — which reported a 3-stop tab order on a page
  // that has about sixty.
  await page.evaluate(() => {
    (window as unknown as { __tabStops?: Element[] }).__tabStops = [];
    (document.activeElement as HTMLElement | null)?.blur?.();
  });
  const bad = new Set<string>();
  let stops = 0;
  for (let i = 0; i < 200; i += 1) {
    await page.keyboard.press('Tab');
    const stop = await page.evaluate(() => {
      const seen = (window as unknown as { __tabStops: Element[] }).__tabStops;
      const el = document.activeElement as HTMLElement | null;
      // Focus has left the document's focusable set — Tab from the LAST tab
      // stop lands on <body>. That is not the end of the walk: the next Tab
      // re-enters at the top. Returning `null` (rather than stopping) is what
      // lets a walk that starts mid-document — which it does at the end of the
      // drive, because the last thing clicked was a <summary> near the bottom —
      // still reach every stop. Stopping here reported a 4-stop tab order.
      if (!el || el === document.body || el === document.documentElement) return 'edge';
      if (seen.includes(el)) return 'wrapped';
      seen.push(el);
      const cs = getComputedStyle(el);
      const w = parseFloat(cs.outlineWidth || '0');
      const drawn =
        (cs.outlineStyle !== 'none' && (cs.outlineStyle === 'auto' || w > 0)) ||
        (!!cs.boxShadow && cs.boxShadow !== 'none');
      const cls = (el.getAttribute('class') ?? '').trim().split(/\s+/).filter(Boolean).join('.');
      return {
        id: `${el.tagName.toLowerCase()}${el.id ? '#' + el.id : ''}${cls ? '.' + cls : ''}`,
        drawn,
        detail: `outline ${cs.outlineStyle} ${cs.outlineWidth}, box-shadow ${cs.boxShadow}`,
      } as const;
    });
    if (stop === 'edge') continue;
    if (stop === 'wrapped') break;
    stops += 1;
    if (!stop.drawn) bad.add(`${stop.id} — ${stop.detail}`);
  }
  await page.evaluate(() => {
    delete (window as unknown as { __tabStops?: Element[] }).__tabStops;
  });
  expect(stops, `tab order must have stops in state: ${label}`).toBeGreaterThan(10);
  expect(
    Array.from(bad),
    `tab stops with no visible focus indicator in state: ${label}`
  ).toEqual([]);
}

/**
 * When `A11Y_COLLECT` is set, `scan` records failures instead of throwing.
 *
 * A strict gate reports the first failing assertion in the first failing state
 * and stops, so a page with defects in several states needs one full run per
 * defect to enumerate them. The collection pass turns that into a single run. It
 * is a debugging aid only: `A11Y_COLLECT` is never set in CI or in the committed
 * workflow, and a run with it set prints every finding as it happens and then
 * fails at the end, so a green collection run cannot be mistaken for a green
 * gate.
 */
const COLLECTING = !!process.env.A11Y_COLLECT;
const collected: string[] = [];

function record(entry: string): void {
  collected.push(entry);
  // Printed as it happens, not only at the end: a hard assertion later in the
  // drive would otherwise abort the test before anything collected so far was
  // ever shown.
  console.log(`\n[A11Y_COLLECT #${collected.length}] ${entry}`);
}

export function softExpect(actual: unknown, message: string, expected: unknown): void {
  if (!COLLECTING) {
    expect(actual, message).toEqual(expected);
    return;
  }
  try {
    expect(actual, message).toEqual(expected);
  } catch {
    record(`${message}\n  ${JSON.stringify(actual, null, 2)}`);
  }
}

/** Same, for the assertions that live inside an async page probe. */
async function soft(fn: () => Promise<void>): Promise<void> {
  if (!COLLECTING) return fn();
  try {
    await fn();
  } catch (e) {
    record(String(e).slice(0, 900));
  }
}

/**
 * Fail the test if the collection pass recorded anything. Without this a
 * collection run would end green, and a green collection run is
 * indistinguishable from a green gate — which is the exact confusion the whole
 * exercise exists to remove.
 */
export function reportCollected(): void {
  if (!COLLECTING) return;
  expect(collected, `A11Y_COLLECT recorded ${collected.length} failure(s)`).toEqual([]);
}

/**
 * Scan the page as it currently stands.
 *
 * Seven assertions, because axe's `violations` array alone is not a complete
 * oracle. (The eighth, WCAG 2.4.7, is `expectFocusVisibleThroughTabOrder`,
 * which is not called from here because it MOVES focus — walking the whole tab
 * order inside every scan would change the state being scanned. It is driven
 * separately at the two states where this page's tab order differs: first
 * paint, and the fully populated page.)
 *
 *  - reduced-motion end state — see `expectNotBlank`.
 *  - `violations` — the usual WCAG A/AA rule failures, plus four landmark
 *    best-practice rules `withTags` does not run on its own.
 *  - `incomplete` — axe's "could not decide" bucket, which never reaches the
 *    violations array. The one rule id allowed to remain incomplete is
 *    `color-contrast`, and only because the next assertion computes those
 *    ratios arithmetically — which matters more here than in most labs, since
 *    every `.exhibit` is a `linear-gradient` and axe therefore declines to
 *    resolve the surface under any prose on the page. Everything else in that
 *    bucket is a real result axe simply could not finish — including
 *    `aria-prohibited-attr`, which is where an `aria-label` on a role-less
 *    element hides, a defect that never reaches the violations array at all.
 *  - arithmetic contrast — composite-aware WCAG 1.4.3 over every text node.
 *  - non-text contrast — SC 1.4.11 control boundaries AND `::before`/`::after`
 *    generated content, neither of which axe has a rule for and neither of
 *    which the text walk can reach. See `e2e/nontext.ts`.
 *  - keyboard reachability of scrolling regions — WCAG 2.1.1.
 *  - reflow — WCAG 1.4.10, which axe has no rule for at all.
 */
export async function scan(page: Page, label: string): Promise<void> {
  await settle(page);
  await expectNotBlank(page, label);
  // TWO axe runs, deliberately, and this is not a style choice.
  //
  // `AxeBuilder.withTags()` and `AxeBuilder.withRules()` both write `runOnly`,
  // so the second call SILENTLY REPLACES the first — the axe-core/playwright
  // source says so in as many words ("Cannot be used with AxeBuilder#withTags").
  // Chained as `.withTags(TAGS).withRules([...])`, which is the form this gate
  // was copied from, axe ran those best-practice rules and NOT ONE WCAG RULE.
  // Measured on the source repo: the chained form executes 4 rules where
  // `withTags` alone executes 63. A green result meant "no duplicate landmarks"
  // and nothing whatever about WCAG A/AA, while reading like a full pass.
  //
  // Running the two sets separately and merging is the only way to have both.
  // The landmark rules are wanted because they are best-practice rather than
  // WCAG-tagged, so `withTags` alone does not reach them.
  const wcag = await new AxeBuilder({ page }).withTags(TAGS).analyze();
  const landmarks = await new AxeBuilder({ page }).withRules([ 'landmark-no-duplicate-banner', 'landmark-unique', 'landmark-one-main', 'landmark-complementary-is-top-level', ]).analyze();
  const results = {
    violations: [...wcag.violations, ...landmarks.violations],
    incomplete: [...wcag.incomplete, ...landmarks.incomplete],
  };

  const violations = results.violations.map((v) => ({
    state: label,
    id: v.id,
    impact: v.impact,
    help: v.help,
    nodes: v.nodes.map((n) => n.target.join(' ')).slice(0, 8),
  }));
  softExpect(violations, `axe violations in state: ${label}`, []);

  const unexplainedIncomplete = results.incomplete
    .filter((v) => v.id !== 'color-contrast')
    .map((v) => ({
      state: label,
      id: v.id,
      nodes: v.nodes.map((n) => n.target.join(' ')).slice(0, 8),
    }));
  softExpect(unexplainedIncomplete, `axe incomplete results in state: ${label}`, []);

  const contrast = Array.from(new Set(formatContrastFailures(await auditContrast(page))));
  softExpect(contrast, `measured contrast failures in state: ${label}`, []);

  // A 1.4.11 oracle that silently measured nothing would be the same failure
  // this sweep exists to remove, so the population is asserted before the
  // verdict: `#app` must really contain visible controls to judge.
  expect(
    await page.locator('#app button:visible, #app select:visible').count(),
    `no controls found to measure in state: ${label}`
  ).toBeGreaterThan(0);
  softExpect(
    Array.from(new Set(formatNonTextFailures(await auditNonText(page)))),
    `non-text contrast failures (SC 1.4.11 / generated content) in state: ${label}`,
    []
  );

  await soft(() => expectScrollersReachable(page, label));
  await soft(() => expectNoHorizontalOverflow(page, label));
}

// ── The drive ───────────────────────────────────────────────────────────────

/**
 * Open every VISIBLE shut disclosure by clicking its summary.
 *
 * `:visible` is defensive rather than load-bearing — all twelve `<details>` on
 * this page are in the static template and none of them is behind a hidden
 * panel — but forcing `.open` from script is exactly what the gate this
 * replaces did to all twelve, and the point of clicking is that the summary is
 * the control a reader operates.
 */
async function openAllDisclosures(page: Page, expectSome = true): Promise<void> {
  const shut = page.locator('details:not([open]) > summary:visible');
  let opened = 0;
  for (let i = await shut.count(); i > 0 && opened < 40; i = await shut.count()) {
    await shut.first().click();
    opened += 1;
  }
  await expect(page.locator('details:not([open]) > summary:visible')).toHaveCount(0);
  if (expectSome) {
    expect(opened, 'no shut disclosure was found where one was expected').toBeGreaterThan(0);
  }
}

/**
 * Drive the lab through the states that render content, scanning each.
 *
 * Six things shape this drive:
 *
 *  - THE SHIPPED STATE IS SCANNED FIRST, AND IT IS HALF-EMPTY. Exhibits 1, 2, 3
 *    and 5 render from their defaults on mount; Exhibit 4 renders nothing at all
 *    until Generate is pressed, its one auxiliary control ships `disabled`, and
 *    its transcript, meter and 13-column matrix are absent. That asymmetry is
 *    the first thing a reader meets and it is asserted in `boot`.
 *
 *  - EVERY ERROR AND EMPTY STATE IS DRIVEN, NOT JUST THE HAPPY PATH. Two exist
 *    and neither had ever been rendered by a gate: pressing Run Attack before
 *    Generate replaces the transcript with "Generate an LWE instance first.",
 *    and a non-square or non-integer LLL matrix replaces the step log with
 *    "Input error: …". Both are prose the reader is meant to read, in surfaces
 *    nothing else on the page uses.
 *
 *  - THE EXTREMES OF EVERY SLIDER, NOT THE DEFAULTS. `#lll-delta` at 0.999 is
 *    the only route to a long swap cascade, which is the only thing that makes
 *    `#lll-log` (capped at 420px) actually overflow — and whether it can then
 *    be scrolled from a keyboard is a WCAG 2.1.1 question that only exists in a
 *    state a drive has to go and build. `#explore-n` at 256 is the only route to
 *    the Kyber-scale column, and `#lwe-n` at 12 is the only route to the
 *    "n too large for brute force" branch.
 *
 *  - BOTH VERDICTS OF THE ATTACK. Exhibit 4 exists to separate a genuine
 *    lattice break from the brute-force baseline, and the two verdicts are
 *    rendered in different words with a different meter colour: `good` (a green
 *    gradient) for a recovered secret and `bad` (a red one) for a failure. A
 *    gate that only ever drives the seeded success case scans one of them.
 *
 *  - THE COLLAPSE REPLAY'S END STATE. `.emerge-hot` (a green outline on the
 *    secret block) and `.embed-resultrow-done` are the non-colour-plus-colour
 *    cues the whole embedding exhibit exists to deliver, and they only exist
 *    after "Show the short vector emerge" has run to completion.
 *
 *  - NO FIXED TIMEOUTS. Every step has a DOM completion signal — a transcript
 *    gaining a line, a progress value, a button returning from `disabled`, a
 *    step counter — and the drive waits on those.
 */
export async function driveAllStates(page: Page, theme: string): Promise<void> {
  const scanAt = (s: string): Promise<void> => scan(page, `${theme} / ${s}`);

  await expectHiddenActuallyHides(page);
  await scanAt('first paint: Exhibit 4 empty, its collapse button locked, 12 disclosures shut');

  // ── Both skip links, focused ─────────────────────────────────────────────
  await page.evaluate(() => (document.activeElement as HTMLElement | null)?.blur?.());
  await page.keyboard.press('Tab');
  await expect(page.locator('a.cl-skip-link')).toBeFocused();
  await scanAt('shared-header skip link focused');

  // The lab's own skip link is the sixth tab stop: after the shared bar's skip
  // link, brand, Menu, GitHub and theme toggle. Tabbing to it rather than
  // calling focus() is what proves it is reachable at all.
  for (let i = 0; i < 8; i += 1) {
    if (await page.locator('#app a.skip-link').evaluate((el) => el === document.activeElement)) break;
    await page.keyboard.press('Tab');
  }
  await expect(page.locator('#app a.skip-link')).toBeFocused();
  await scanAt("the lab's own skip link focused");

  // Ordered after the two skip-link steps deliberately: this walks the whole
  // tab order with ~60 Tab presses, and the sequential focus navigation
  // starting point it leaves behind is NOT reset by `blur()` — running it first
  // makes "one Tab from the top lands on the skip link" untestable.
  await soft(() => expectFocusVisibleThroughTabOrder(page, `${theme} / first paint`));

  // ── Exhibit 1 ────────────────────────────────────────────────────────────
  await page.click('#det-help');
  await expect(page.locator('#det-tip')).toBeVisible();
  await expect(page.locator('#det-help')).toHaveAttribute('aria-expanded', 'true');
  await scanAt('determinant tooltip revealed');

  // Both sliders at their extremes: b1 and b2 collinear, which is the only
  // input that drives the determinant to zero — a degenerate basis, and the one
  // reading of this exhibit that is not a lattice at all.
  await page.locator('#b1-angle').fill('0');
  await page.locator('#b1-len').fill('10');
  await page.locator('#b2-angle').fill('0');
  await page.locator('#b2-len').fill('1');
  await expect(page.locator('#det-label')).toHaveText('det(Lambda) = 0.000');
  await scanAt('Exhibit 1 sliders at their extremes, determinant collapsed to zero');

  await page.locator('#b2-angle').fill('180');
  await expect(page.locator('#det-label')).toHaveText('det(Lambda) = 0.000');
  await page.locator('#b2-angle').fill('95');
  await page.locator('#b1-angle').fill('25');
  await page.locator('#b1-len').fill('5');
  await page.locator('#b2-len').fill('4');
  await expect(page.locator('#det-label')).toHaveText('det(Lambda) = 18.794');
  await page.click('#same-lattice-btn');
  await expect(page.locator('#same-lattice-label')).toHaveText(
    'Same lattice. Same det. Uglier basis.'
  );
  await scanAt('unimodular re-basis applied, determinant unchanged');

  // ── Exhibit 2: the Gram-Schmidt step-through, and BOTH Lovasz verdicts ───
  // The panel is a three-state cycle (b1*, then mu21 and b2*, then the Lovasz
  // check) that editing any matrix cell rewinds to state 1, so every state has
  // to be reached by clicking Next the right number of times from a known
  // start. `boot` asserted the start.
  await expect(page.locator('#gso-text')).toContainText('Step 1');
  await page.click('#gso-next');
  await expect(page.locator('#gso-text')).toContainText('mu21 =');
  await scanAt('Gram-Schmidt step 2: the mu coefficient and b2*');
  await page.click('#gso-next');
  // The SHIPPED basis fails the check: [3,1],[1,2] gives 0.75*||b1*||^2 = 7.5
  // against ||b2||^2 = 5. That is the verdict this exhibit ships showing, and
  // it was worth asserting rather than assuming — the other one is reachable
  // only from a basis nobody arrives at.
  await expect(page.locator('#gso-text')).toContainText('Condition fails: swap needed.');
  await scanAt('Gram-Schmidt step 3: the Lovasz condition failing, as shipped');

  // The other verdict. The check reduces to 0.75*||b1||^2 vs ||b2||^2, so it
  // holds only when b2 is the longer vector — "Condition holds: advance." is a
  // line of prose no other state on the page renders.
  await page.locator('#gso-a11').fill('1');
  await page.locator('#gso-a12').fill('0');
  await page.locator('#gso-a21').fill('0');
  await page.locator('#gso-a22').fill('3');
  await expect(page.locator('#gso-text')).toContainText('Step 1');
  await page.click('#gso-next');
  await page.click('#gso-next');
  await expect(page.locator('#gso-text')).toContainText('Condition holds: advance.');
  await scanAt('Gram-Schmidt step 3: the Lovasz condition holding');

  await page.locator('#gso-a11').fill('3');
  await page.locator('#gso-a12').fill('1');
  await page.locator('#gso-a21').fill('1');
  await page.locator('#gso-a22').fill('2');
  await expect(page.locator('#gso-text')).toContainText('Step 1');

  // ── Exhibit 3: presets, both dimensions, the extremes of delta ───────────
  for (const preset of ['random2d', 'classic', 'near', 'challenge3d', 'qary'] as const) {
    await page.click(`.preset[data-preset="${preset}"]`);
    await expect(page.locator('#lll-metrics')).toContainText('Step count: 1 /');
    await scanAt(`LLL preset "${preset}" loaded`);
  }

  // delta at its maximum on a Fibonacci-style near-parallel basis is the only
  // route to a long swap cascade, which is the only thing that makes #lll-log
  // overflow its 420px cap and become a 2.1.1 question at all.
  await page.locator('#lll-dim').selectOption('2');
  await page.locator('#lll-matrix').fill('610 377\n377 233');
  await page.locator('#lll-delta').fill('0.999');
  await page.click('#lll-reset');
  await expect(page.locator('#lll-metrics')).toContainText('Step count: 1 /');
  await page.click('#lll-auto');
  // The completion signal the code itself defines: the cursor stops advancing
  // once it reaches the last step, so wait for the counter to settle at N / N.
  await expect(page.locator('#lll-metrics')).toContainText(/Step count: (\d+) \/ \1$/m, {
    timeout: 60_000,
  });
  const logLines = (await page.locator('#lll-log').textContent())?.split('\n').length ?? 0;
  expect(logLines, 'the swap cascade must be long enough to overflow #lll-log').toBeGreaterThan(12);
  await scanAt('LLL auto-run to convergence at delta=0.999, step log overflowing');

  await page.locator('#lll-delta').fill('0.5');
  await page.click('#lll-reset');
  await expect(page.locator('#lll-metrics')).toContainText('Step count: 1 /');
  await scanAt('LLL reset at delta=0.5, the weakest Lovasz bound');

  // The input-error state: a non-square matrix. The only route to it, and the
  // only state where #lll-log carries prose rather than a step trace.
  await page.locator('#lll-matrix').fill('1 2 3\n4 5 6');
  await page.click('#lll-reset');
  await expect(page.locator('#lll-log')).toHaveText('Input error: Matrix must be square');
  await scanAt('LLL input error rendered in the step log');

  await page.locator('#lll-matrix').fill('19 2\n7 1');
  await page.click('#lll-reset');
  await expect(page.locator('#lll-metrics')).toContainText('Step count: 1 /');
  await scanAt('LLL reset back to the shipped basis');

  // ── Exhibit 4: the empty-prerequisite state first ────────────────────────
  await expect(page.locator('#embed-collapse')).toBeDisabled();
  await page.click('#lwe-attack');
  await expect(page.locator('#lwe-output')).toHaveText('Generate an LWE instance first.');
  await scanAt('Run Attack pressed with no instance: the prerequisite message');

  // A seeded instance, so the verdict below is reproducible rather than a coin
  // flip. 1234 / n=4 / q=71 / sigma=2 is the lab's own "toy LWE success" lab.
  await page.locator('#lab-seed').fill('1234');
  await page.click('#lab-apply-seed');
  await expect(page.locator('#lab-status')).toContainText('Mode: seeded (seed=1234)');
  await scanAt('a seed applied, the reproducible-mode status line');

  await page.click('#lwe-generate');
  await expect(page.locator('#lwe-output')).toContainText('Secret s =');
  await expect(page.locator('#embed-grid tbody tr').first()).toBeVisible();
  await expect(page.locator('#embed-collapse')).toBeEnabled();
  await expect(page.locator('#lwe-progress')).toHaveJSProperty('value', 40);
  await scanAt('LWE instance generated, the embedding matrix painted by block');

  // The collapse replay. Under reduced motion every wait() inside it is 0, so
  // the end state is reached synchronously; the button returning from
  // `disabled` is the completion signal the code itself defines.
  await page.click('#embed-collapse');
  await expect(page.locator('#embed-collapse')).toBeEnabled();
  await expect(page.locator('#embed-result-row')).toHaveClass(/embed-resultrow-done/);
  await expect(page.locator('#embed-result-row td.emerge-hot').first()).toBeVisible();
  await expect(page.locator('.embed-hint')).toContainText('collapsed: secret block is -s');
  await scanAt('short vector emerged: secret block outlined, result row marked done');

  await page.click('#lwe-attack');
  await expect(page.locator('#lwe-output')).toContainText('Lattice attack result:');
  await expect(page.locator('#lwe-progress')).toHaveJSProperty('value', 100);
  await expect(page.locator('#lwe-meter-text')).toContainText('Norm-gap confidence:');
  await scanAt('lattice attack run on the seeded toy instance');

  // The failing verdict, and with n = 10 also the "n too large for brute force"
  // baseline branch (the threshold in `main.ts` is n > 8). Both are the red
  // meter, which nothing above paints. beta stays at 2 — the point here is the
  // FAILED verdict at high sigma, and beta=6 is already driven by the
  // "bkz-helps" canonical lab below.
  await page.locator('#lwe-n').fill('10');
  await page.locator('#lwe-sigma').fill('10');
  await page.click('#lwe-generate');
  await expect(page.locator('#lwe-output')).toContainText('Secret s =');
  await page.click('#lwe-attack');
  await expect(page.locator('#lwe-output')).toContainText('Lattice attack result: FAILED');
  await expect(page.locator('#lwe-output')).toContainText(
    'Teaching baseline: skipped (n too large for brute force)'
  );
  await scanAt('attack FAILED at n=10, sigma=10 — the red norm-gap meter');

  await page.click('#kyber-try');
  await expect(page.locator('#lwe-output')).toContainText('ML-KEM-512');
  await expect(page.locator('#lwe-meter-text')).toContainText('Kyber-scale regime');
  await scanAt('Kyber-512 parameters appended to the transcript');

  // ── Exhibit 5 ────────────────────────────────────────────────────────────
  await page.click('#beta-help');
  await expect(page.locator('#beta-tip')).toBeVisible();
  await scanAt('beta-formula tooltip revealed');

  // The slider drives the RIGHT column only — the left one is a fixed toy
  // reference clamped at n<=16, which is why it is not the thing asserted here.
  await page.locator('#explore-n').fill('256');
  await expect(page.locator('#explore-right')).toContainText('n=256, q=3329');
  await expect(page.locator('#explore-right')).toContainText('out of reach');
  await scanAt('parameter explorer at n=256, the Kyber-scale column');
  // n > 20 flips the verdict block but stays under the two accent thresholds
  // the slider paints itself with (<50 red, <100 amber, else green), so this is
  // the only state that renders the middle one.
  await page.locator('#explore-n').fill('60');
  await expect(page.locator('#explore-right')).toContainText('n=60, q=3329');
  await scanAt('parameter explorer at n=60, the middle slider accent');
  await page.locator('#explore-n').fill('4');
  await expect(page.locator('#explore-right')).toContainText('n=4, q=3329');
  await expect(page.locator('#explore-right')).toContainText('Status: still toy scale');
  await scanAt('parameter explorer at its minimum n');

  // ── The canonical labs, which re-drive several exhibits at once ──────────
  for (const lab of [
    'swap-cascade',
    'delta-sensitivity',
    'toy-lwe-success',
    'bkz-helps',
    'kyber-failure',
  ] as const) {
    await page.locator('#lab-preset').selectOption(lab);
    await expect(page.locator('#lab-status')).toContainText('Mode: seeded');
    await scanAt(`canonical lab "${lab}" loaded`);
  }

  // The challenge buttons drive the same machinery through a different route —
  // they set the <select> and dispatch its change event — so one is enough to
  // prove the route, and it is the one that lands in a different exhibit.
  await page.click('.challenge-setup[data-lab="bkz-helps"]');
  await expect(page.locator('#lab-preset')).toHaveValue('bkz-helps');
  await expect(page.locator('#lwe-output')).toContainText('Secret s =');
  await scanAt('challenge 6 set up through its own button');

  await page.click('#lab-copy-link');
  await expect(page.locator('#lab-status')).toContainText('Lab link copied to clipboard.');
  await scanAt('lab link copied, confirmation in the status line');

  // ── Everything open ──────────────────────────────────────────────────────
  await openAllDisclosures(page);
  await expect(page.locator('details[open]')).toHaveCount(12);
  await scanAt('the finished page with all twelve disclosures open');
  // The tab order is at its longest here — every disclosure open, every output
  // region populated — so this is the state where a missing focus indicator on
  // a scroll region is reachable.
  await soft(() => expectFocusVisibleThroughTabOrder(page, `${theme} / fully populated`));
}
