import { defineConfig } from '@playwright/test';

/**
 * Accessibility gate. Unlike the smoke config (which runs against the Vite dev
 * server), this scans the *production build* served by `vite preview`, so what
 * passes here is exactly what ships to GitHub Pages. Run `npm run build` first
 * (CI does). Single Chromium project; the spec drives both themes itself.
 */
const PORT = 4372;
const BASE = '/crypto-lab-lll-break/';

export default defineConfig({
  testDir: './e2e',
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 1 : 0,
  reporter: process.env.CI ? 'line' : 'list',
  use: {
    baseURL: `http://localhost:${PORT}${BASE}`,
    // The gate runs in the default (dark) theme; the light-theme test toggles.
    colorScheme: 'dark',
    trace: 'on-first-retry',
  },
  projects: [
    {
      name: 'a11y-chromium',
      use: {},
    },
  ],
  webServer: {
    command: `npm run preview -- --port ${PORT} --strictPort`,
    url: `http://localhost:${PORT}${BASE}`,
    reuseExistingServer: !process.env.CI,
    timeout: 120_000,
  },
});
