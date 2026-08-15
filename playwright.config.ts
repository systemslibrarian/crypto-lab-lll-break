import { defineConfig, devices } from '@playwright/test';

// Smoke + accessibility tests run against the *production build* served by
// `vite preview`, which Playwright builds, starts and stops automatically, so
// what passes here is exactly what ships to GitHub Pages. This used to run the
// Vite dev server: that bundle differs from the shipped one in minification,
// base-path handling and asset URLs, so a green run proved nothing about it.
// The port is unique across the crypto-lab fleet — with `reuseExistingServer`
// a shared port lets the suite silently scan a *different* lab's preview.
const PORT = 4606;
// Vite serves under a GitHub Pages base path; tests must hit it, not '/'.
const BASE = '/crypto-lab-lll-break/';

export default defineConfig({
  testDir: './tests',
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 1 : 0,
  reporter: process.env.CI ? 'line' : 'list',
  use: {
    baseURL: `http://localhost:${PORT}${BASE}`,
    trace: 'on-first-retry',
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] },
    },
  ],
  webServer: {
    // Build before serving. `preview` only serves whatever is already in
    // dist/; without the build in front, a failing build leaves the previous
    // good bundle on disk and the gate passes green against code that no
    // longer compiles — silently invalidating mutation checks.
    command: `npm run build && npm run preview -- --port ${PORT} --strictPort`,
    url: `http://localhost:${PORT}${BASE}`,
    reuseExistingServer: !process.env.CI,
    timeout: 120_000,
  },
});
