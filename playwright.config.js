import { defineConfig, devices } from '@playwright/test'

const PORT = 8139

export default defineConfig({
  testDir: 'tests/e2e',
  forbidOnly: !!process.env.CI,
  reporter: process.env.CI ? 'github' : 'list',
  use: { baseURL: `http://localhost:${PORT}` },
  projects: [{ name: 'chromium', use: { ...devices['Desktop Chrome'] } }],
  webServer: { command: 'node tests/e2e/serve.mjs', port: PORT, env: { PORT: String(PORT) }, reuseExistingServer: false },
})
