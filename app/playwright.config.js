import { defineConfig } from "@playwright/test";
export default defineConfig({
  testDir: "./test/browser",
  use: {
    baseURL: process.env.APP_URL || "http://127.0.0.1:5173",
    viewport: { width: 1440, height: 1000 },
    launchOptions: {
      executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
      args: [
        "--no-sandbox",
        "--use-gl=angle",
        "--use-angle=swiftshader",
        "--enable-unsafe-swiftshader",
      ],
    },
  },
  workers: 1,
  timeout: 120000,
});
