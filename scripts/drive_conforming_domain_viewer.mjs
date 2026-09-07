// Drive the built workbench in headless Chromium and record what the conforming-domain
// views actually render, plus every console message and failed request.
//
//   node scripts/drive_conforming_domain_viewer.mjs [--url http://127.0.0.1:8765] [--out DIR]
import { createRequire } from "node:module";
import { mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const require = createRequire(resolve(root, "app/package.json"));
const { chromium } = require("playwright");

const argv = process.argv.slice(2);
const option = (name, fallback) => {
  const index = argv.indexOf(name);
  return index >= 0 ? argv[index + 1] : fallback;
};
const url = option("--url", "http://127.0.0.1:8765");
const out = resolve(option("--out", "test-results/conforming-domains"));
mkdirSync(out, { recursive: true });

const browser = await chromium.launch({
  executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
  args: ["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
});
const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
const console_ = [];
const failed = [];
page.on("console", (m) => console_.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => console_.push(`[pageerror] ${e.message}`));
page.on("requestfailed", (r) => failed.push(`${r.url()} :: ${r.failure()?.errorText}`));
page.on("response", (r) => { if (r.status() >= 400) failed.push(`${r.url()} :: HTTP ${r.status()}`); });

const report = { url, views: [] };
await page.goto(url, { waitUntil: "load" });
await page.waitForFunction(() => document.getElementById("model")?.value === "ihm-body", null, { timeout: 120000 });
await page.waitForSelector("#regional-study option[value^='conforming-domain-']", { state: "attached", timeout: 60000 });
// The Experiments pane is not one of the default monitor cards; add it the way a user would.
if (!(await page.locator("#regional-study").isVisible())) {
  await page.locator("button.new-pane").click();
  await page.getByRole("button", { name: "Add Experiments" }).click();
}
await page.waitForSelector("#regional-study", { state: "visible", timeout: 120000 });
report.options = await page.locator("#regional-study option").evaluateAll((os) =>
  os.map((o) => ({ value: o.value, label: o.textContent })));

for (const kind of (option("--domains", "hand-reflexive-grip-coarse,whole-body-0.01m").split(","))) {
  const value = "conforming-domain-" + kind;
  const started = Date.now();
  await page.selectOption("#regional-study", value);
  await page.waitForFunction(
    (v) => document.getElementById("render-count")?.textContent.includes("boundary triangles"),
    value, { timeout: 180000 });
  await page.waitForTimeout(2500);
  const view = {
    kind,
    load_seconds: (Date.now() - started) / 1000,
    view_title: await page.locator("#view-title").textContent(),
    render_count: await page.locator("#render-count").textContent(),
    regional_note: await page.locator("#regional-note").textContent(),
    flow_legend: await page.locator("#flow-legend").textContent(),
    roles: await page.locator("#domain-roles label").allTextContents(),
    detail_head: await page.locator("#details").textContent(),
    // What is actually on the GPU, read back from the live three.js scene.
    drawn: await page.evaluate(() => {
      const canvas = document.querySelector("#viewport canvas");
      const context = canvas?.getContext("webgl2") || canvas?.getContext("webgl");
      return { canvas: canvas ? [canvas.width, canvas.height] : null, context: !!context };
    }),
    // Non-background pixel fraction inside the viewport tells us something was rasterised.
    coverage: await page.evaluate(() => {
      const canvas = document.querySelector("#viewport canvas");
      if (!canvas) return null;
      const copy = document.createElement("canvas");
      copy.width = 320; copy.height = 200;
      const ctx = copy.getContext("2d");
      ctx.drawImage(canvas, 0, 0, copy.width, copy.height);
      const data = ctx.getImageData(0, 0, copy.width, copy.height).data;
      const counts = new Map();
      for (let i = 0; i < data.length; i += 4) {
        const key = `${data[i] >> 4},${data[i + 1] >> 4},${data[i + 2] >> 4}`;
        counts.set(key, (counts.get(key) || 0) + 1);
      }
      const total = data.length / 4;
      const sorted = [...counts.entries()].sort((a, b) => b[1] - a[1]);
      return { distinct_colors: counts.size, background_fraction: sorted[0][1] / total,
               top: sorted.slice(0, 5).map(([k, n]) => [k, +(n / total).toFixed(4)]) };
    }),
  };
  await page.screenshot({ path: `${out}/${kind}.png` });
  // Click the middle of the viewport to exercise owner identification.
  const box = await page.locator("#viewport").boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.waitForTimeout(400);
  view.picked = await page.locator("#domain-selected").textContent();
  // Reveal the enclosing owner too, then capture again.
  const envelope = page.locator('#domain-roles input[data-role="body_envelope"], #domain-roles input[data-role="skin"]').first();
  if (await envelope.count()) {
    await envelope.click();
    await page.waitForTimeout(1200);
    await page.screenshot({ path: `${out}/${kind}-enclosing.png` });
    view.enclosing_shown = true;
  }
  report.views.push(view);
}
await page.selectOption("#regional-study", "body");
await page.waitForTimeout(1500);
await page.screenshot({ path: `${out}/back-to-body.png` });
report.console = console_;
report.failed_requests = failed;
await browser.close();
console.log(JSON.stringify(report, null, 2));
