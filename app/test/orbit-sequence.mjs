// Screenshot an orbit sequence around one body plane, one frame per drag, so a
// stuck arc is visible rather than argued about.
//
//   node app/test/orbit-sequence.mjs [--url http://127.0.0.1:8765] [--tag after]
//
// Two sweeps, both starting from the coronal snap: a vertical drag travels the
// coronal plane (up over the head, down the back), a horizontal drag travels
// the transverse plane (around the body). Each frame records where the camera
// actually ended up, so the contact sheet and the numbers agree.
import { createRequire } from "node:module";
import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "../..");
const require = createRequire(resolve(root, "app/package.json"));
const { chromium } = require("playwright");

const argv = process.argv.slice(2);
const option = (name, fallback) => {
  const index = argv.indexOf(name);
  return index >= 0 ? argv[index + 1] : fallback;
};
const url = option("--url", "http://127.0.0.1:8765");
const tag = option("--tag", "after");
const out = resolve(option("--out", `test-results/orbit-sequence/${tag}`));
mkdirSync(out, { recursive: true });

const deg = (r) => (r * 180) / Math.PI;
const angle = (a, b) => deg(Math.acos(Math.min(1, Math.max(-1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]))));

const browser = await chromium.launch({
  executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
  args: ["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const console_messages = [], network_errors = [];
page.on("console", (m) => console_messages.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => console_messages.push(`[pageerror] ${e.message}`));
page.on("requestfailed", (r) => network_errors.push(`${r.url()} :: ${r.failure()?.errorText}`));
page.on("response", (r) => { if (r.status() >= 400) network_errors.push(`${r.url()} :: HTTP ${r.status()}`); });

await page.goto(url, { waitUntil: "load" });
await page.waitForSelector("#materialization", { timeout: 120000 });
await page.waitForFunction(() => document.getElementById("scene-status")?.textContent === "", null, { timeout: 240000 });
await page.waitForFunction(() => !!globalThis.__ihmCamera, null, { timeout: 30000 });
await page.waitForTimeout(2500);

const box = await page.locator("#scene").boundingBox();
const centre = { x: box.x + box.width * 0.5, y: box.y + box.height * 0.5 };
const pose = () => page.evaluate(() => globalThis.__ihmCamera.pose());
// The same update the render loop runs, repeated until damping is spent: under
// swiftshader wall-clock waiting measures the rasteriser, not the orbit.
const stop = () => page.evaluate(() => globalThis.__ihmCamera.settle(900));
async function drag(dx, dy) {
  await page.mouse.move(centre.x, centre.y);
  await page.mouse.down();
  await page.mouse.move(centre.x + dx, centre.y + dy, { steps: 16 });
  await page.mouse.up();
  return stop();
}

const report = { url, tag, sweeps: [] };
// 12 frames of 84 px is 12 * 30.24 = 362.9 degrees on a 1000 px canvas: one
// complete lap, if the orbit lets it be completed.
const STEP = 84, FRAMES = 12;
for (const [name, dx, dy] of [["coronal", 0, -STEP], ["transverse", STEP, 0]]) {
  await page.evaluate(() => globalThis.__ihmCamera.snap("coronal"));
  let previous = await stop();
  const frames = [{ frame: 0, direction: previous.direction.map((v) => +v.toFixed(3)), from_start_deg: 0, file: `${name}-00.png` }];
  await page.screenshot({ path: resolve(out, `${name}-00.png`) });
  const start = previous;
  let travelled = 0;
  for (let i = 1; i <= FRAMES; i++) {
    const now = await drag(dx, dy);
    travelled += angle(previous.direction, now.direction);
    previous = now;
    const file = `${name}-${String(i).padStart(2, "0")}.png`;
    await page.screenshot({ path: resolve(out, file) });
    frames.push({ frame: i, direction: now.direction.map((v) => +v.toFixed(3)), from_start_deg: +angle(start.direction, now.direction).toFixed(2), travelled_deg: +travelled.toFixed(2), file });
  }
  report.sweeps.push({ plane: name, drag_px: [dx, dy], requested_deg: +((360 * Math.hypot(dx, dy) * FRAMES) / 1000).toFixed(1), travelled_deg: +travelled.toFixed(2), frames });
}
report.console_messages = console_messages;
report.network_errors = network_errors;
writeFileSync(resolve(out, "orbit-sequence.json"), JSON.stringify(report, null, 2));
await browser.close();
console.log(JSON.stringify(report.sweeps.map((s) => ({ plane: s.plane, requested_deg: s.requested_deg, travelled_deg: s.travelled_deg, path: s.frames.map((f) => f.direction) })), null, 1));
console.log("console:", JSON.stringify(console_messages));
console.log("network:", JSON.stringify(network_errors));
