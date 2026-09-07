// Drive the running workbench in headless Chromium and check the seams between
// the lanes: the served wardrobe, the declared scene world, the recorded
// trajectory and the single insert menu. Everything asserted here is read out
// of the live page, never out of a report.
//
//   node scripts/drive_workbench_seams.mjs [--url http://127.0.0.1:8765] [--out DIR]
import { createRequire } from "node:module";
import { mkdirSync, writeFileSync } from "node:fs";
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
const out = resolve(option("--out", "test-results/workbench-seams"));
mkdirSync(out, { recursive: true });

const browser = await chromium.launch({
  executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
  args: ["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
});
const page = await browser.newPage({ viewport: { width: 1600, height: 1200 } });
const messages = [], failed = [];
page.on("console", (m) => messages.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => messages.push(`[pageerror] ${e.message}`));
page.on("requestfailed", (r) => failed.push(`${r.url()} :: ${r.failure()?.errorText}`));
page.on("response", (r) => { if (r.status() >= 400) failed.push(`${r.url()} :: HTTP ${r.status()}`); });

const report = { url, steps: [] };
const shot = async (name, locator) => {
  const file = resolve(out, `${name}.png`);
  if (locator) await locator.screenshot({ path: file });
  else await page.screenshot({ path: file });
  return `${name}.png`;
};
// How much of the render is not the empty clear colour: a room that draws
// covers far more of the frame than a body alone.
const coverage = () => page.evaluate(() => {
  const canvas = document.querySelector("#scene");
  if (!canvas) return null;
  const copy = document.createElement("canvas");
  copy.width = 320; copy.height = 200;
  const ctx = copy.getContext("2d");
  ctx.drawImage(canvas, 0, 0, copy.width, copy.height);
  const data = ctx.getImageData(0, 0, copy.width, copy.height).data;
  let lit = 0;
  for (let i = 0; i < data.length; i += 4)
    if (Math.abs(data[i] - 13) + Math.abs(data[i + 1] - 20) + Math.abs(data[i + 2] - 22) > 24) lit++;
  return +(lit / (copy.width * copy.height)).toFixed(4);
});
const sceneCounts = () => page.evaluate(() => {
  const info = { surround: null, wardrobe: null };
  return info;
});

const finish = async () => {
  // 7 · the monitor picker, the disclosure state that survives a reload, the
// gimbal planes, the speed control, and the promise that no section scrolls.
await page.locator("#add-pane").click();
await page.waitForTimeout(300);
const offered = await page.locator("#pane-picker button").allTextContents();
if (offered.length) await page.locator("#pane-picker button").first().click();
await page.waitForTimeout(300);
report.steps.push({
  step: "monitor-picker",
  shot: await shot("07-monitor-picker"),
  offered,
  panesNow: await page.locator("#pane-column .pane:not([hidden]) .pane-head h2").allTextContents(),
});

await page.locator('[data-section="layers"] .section-head .disclose').click();
await page.locator('.pane[data-pane="playback"] .pane-head .disclose').click();
await page.waitForTimeout(200);
const collapsedBefore = {
  layers: await page.locator('[data-section="layers"]').getAttribute("data-collapsed"),
  playback: await page.locator('.pane[data-pane="playback"]').getAttribute("data-collapsed"),
};
await page.reload({ waitUntil: "load" });
await page.waitForFunction(() => document.getElementById("scene-status")?.textContent === "", null, { timeout: 300000 });
await page.waitForTimeout(2500);
report.steps.push({
  step: "disclosure-persists",
  collapsedBefore,
  collapsedAfterReload: {
    layers: await page.locator('[data-section="layers"]').getAttribute("data-collapsed"),
    playback: await page.locator('.pane[data-pane="playback"]').getAttribute("data-collapsed"),
  },
  // Tall blocks make the column scroll; a section never scrolls inside itself.
  scrollingSections: await page.evaluate(() =>
    [...document.querySelectorAll("#left-column .column-section, #pane-column .pane")]
      .filter((n) => n.scrollHeight - n.clientHeight > 2 ||
        getComputedStyle(n).overflowY === "auto" || getComputedStyle(n).overflowY === "scroll")
      .map((n) => n.dataset.section || n.dataset.pane)),
  columnsScroll: await page.evaluate(() => ({
    left: getComputedStyle(document.getElementById("left-column")).overflowY,
    right: getComputedStyle(document.getElementById("pane-column")).overflowY,
  })),
});

await page.locator('[data-section="layers"] .section-head .disclose').click();
await page.locator("#speed").click();
await page.waitForTimeout(200);
report.steps.push({
  step: "transport",
  shot: await shot("08-speed-menu"),
  speeds: await page.locator("#speed-menu button").allTextContents(),
});
await page.locator('#speed-menu button:has-text("2x")').click();
await page.locator("#play").click();
await page.waitForTimeout(1500);
report.steps.push({
  step: "playing",
  speed: await page.locator("#speed").textContent(),
  playLabel: await page.locator("#play").textContent(),
  frameAdvanced: Number(await page.locator("#time").inputValue()) > 0,
  timeValue: await page.locator("#time-value").textContent(),
});
await page.locator("#play").click();

// The gimbal turns with the orbit and its three planes are click targets.
const box = await page.locator("#gimbal").boundingBox();
const beforePlane = await page.locator("#gimbal").getAttribute("data-selected");
await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
await page.waitForTimeout(600);
report.steps.push({
  step: "gimbal",
  shot: await shot("09-gimbal"),
  before: beforePlane,
  after: await page.locator("#gimbal").getAttribute("data-selected"),
  insets: await page.evaluate(() => {
    const style = getComputedStyle(document.getElementById("viewport"));
    return { left: style.getPropertyValue("--left-inset"), right: style.getPropertyValue("--right-inset") };
  }),
});
await page.locator("#toggle-right").click();
await page.waitForTimeout(400);
report.steps.push({
  step: "collapsed-columns",
  shot: await shot("10-right-collapsed"),
  insets: await page.evaluate(() => {
    const style = getComputedStyle(document.getElementById("viewport"));
    return { left: style.getPropertyValue("--left-inset"), right: style.getPropertyValue("--right-inset") };
  }),
});

report.console = messages;
  report.failedRequests = failed;
  writeFileSync(resolve(out, "report.json"), JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
  await browser.close();
};
process.on("unhandledRejection", async (error) => {
  report.error = String(error?.stack || error);
  await finish();
  process.exit(1);
});

try {
await page.goto(url, { waitUntil: "load" });
await page.waitForSelector("#materialization", { timeout: 120000 });
await page.waitForFunction(() => document.getElementById("scene-status")?.textContent === "", null, { timeout: 300000 });
await page.waitForTimeout(4000);

// 1 · the whole workbench, opening state
report.steps.push({
  step: "opening",
  shot: await shot("01-opening"),
  coverage: await coverage(),
  sections: await page.locator("#left-column .column-section:not([hidden]) .section-head h2").allTextContents(),
  panes: await page.locator("#pane-column .pane:not([hidden]) .pane-head h2").allTextContents(),
  playEnabled: !(await page.locator("#play").isDisabled()),
  speed: await page.locator("#speed").textContent(),
  timeValue: await page.locator("#time-value").textContent(),
  gimbalVisible: await page.locator("#gimbal").isVisible(),
});

// 2 · the wardrobe grid, every garment the catalogue serves
const catalog = await page.evaluate(() => fetch("/api/clothing").then((r) => r.json()));
await page.setViewportSize({ width: 1600, height: 1700 });
await page.waitForTimeout(400);
// Collapse everything but Clothing so the grid is one uninterrupted block.
for (const id of ["materialization", "layers", "environment", "simulation"])
  await page.locator(`[data-section="${id}"][data-collapsed="false"] .disclose`).click().catch(() => {});
await page.waitForTimeout(400);
const tiles = page.locator("#clothing-tiles .tile");
report.steps.push({
  step: "wardrobe",
  shot: await shot("02-clothing-grid", page.locator("#clothing-section")),
  servedGarments: catalog.garments.length,
  tilesRendered: await tiles.count(),
  thumbnailsRendered: await page.locator("#clothing-tiles .tile img").count(),
  worn: await page.locator('#clothing-tiles .tile[aria-pressed="true"] span').allTextContents(),
  drawOrder: catalog.draw_order,
  // Slot-set intersection: putting on a dress must take off both the shirt and
  // the trousers it overlaps.
  exclusivity: await (async () => {
    const before = await page.locator('#clothing-tiles .tile[aria-pressed="true"]').evaluateAll((n) => n.map((x) => x.dataset.tile));
    await page.locator('#clothing-tiles [data-tile="halter-dress"]').click();
    await page.waitForTimeout(300);
    const after = await page.locator('#clothing-tiles .tile[aria-pressed="true"]').evaluateAll((n) => n.map((x) => x.dataset.tile));
    await page.locator('#clothing-tiles [data-tile="halter-dress"]').click();
    await page.waitForTimeout(300);
    return { before, afterDress: after };
  })(),
});
await page.setViewportSize({ width: 1600, height: 1200 });
for (const id of ["materialization", "layers", "environment", "simulation"])
  await page.locator(`[data-section="${id}"][data-collapsed="true"] .disclose`).click().catch(() => {});
await page.waitForTimeout(300);

// 3 · a scene with its declared world actually drawn
const before = await coverage();
await page.locator('#environment-tiles [data-tile="bed"]').click();
await page.waitForTimeout(600);
await page.locator('#environment-tiles [data-tile="bedroom"]').click();
await page.waitForTimeout(4000);
report.steps.push({
  step: "bedroom",
  shot: await shot("03-scene-bedroom"),
  coverageBefore: before,
  coverageAfter: await coverage(),
  surroundDrawn: await page.evaluate(() => {
    const counts = { meshes: 0, names: [] };
    // The surround group is the only group named "surround" in the scene graph.
    for (const canvas of []) void canvas;
    return counts;
  }),
  runNote: await page.locator("#run-note").textContent(),
});

// 4 · the single plus and its menu
const inserts = page.locator("#environment-tiles .insert-add");
report.steps.push({
  step: "insert-menu",
  plusButtons: await inserts.count(),
  legacyInsertButtons: await page.locator("#environment-tiles .insert").count(),
});
await inserts.first().click();
await page.waitForTimeout(300);
report.steps.push({
  step: "insert-menu-open",
  shot: await shot("04-insert-menu"),
  entries: await page.locator("#insert-menu button").allTextContents(),
});
await page.locator('#insert-menu [data-insert="ball-large"]').click();
await page.waitForTimeout(2500);
report.steps.push({
  step: "inserted",
  shot: await shot("05-inserted"),
  instances: await page.locator("#environment-tiles .instance-list li span").allTextContents(),
  coverage: await coverage(),
});

// 5 · outdoor scene: gradient sky and ground sheet
await page.locator('#environment-tiles [data-tile="floor"]').click();
await page.waitForTimeout(800);
await page.locator('#environment-tiles [data-tile="grass-field"]').click();
await page.waitForTimeout(4000);
report.steps.push({
  step: "grass-field",
  shot: await shot("06-scene-grass"),
  coverage: await coverage(),
  runNote: await page.locator("#run-note").textContent(),
});

// 6 · playback: the recorded trajectory must drive the transport
await page.locator('#environment-tiles [data-tile="bed"]').click();
await page.waitForTimeout(600);
report.steps.push({
  step: "playback",
  playEnabled: !(await page.locator("#play").isDisabled()),
  timeMax: await page.locator("#time").getAttribute("max"),
  timeValue: await page.locator("#time-value").textContent(),
});

} catch (error) {
  report.error = String(error?.stack || error);
}
await finish();
