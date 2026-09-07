// Drive the rebuilt workbench in headless Chromium and record what the new
// left column, floating panes, gimbal and transport actually render, plus every
// console message and failed request.
//
//   node scripts/drive_workbench_ui.mjs [--url http://127.0.0.1:8765] [--out DIR]
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
const out = resolve(option("--out", "test-results/workbench-ui"));
mkdirSync(out, { recursive: true });

const browser = await chromium.launch({
  executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
  args: ["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
});
const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
const messages = [], failed = [];
page.on("console", (m) => messages.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => messages.push(`[pageerror] ${e.message}`));
page.on("requestfailed", (r) => failed.push(`${r.url()} :: ${r.failure()?.errorText}`));
page.on("response", (r) => { if (r.status() >= 400) failed.push(`${r.url()} :: HTTP ${r.status()}`); });

const report = { url, steps: [] };
const shot = async (name) => {
  await page.screenshot({ path: resolve(out, `${name}.png`) });
  return `${name}.png`;
};
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

await page.goto(url, { waitUntil: "load" });
await page.waitForSelector("#materialization", { timeout: 120000 });
await page.waitForFunction(() => document.getElementById("scene-status")?.textContent === "", null, { timeout: 240000 });
await page.waitForTimeout(2500);

// 1 · whole body
report.steps.push({
  step: "whole-body",
  materialization: await page.locator("#materialization").inputValue(),
  materializations: await page.locator("#materialization option").allTextContents(),
  sections: await page.locator("#left-column .column-section h2").allTextContents(),
  clothing_visible: await page.locator("#clothing-section").isVisible(),
  environment_visible: await page.locator("#environment-section").isVisible(),
  clothing_tiles: await page.locator("#clothing-tiles .tile span").allTextContents(),
  environment_slots: await page.locator("#environment-tiles .tile-slot").allTextContents(),
  environment_tiles: await page.locator("#environment-tiles .tile").evaluateAll((n) =>
    n.map((x) => ({ id: x.dataset.tile, slot: x.dataset.slot, on: x.getAttribute("aria-pressed") === "true", available: !x.disabled }))),
  systems: await page.locator("#layers .layer-row > label > span").allTextContents(),
  dynamics: await page.locator("#dynamics .check-row span").allTextContents(),
  gimbal: await page.locator("#gimbal").evaluate((el) => ({
    tag: el.tagName.toLowerCase(), engine: el.dataset.engine || null, selected: el.dataset.selected || null,
    box: [el.clientWidth, el.clientHeight],
  })),
  insets: await page.evaluate(() => {
    const style = getComputedStyle(document.getElementById("viewport"));
    return { left: style.getPropertyValue("--left-inset").trim(), right: style.getPropertyValue("--right-inset").trim() };
  }),
  transport_speed: await page.locator("#speed").textContent(),
  coverage: await coverage(),
  screenshot: await shot("01-whole-body"),
});

// 2 · layers: expand a system and reveal its members
await page.locator("#layer-search").fill("lung");
await page.waitForTimeout(400);
report.steps.push({
  step: "layer-search",
  visible_systems: await page.locator("#layers .layer-row > label > span").allTextContents(),
  members: (await page.locator("#layers .members .member-name").allTextContents()).slice(0, 10),
  screenshot: await shot("02-layer-search"),
});
await page.locator("#layer-search").fill("");
await page.locator("#layers .disclose").first().click();
await page.waitForTimeout(300);
report.steps.push({
  step: "layer-expand",
  members: (await page.locator("#layers .members .member-name").allTextContents()).slice(0, 6),
  screenshot: await shot("03-layer-expand"),
});
await page.locator("#layers .disclose").first().click();

// 3 · click-to-inspect over busy geometry, pane floating on the canvas
const box = await page.locator("#scene").boundingBox();
await page.mouse.click(box.x + box.width * 0.5, box.y + box.height * 0.42);
await page.waitForTimeout(2500);
report.steps.push({
  step: "click-inspect",
  selection: (await page.locator("#details").textContent()).slice(0, 400),
  provenance: await page.locator("#details .structure-provenance").count(),
  screenshot: await shot("04-selection-pane"),
});

// 2b · the three interaction kinds: tiles for environments and scenes,
// contextual dropdowns for configuration, insert actions for objects.
const envState = () => page.evaluate(() => ({
  tiles: [...document.querySelectorAll("#environment-tiles .tile")].map((t) => ({
    id: t.dataset.tile, slot: t.dataset.slot,
    on: t.getAttribute("aria-pressed") === "true", available: !t.disabled })),
  configuration: [...document.querySelectorAll("#environment-tiles select")].map((s) => ({
    slot: s.dataset.slot, value: s.value,
    options: [...s.options].map((o) => o.value).filter(Boolean) })),
  inserts: [...document.querySelectorAll("#environment-tiles [data-insert]")].map((b) => b.dataset.insert),
  instances: [...document.querySelectorAll("#environment-tiles .instance-list li")].map((l) => l.dataset.instance),
}));
report.steps.push({ step: "environment-default", ...(await envState()), screenshot: await shot("02b-environment-default") });

await page.locator('#environment-tiles select[data-slot="bed_support_model"]').selectOption("bed-support-skin-quadrature");
await page.waitForTimeout(300);
report.steps.push({ step: "environment-configuration-appears", ...(await envState()) });

for (const id of ["ball-small", "ball-small", "pillow"])
  await page.locator(`#environment-tiles [data-insert="${id}"]`).click();
await page.waitForTimeout(300);
report.steps.push({ step: "objects-inserted", ...(await envState()), screenshot: await shot("02c-environment-objects") });

await page.locator('#environment-tiles [data-tile="studio"]').click();
await page.waitForTimeout(300);
report.steps.push({ step: "environment-cascade", ...(await envState()) });
await page.locator('#environment-tiles [data-tile="bed"]').click();
await page.waitForTimeout(300);

// 3b · zoom until geometry fills the frame, so the panes float over busy 3D
await page.mouse.move(box.x + box.width * 0.5, box.y + box.height * 0.5);
for (let i = 0; i < 9; i++) { await page.mouse.wheel(0, -220); await page.waitForTimeout(80); }
await page.waitForTimeout(1500);
report.steps.push({
  step: "panes-over-busy-geometry",
  coverage: await coverage(),
  screenshot: await shot("04b-panes-over-geometry"),
});

// 4 · the orientation gizmo. It is a real 3D widget: it turns with the main
// camera, and its planes are picked by actual hit testing, which the widget
// reports through data-hover as the pointer crosses them.
async function probeGimbal() {
  const rect = await page.locator("#gimbal").boundingBox();
  const found = {};
  const N = 9;
  outer: for (let i = 1; i < N; i++)
    for (let j = 1; j < N; j++) {
      const x = rect.x + (rect.width * i) / N, y = rect.y + (rect.height * j) / N;
      await page.mouse.move(x, y);
      const hover = await page.locator("#gimbal").getAttribute("data-hover");
      if (hover && !found[hover]) found[hover] = { x, y };
      if (Object.keys(found).length === 3) break outer;
    }
  await page.mouse.move(rect.x - 20, rect.y - 20);
  return found;
}
// Orbit the main scene and record that the widget follows it.
const orbit = async (dx, dy) => {
  await page.mouse.move(box.x + box.width * 0.5, box.y + box.height * 0.5);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.5 + dx, box.y + box.height * 0.5 + dy, { steps: 12 });
  await page.mouse.up();
  await page.waitForTimeout(900);
};
const gimbalPixels = () => page.locator("#gimbal").screenshot().then((b) => b.length);
for (const [name, dx, dy] of [["front", 0, 0], ["three-quarter", 220, 60], ["profile", 220, 0]]) {
  await orbit(dx, dy);
  report.steps.push({
    step: `gimbal-orbit-${name}`,
    hoverable_planes: Object.keys(await probeGimbal()),
    gimbal_bytes: await gimbalPixels(),
    screenshot: await shot(`05-gimbal-orbit-${name}`),
  });
}
for (const plane of ["sagittal", "coronal", "transverse"]) {
  // The widget turns after every snap, so where a plane sits has to be found
  // again before each click rather than reused from an earlier orientation.
  const point = (await probeGimbal())[plane];
  if (!point) { report.steps.push({ step: `gimbal-${plane}`, error: "no visible hit area at this orbit" }); continue; }
  await page.mouse.move(point.x, point.y);
  await page.waitForTimeout(250);
  const hover = await page.locator("#gimbal").getAttribute("data-hover");
  const shotName = plane === "sagittal" ? await shot("05-gimbal-hover-sagittal") : null;
  await page.mouse.click(point.x, point.y);
  await page.waitForTimeout(1200);
  report.steps.push({
    step: `gimbal-${plane}`,
    hover_reported: hover,
    selected: await page.locator("#gimbal").getAttribute("data-selected"),
    coverage: await coverage(),
    hover_screenshot: shotName,
    screenshot: await shot(`05-gimbal-${plane}`),
  });
}

// 5 · transport speed menu
await page.locator("#speed").click();
await page.waitForTimeout(200);
const rates = await page.locator("#speed-menu button").allTextContents();
await page.locator("#speed-menu button").getByText("2x", { exact: true }).click();
report.steps.push({
  step: "transport",
  rates,
  speed: await page.locator("#speed").textContent(),
  play_enabled: await page.locator("#play").isEnabled(),
  screenshot: await shot("06-transport"),
});

// 6 · collapsed columns
await page.locator("#toggle-left").click();
await page.locator("#toggle-right").click();
await page.waitForTimeout(500);
report.steps.push({
  step: "collapsed",
  left_visible: await page.locator("#left-column").isVisible(),
  pane_visible: await page.locator("#pane-column").isVisible(),
  left_toggle_expanded: await page.locator("#toggle-left").getAttribute("aria-expanded"),
  right_toggle_expanded: await page.locator("#toggle-right").getAttribute("aria-expanded"),
  transport_visible: await page.locator("#transport").isVisible(),
  gimbal_visible: await page.locator("#gimbal").isVisible(),
  insets: await page.evaluate(() => {
    const style = getComputedStyle(document.getElementById("viewport"));
    return { left: style.getPropertyValue("--left-inset").trim(), right: style.getPropertyValue("--right-inset").trim() };
  }),
  screenshot: await shot("07-collapsed"),
});
await page.locator("#toggle-left").click();
await page.locator("#toggle-right").click();
await page.waitForTimeout(400);

// 6b · the monitor picker
await page.locator("#add-pane").click();
await page.waitForTimeout(300);
report.steps.push({
  step: "monitor-picker",
  available: await page.locator("#pane-picker button").allTextContents(),
  screenshot: await shot("07b-monitor-picker"),
});
await page.locator("#add-pane").click();
await page.waitForTimeout(200);

// 6c · panes collapsed to their titles
for (const id of ["playback", "live", "scene"]) await page.locator("#add-pane").click().catch(() => {});
await page.locator("#pane-picker").waitFor({ state: "visible" }).catch(() => {});
for (const title of ["Live body", "Body interaction", "Motor & skin inputs"]) {
  const entry = page.locator("#pane-picker button").getByText(title, { exact: true });
  if (await entry.count()) await entry.click();
}
await page.locator("#add-pane").click();
for (const id of ["selection", "live"]) await page.locator(`[data-pane='${id}'] .disclose`).click();
await page.waitForTimeout(400);
report.steps.push({
  step: "panes-collapsed",
  heights: await page.locator("#pane-column .pane:not([hidden])").evaluateAll((n) =>
    n.map((x) => ({ pane: x.dataset.pane, collapsed: x.dataset.collapsed, height: Math.round(x.getBoundingClientRect().height) }))),
  screenshot: await shot("07c-panes-collapsed"),
});

// 7 · a conforming tetrahedral domain: clothing and environment disappear
const domain = option("--domain", "conforming-domain-whole-body-0.01m");
await page.locator("#materialization").selectOption(domain);
await page.waitForFunction(() => document.getElementById("scene-status")?.textContent === "", null, { timeout: 240000 });
await page.waitForTimeout(3000);
await page.mouse.click(box.x + box.width * 0.5, box.y + box.height * 0.45);
await page.waitForTimeout(800);
report.steps.push({
  step: "conforming-domain",
  materialization: await page.locator("#materialization").inputValue(),
  clothing_visible: await page.locator("#clothing-section").isVisible(),
  environment_visible: await page.locator("#environment-section").isVisible(),
  note: await page.locator("[data-section='materialization'] .note").first().textContent(),
  roles: await page.locator("#domain-roles .check-row span").allTextContents(),
  picked: await page.locator("#domain-selected").textContent(),
  coverage: await coverage(),
  screenshot: await shot("08-conforming-domain"),
});

report.console = messages;
report.failed_requests = failed;
writeFileSync(resolve(out, "report.json"), JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
await browser.close();
