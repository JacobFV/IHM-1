// Measure how far an equal-length drag turns the camera, from several starting
// orientations, in degrees per pixel. "Feels stuck" is not a measurement; this
// reports the numbers behind it.
//
//   node app/test/orbit-isotropy.mjs [--url http://127.0.0.1:8765] [--out DIR] [--tag before]
//
// For each start orientation it drags horizontally and vertically by the same
// pixel count and records the change in camera *direction* (where the camera
// sits relative to the target) and in full camera *orientation* (which includes
// roll). It then asks the question the owner actually asked: can a sustained
// vertical drag carry the camera over the head, around the coronal plane, and
// back? A turntable cannot, and the sweep test reports by how much it falls
// short of the rotation the drag asked for.
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
const tag = option("--tag", "measure");
const out = resolve(option("--out", `test-results/orbit-isotropy/${tag}`));
const DRAG = Number(option("--drag", "120"));
mkdirSync(out, { recursive: true });

const deg = (r) => (r * 180) / Math.PI;
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.hypot(...a);
const unit = (a) => { const n = norm(a) || 1; return a.map((x) => x / n); };
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const angleBetween = (a, b) => deg(Math.acos(Math.min(1, Math.max(-1, dot(unit(a), unit(b))))));
// Relative rotation q1 * inverse(q0), reported as angle plus world axis.
function relativeRotation(q0, q1) {
  const inv = [-q0[0], -q0[1], -q0[2], q0[3]];
  const [x1, y1, z1, w1] = q1, [x2, y2, z2, w2] = inv;
  const q = [
    w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
    w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
    w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
  ];
  const s = Math.hypot(q[0], q[1], q[2]);
  const angle = deg(2 * Math.atan2(s, Math.abs(q[3])));
  const sign = q[3] < 0 ? -1 : 1;
  const axis = s < 1e-9 ? [0, 0, 0] : [q[0], q[1], q[2]].map((v) => (sign * v) / s);
  return { angle, axis };
}

const browser = await chromium.launch({
  executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
  args: ["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const console_messages = [], network_errors = [];
page.on("console", (m) => console_messages.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => console_messages.push(`[pageerror] ${e.message}\n${e.stack || ""}`));
page.on("requestfailed", (r) => network_errors.push(`${r.url()} :: ${r.failure()?.errorText}`));
page.on("response", (r) => { if (r.status() >= 400) network_errors.push(`${r.url()} :: HTTP ${r.status()}`); });

await page.goto(url, { waitUntil: "load" });
await page.waitForSelector("#materialization", { timeout: 120000 });
await page.waitForFunction(() => document.getElementById("scene-status")?.textContent === "", null, { timeout: 240000 });
await page.waitForFunction(() => !!globalThis.__ihmCamera, null, { timeout: 30000 });
await page.waitForTimeout(2000);
process.stderr.write("loaded\n");

const box = await page.locator("#scene").boundingBox();
const centre = { x: box.x + box.width * 0.5, y: box.y + box.height * 0.5 };
// A drag of one canvas height is one full turn, which is what the requested
// rotation below is measured against.
const canvasHeight = await page.locator("#scene").evaluate((node) => node.clientHeight);
const pose = () => page.evaluate(() => globalThis.__ihmCamera.pose());
const snap = (plane) => page.evaluate((p) => globalThis.__ihmCamera.snap(p), plane);
const setPose = (p) => page.evaluate((x) => globalThis.__ihmCamera.setPose(x), p);

// Damping spreads one drag over many frames, and under swiftshader the render
// loop runs at a few frames a second, so waiting in wall-clock time measures the
// software rasteriser rather than the orbit. `settle` runs the very update the
// render loop runs, as many times as damping needs, and returns a camera that
// has finished moving.
const settle = () => page.evaluate(() => globalThis.__ihmCamera.settle(900));
async function drag(dx, dy, steps = 8) {
  await page.mouse.move(centre.x, centre.y);
  await page.mouse.down();
  await page.mouse.move(centre.x + dx, centre.y + dy, { steps });
  await page.mouse.up();
  return settle();
}

// Four starts: the three anatomical snaps, plus a general oblique that no snap
// produces, plus a pose deliberately close to the superior pole.
const starts = [];
{
  const coronal = await snap("coronal");
  await settle();
  const t = coronal.target, p = coronal.position;
  const radius = norm([p[0] - t[0], p[1] - t[1], p[2] - t[2]]);
  const at = (d, up) => ({ target: t, up, position: [t[0] + d[0] * radius, t[1] + d[1] * radius, t[2] + d[2] * radius] });
  starts.push(
    { name: "coronal-snap", kind: "snap", plane: "coronal" },
    { name: "sagittal-snap", kind: "snap", plane: "sagittal" },
    { name: "transverse-snap", kind: "snap", plane: "transverse" },
    { name: "oblique", kind: "pose", pose: at(unit([0.6, 0.5, 0.62]), [0, 1, 0]) },
    // 12 degrees off the superior pole: the arc over the head, where a
    // turntable's azimuth degenerates.
    { name: "near-superior-pole", kind: "pose", pose: at([0, Math.cos(0.209), Math.sin(0.209)], [0, 1, 0]) },
  );
}

const started = Date.now();
const mark = (what) => process.stderr.write(`${((Date.now() - started) / 1000).toFixed(1)}s ${what}\n`);
const report = { url, tag, drag_px: DRAG, canvas_height_px: canvasHeight, expected_deg_per_px: +(360 / canvasHeight).toFixed(5), starts: [] };
for (const start of starts) {
  const place = async () => {
    if (start.kind === "snap") await snap(start.plane); else await setPose(start.pose);
    await settle();
    return pose();
  };
  const base = await place();
  const horizontal = await drag(DRAG, 0);
  const afterH = { direction: angleBetween(base.direction, horizontal.direction), ...relativeRotation(base.quaternion, horizontal.quaternion) };
  await place();
  const vertical = await drag(0, DRAG);
  const afterV = { direction: angleBetween(base.direction, vertical.direction), ...relativeRotation(base.quaternion, vertical.quaternion) };
  const per = (v) => +(v / DRAG).toFixed(5);
  const ratio = Math.max(afterH.direction, afterV.direction) / Math.max(1e-6, Math.min(afterH.direction, afterV.direction));
  report.starts.push({
    name: start.name,
    start_direction: base.direction.map((v) => +v.toFixed(4)),
    start_up: base.up.map((v) => +v.toFixed(4)),
    horizontal: { deg: +afterH.direction.toFixed(3), deg_per_px: per(afterH.direction), orientation_deg: +afterH.angle.toFixed(3), axis: afterH.axis.map((v) => +v.toFixed(3)) },
    vertical: { deg: +afterV.direction.toFixed(3), deg_per_px: per(afterV.direction), orientation_deg: +afterV.angle.toFixed(3), axis: afterV.axis.map((v) => +v.toFixed(3)) },
    asymmetry: +ratio.toFixed(3),
  });
  mark(`start ${start.name}`);
  // Written as it goes: a run that is cut short still reports what it measured.
  writeFileSync(resolve(out, "orbit-isotropy.json"), JSON.stringify(report, null, 2));
}

// The owner's complaint, as a measurement: hold a vertical drag and try to
// travel the coronal plane, over the head and down the back. The drag asks for
// 360 * total_px / clientHeight degrees; record how much arrives.
async function sweep(name, place, dx, dy, pulls = 6) {
  await place();
  await settle();
  const first = await pose();
  const path = [];
  let previous = first, travelled = 0;
  for (let i = 0; i < pulls; i++) {
    const now = await drag(dx, dy);
    travelled += angleBetween(previous.direction, now.direction);
    previous = now;
    path.push({ pull: i + 1, direction: now.direction.map((v) => +v.toFixed(4)), from_start_deg: +angleBetween(first.direction, now.direction).toFixed(2), travelled_deg: +travelled.toFixed(2) });
  }
  const requested = (360 * Math.hypot(dx, dy) * pulls) / canvasHeight;
  return { name, requested_deg: +requested.toFixed(1), travelled_deg: +travelled.toFixed(2), delivered_fraction: +(travelled / requested).toFixed(3), path };
}
report.sweeps = [];
for (const [name, dx, dy] of [
  ["coronal-plane · vertical drag over the head", 0, -160],
  ["transverse-plane · horizontal drag around the body", 160, 0],
]) {
  report.sweeps.push(await sweep(name, () => snap("coronal"), dx, dy));
  mark(`sweep ${name}`);
  report.console_messages = console_messages;
  report.network_errors = network_errors;
  writeFileSync(resolve(out, "orbit-isotropy.json"), JSON.stringify(report, null, 2));
}
await browser.close();
console.log(JSON.stringify(report, null, 2));
