// Verify the structure-provenance API from inside a real browser against the running
// workbench: fetch the index and one record per origin family, and record every console
// message and failed request verbatim. Front-end independent on purpose - the inspector
// is being rewritten, so this asserts only the endpoint the new UI will call.
//
//   node scripts/verify_structure_provenance_api.mjs [--url http://127.0.0.1:8765] [--out DIR]
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
const out = resolve(option("--out", "test-results/structure-provenance"));
mkdirSync(out, { recursive: true });

const browser = await chromium.launch({
  executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
  args: ["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
});
const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
const console_ = [];
const failed = [];
page.on("console", (m) => console_.push(`[${m.type()}] ${m.text()} @ ${JSON.stringify(m.location())}`));
page.on("pageerror", (e) => console_.push(`[pageerror] ${e.message}`));
page.on("requestfailed", (r) => failed.push(`${r.url()} :: ${r.failure()?.errorText}`));
page.on("response", (r) => { if (r.status() >= 400) failed.push(`${r.url()} :: HTTP ${r.status()}`); });

const report = { url, structures: [] };
await page.goto(url, { waitUntil: "load" });
await page.waitForFunction(() => document.getElementById("model")?.value === "ihm-body", null, { timeout: 120000 })
  .catch(() => console_.push("[note] the page did not reach a loaded ihm-body model; endpoint checks continue"));

const get = (path) => page.evaluate(async (p) => {
  const response = await fetch(p);
  let body = null;
  try { body = await response.json(); } catch { body = null; }
  return { status: response.status, body };
}, path);

report.index = await get("/api/body/experiments/provenance");

// One case per origin family: canonical BodyParts3D, Z-Anatomy registered into the body,
// an unregistered display atlas surface, a published musculoskeletal model, a case-specific
// vascular surface, a published graph, a repository-built structure, a solver-generated one,
// and a canonical entity that carries no display structure.
const cases = (option("--structures",
  "body-bp3d-FJ1252,body-za-b0e9cdd77cc7f9fe,za-2b3e1cd7bf55f6d1,opensim-rajagopal-bone-pelvis_geom_1," +
  "vascular-aorta-wall,body-published-lymphatic-network,body-dynamic-scalp-hair,body-detail-microvascular," +
  "betse-tissue-cells,body-skin-dermis").split(","));

for (const id of cases) {
  const started = Date.now();
  const record = await get("/api/body/experiments/provenance-" + encodeURIComponent(id));
  const body = record.body || {};
  report.structures.push({
    id, status: record.status, seconds: (Date.now() - started) / 1000,
    tier: body.tier, dataset: body.dataset?.id, license: body.dataset?.license ?? null,
    source_file: body.source_file?.path ?? null, source_hash_verified: body.source_file?.sha256_verified ?? null,
    geometry_hash_verified: body.geometry?.sha256_verified ?? null,
    build_script: body.build?.script, build_commit: body.build?.commit,
    transforms: (body.transforms || []).map((t) => ({ kind: t.kind, residual_m: t.residual?.value ?? null })),
    required_present: body.completeness?.required_present, missing: body.completeness?.missing,
    derived_artifacts: (body.derived_artifacts || []).length,
  });
}

// An unknown id and a traversal attempt must fail closed rather than answer with a guess.
report.unknown_structure = await get("/api/body/experiments/provenance-not-a-structure");
report.traversal_rejected = await get("/api/body/experiments/provenance-" + encodeURIComponent("../../manifest"));

await page.screenshot({ path: resolve(out, "workbench.png"), fullPage: false });
report.console = console_;
report.failed_requests = failed;
writeFileSync(resolve(out, "report.json"), JSON.stringify(report, null, 1));
await browser.close();
console.log(JSON.stringify({
  out,
  index_status: report.index.status,
  answered: report.structures.filter((s) => s.status === 200).length,
  attempted: report.structures.length,
  unknown_structure_status: report.unknown_structure.status,
  traversal_status: report.traversal_rejected.status,
  console_messages: console_.length,
  failed_requests: failed.length,
}, null, 1));
