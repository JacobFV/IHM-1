import test from "node:test";
import assert from "node:assert/strict";
import { filterMonitors } from "../src/monitors.js";
test("monitor search combines words with overlapping category membership", () => {
  assert.deepEqual(
    filterMonitors("  finite LAPLACE ").map((m) => m.id),
    ["signals"],
  );
  assert.deepEqual(
    filterMonitors("contact", "mechanics").map((m) => m.id),
    ["contact"],
  );
  assert.deepEqual(
    filterMonitors("contact", "evidence").map((m) => m.id),
    ["contact"],
  );
  assert.deepEqual(filterMonitors("contact", "vascular"), []);
  assert.deepEqual(filterMonitors("invented medicine"), []);
});
