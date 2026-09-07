import test from "node:test";
import assert from "node:assert/strict";
import { filterStructures, geometryArrays } from "../src/state.js";

test("search and systems never mix unregistered model families", () => {
  const rows = [
    { id: "1", name: "Left heart", system: "circulatory", model_id: "a" },
    { id: "2", name: "Heart", system: "circulatory", model_id: "b" },
    { id: "3", name: "Femur", system: "skeletal", model_id: "a" },
  ];
  assert.deepEqual(
    filterStructures(rows, "a", new Set(["circulatory"]), "heart").map((x) => x.id),
    ["1"],
  );
  assert.equal(filterStructures(rows, "a", new Set(), "").length, 0);
});

test("geometry rejects corrupt coordinates and out-of-range triangles", () => {
  assert.throws(() => geometryArrays({ positions: [0, 0, NaN] }));
  assert.throws(() => geometryArrays({ positions: [0, 0, 0], indices: [0, 1, 2] }));
  assert.equal(
    geometryArrays({ positions: [0, 0, 0, 1, 0, 0, 0, 1, 0], indices: [0, 1, 2] }).positions.length,
    9,
  );
});
