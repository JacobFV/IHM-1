import test from "node:test";
import assert from "node:assert/strict";
import {
  filterStructures,
  geometryArrays,
  chartPath,
  scenarioInput,
} from "../src/state.js";
test("search and systems never mix unregistered model families", () => {
  const rows = [
    { id: "1", name: "Left heart", system: "circulatory", model_id: "a" },
    { id: "2", name: "Heart", system: "circulatory", model_id: "b" },
    { id: "3", name: "Femur", system: "skeletal", model_id: "a" },
  ];
  assert.deepEqual(
    filterStructures(rows, "a", new Set(["circulatory"]), "heart").map(
      (x) => x.id,
    ),
    ["1"],
  );
  assert.equal(filterStructures(rows, "a", new Set(), "").length, 0);
});
test("geometry rejects corrupt coordinates and out-of-range triangles", () => {
  assert.throws(() => geometryArrays({ positions: [0, 0, NaN] }));
  assert.throws(() =>
    geometryArrays({ positions: [0, 0, 0], indices: [0, 1, 2] }),
  );
  assert.equal(
    geometryArrays({
      positions: [0, 0, 0, 1, 0, 0, 0, 1, 0],
      indices: [0, 1, 2],
    }).positions.length,
    9,
  );
});
test("chart supports constant values and skips missing points without joining gaps", () => {
  assert.equal(
    chartPath([0, 1, 2], [4, 4, 4], 100, 50),
    "M0,25 L50,25 L100,25",
  );
  assert.match(chartPath([0, 1, 2], [1, null, 3], 100, 50), /M0,50 M100,0/);
});
test("scenario validates bounded duration before transport", () => {
  assert.throws(() => scenarioInput("baseline", 0));
  assert.throws(() => scenarioInput("baseline", 3601));
  assert.deepEqual(scenarioInput("baseline", 60), {
    seconds: 60,
    patient: "StandardMale",
    sample_hz: 10,
    interventions: [],
  });
});
test("Laplace magnitude uses the selected damping row and preserves physical units", async () => {
  const { spectralSeries } = await import("../src/state.js");
  const run = {
    variables: [
      { id: "p", frequency_hz: [0, 1], psd: [2, 3], psd_unit: "Pa²/Hz" },
    ],
    laplace: {
      frequency_hz: [0, 1],
      sigma_per_s: [0, 1],
      variables: [
        {
          id: "p",
          unit: "Pa*s",
          real: [
            [3, 0],
            [0, 5],
          ],
          imag: [
            [4, 0],
            [2, 12],
          ],
        },
      ],
    },
  };
  assert.deepEqual(spectralSeries(run, "p", "laplace", 1), {
    x: [0, 1],
    y: [2, 13],
    unit: "Pa*s",
  });
  assert.deepEqual(spectralSeries(run, "p", "psd", 0).y, [2, 3]);
});
