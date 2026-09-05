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
test("hemorrhage and saline protocol retains patient, native units, ordered stops and solver clock", () => {
  const result = scenarioInput("hemorrhage_saline", 2, "StandardFemale");
  assert.equal(result.patient, "StandardFemale");
  assert.deepEqual(
    result.interventions.map((x) => [x.kind, x.value]),
    [
      ["hemorrhage", 10],
      ["hemorrhage", 0],
      ["saline", 20],
      ["saline", 0],
    ],
  );
  assert.ok(
    result.interventions.every(
      (x, i, a) =>
        Math.abs(x.time_s * 50 - Math.round(x.time_s * 50)) < 1e-9 &&
        (!i || x.time_s > a[i - 1].time_s),
    ),
  );
  assert.throws(() => scenarioInput("unknown", 2));
  assert.throws(() => scenarioInput("baseline", 1.01));
  assert.throws(() => scenarioInput("baseline", 2, "../evil"));
});
test("source reproductive trajectory preserves days, channel units and prescribed-input provenance", async () => {
  const { trajectoryFromChannels } = await import("../src/state.js");
  const data = trajectoryFromChannels({
    time_days: [0, 10],
    time_s: [0, 864000],
    source_kind: "source_model_simulation",
    channels: [
      {
        id: "E2",
        unit: "microg_l",
        kind: "prescribed_time_input",
        values: [1, 2],
      },
    ],
    limitations: ["Inputs prescribed"],
  });
  assert.deepEqual(data.time_axis, [0, 10]);
  assert.equal(data.time_unit, "day");
  assert.equal(data.units.E2, "microg_l");
  assert.equal(data.channel_types.E2, "prescribed_time_input");
  assert.equal(data.metadata.source_kind, "source_model_simulation");
});

test("scalar cells preserve values and use a finite midpoint for constant fields", async () => {
  const { scalarCoordinates } = await import("../src/state.js");
  assert.deepEqual(scalarCoordinates([0, 1, 0], [10, 20], [10, 20]), [0, 1, 0]);
  assert.deepEqual(scalarCoordinates([0], [135], [135, 135]), [0.5]);
  assert.throws(() => scalarCoordinates([2], [10], [0, 20]));
});

test("environment overrides are optional, bounded and preserve explicit zero clothing", () => {
  assert.equal(
    scenarioInput("baseline", 2, "StandardMale", {
      ambient_temperature_c: "",
      clothing_clo: "",
    }).ambient_temperature_c,
    undefined,
  );
  assert.equal(
    scenarioInput("baseline", 2, "StandardMale", {
      ambient_temperature_c: "22",
      clothing_clo: "0",
    }).clothing_clo,
    0,
  );
  assert.throws(() =>
    scenarioInput("baseline", 2, "StandardMale", { ambient_temperature_c: 36 }),
  );
  assert.throws(() =>
    scenarioInput("baseline", 2, "StandardMale", { clothing_clo: -1 }),
  );
});

test("engine variant is explicit and restricted to verified source implementations", () => {
  assert.equal(
    scenarioInput("baseline", 2, "StandardFemale", {
      engine_variant: "saturation_bounds_heatflux",
    }).engine_variant,
    "saturation_bounds_heatflux",
  );
  assert.throws(() =>
    scenarioInput("baseline", 2, "StandardMale", { engine_variant: "unknown" }),
  );
});
