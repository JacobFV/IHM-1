export function filterStructures(rows, model, systems, query = "") {
  const q = query.trim().toLowerCase();
  return rows.filter(
    (x) =>
      x.model_id === model &&
      systems.has(x.system) &&
      `${x.name} ${x.id}`.toLowerCase().includes(q),
  );
}
export function geometryArrays(g) {
  const p = g.positions?.flat?.() ?? g.positions;
  if (!p?.length || p.length % 3 || !Array.from(p).every(Number.isFinite))
    throw Error("Invalid geometry coordinates");
  const i = g.indices?.flat?.() ?? g.indices;
  if (
    i &&
    (!Array.from(i).every(
      (x) => Number.isInteger(x) && x >= 0 && x < p.length / 3,
    ) ||
      i.length % 3)
  )
    throw Error("Invalid geometry indices");
  return {
    positions: new Float32Array(p),
    indices: i ? new Uint32Array(i) : null,
  };
}
export function chartPath(times, values, w = 600, h = 120) {
  const valid = values.filter((x) => x !== null && Number.isFinite(x));
  if (!valid.length) return "";
  const low = Math.min(...valid),
    high = Math.max(...valid),
    t0 = times[0],
    span = times.at(-1) - t0 || 1;
  let pen = false;
  return values
    .map((v, i) => {
      if (v === null || !Number.isFinite(v)) {
        pen = false;
        return "";
      }
      const x = ((times[i] - t0) / span) * w,
        y = high === low ? h / 2 : h - ((v - low) / (high - low)) * h;
      const p = `${pen ? "L" : "M"}${+x.toFixed(2)},${+y.toFixed(2)}`;
      pen = true;
      return p;
    })
    .filter(Boolean)
    .join(" ");
}
export function scenarioInput(
  scenario,
  seconds,
  patient = "StandardMale",
  environment = {},
) {
  seconds = Number(seconds);
  if (!Number.isFinite(seconds) || seconds < 1 || seconds > 600)
    throw Error("Duration must be between 1 and 600 seconds");
  if (Math.abs(seconds * 50 - Math.round(seconds * 50)) > 1e-7)
    throw Error("Duration must align to the 0.02 s native solver step");
  if (!/^[A-Za-z0-9_]+$/.test(patient)) throw Error("Invalid upstream patient");
  if (!["baseline", "exercise", "hemorrhage_saline"].includes(scenario))
    throw Error("Unknown protocol");
  const at = (f) => Math.round(seconds * f * 50) / 50;
  const interventions =
    scenario === "exercise"
      ? [
          { time_s: at(0.2), kind: "exercise", value: 0.1 },
          { time_s: at(0.7), kind: "exercise", value: 0 },
        ]
      : scenario === "hemorrhage_saline"
        ? [
            { time_s: at(0.2), kind: "hemorrhage", value: 10 },
            { time_s: at(0.4), kind: "hemorrhage", value: 0 },
            { time_s: at(0.45), kind: "saline", value: 20 },
            { time_s: at(0.7), kind: "saline", value: 0 },
          ]
        : [];
  const overrides = {};
  for (const [key, low, high] of [
    ["ambient_temperature_c", 10, 35],
    ["clothing_clo", 0, 3],
  ]) {
    const raw = environment[key];
    if (raw === undefined || raw === null || raw === "") continue;
    const value = Number(raw);
    if (!Number.isFinite(value) || value < low || value > high)
      throw Error(`${key} must be between ${low} and ${high}`);
    overrides[key] = value;
  }
  return { seconds, patient, sample_hz: 10, interventions, ...overrides };
}
export function trajectoryFromChannels(data) {
  return {
    time_s: data.time_s,
    time_axis: data.time_days || data.time_s,
    time_unit: data.time_days ? "day" : "s",
    values: Object.fromEntries(data.channels.map((c) => [c.id, c.values])),
    units: Object.fromEntries(data.channels.map((c) => [c.id, c.unit])),
    channel_types: Object.fromEntries(data.channels.map((c) => [c.id, c.kind])),
    metadata: {
      source_kind: data.source_kind,
      source: data.source,
      limitations: data.limitations || [],
    },
  };
}
export function spectralSeries(run, id, mode = "psd", sigma = 0) {
  if (mode === "laplace") {
    const variable = run?.laplace?.variables.find((v) => v.id === id);
    return {
      x: run?.laplace?.frequency_hz || [],
      y: (variable?.real[sigma] || []).map((r, i) =>
        Math.hypot(r, variable.imag[sigma][i]),
      ),
      unit: variable?.unit || "",
    };
  }
  const v = run?.variables.find((v) => v.id === id);
  return { x: v?.frequency_hz || [], y: v?.psd || [], unit: v?.psd_unit || "" };
}

export function scalarCoordinates(cellIds, values, range) {
  const [low, high] = range;
  if (![low, high].every(Number.isFinite) || high < low)
    throw Error("Invalid scalar range");
  return cellIds.map((id) => {
    if (!Number.isInteger(id) || id < 0 || !Number.isFinite(values[id]))
      throw Error("Invalid scalar cell");
    return high === low
      ? 0.5
      : Math.max(0, Math.min(1, (values[id] - low) / (high - low)));
  });
}
