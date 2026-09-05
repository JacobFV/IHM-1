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
export function scenarioInput(scenario, seconds) {
  seconds = Number(seconds);
  if (!Number.isFinite(seconds) || seconds < 1 || seconds > 600)
    throw Error("Duration must be between 1 and 600 seconds");
  return {
    seconds,
    patient: "StandardMale",
    sample_hz: 10,
    interventions:
      scenario === "exercise"
        ? [
            { time_s: seconds * 0.2, kind: "exercise", value: 0.1 },
            { time_s: seconds * 0.7, kind: "exercise", value: 0 },
          ]
        : [],
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
