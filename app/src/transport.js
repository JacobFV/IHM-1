// The transport, and the one distinction it exists to make.
//
// The simulation is the app. It is always live, and the transport is a single
// play/pause that runs or holds *it*. A recording is a different noun: you open
// one, and while it is open the transport scrubs the recording and the header
// says which. That is the difference between a video player and a camera, and
// it is carried by what the control is attached to rather than by a tooltip.
//
// Before this there were two buttons that looked alike -- a play that replayed a
// stored trajectory and said so only in its title text ("This does not run the
// simulation"), and a Start body that posted to /api/embodied/sessions. Nothing
// about their appearance said they were different nouns.
//
// This module is the state, with no DOM in it, so what the control means is a
// value a test can read.

/** Which noun the transport is attached to right now. */
export const SIMULATION = "simulation";
export const RECORDING = "recording";

/** What the one transport control is, given everything that bears on it.
 *
 * `open` is the whole switch: a recording is open or it is not. There is no
 * mode toggle, because opening a recording is not a mode -- it is a document.
 */
export function transportState({
  open = false, recording = null, frames = 0, index = 0, playing = false,
  recordingError = "", started = false, running = false, busy = false,
  staticVolume = false, sessionError = "",
} = {}) {
  if (open) {
    const unavailable = recordingError
      ? recordingError
      : frames < 2 ? "This recording holds fewer than two frames; there is nothing to scrub."
      : "";
    return {
      owner: RECORDING,
      label: "Recording",
      name: recording?.label || "Recording",
      symbol: playing ? "❚❚" : "▶",
      action: playing ? "hold" : "play",
      disabled: !!unavailable,
      title: unavailable || (playing
        ? "Hold the recording."
        : "Play the recording. The simulation is not running while a recording is open."),
      scrub: true,
      frames,
      index: Math.min(Math.max(0, index), Math.max(0, frames - 1)),
      speed: true,
      note: unavailable || "Scrubbing a stored trajectory. Close the recording to return to the live body.",
    };
  }
  const unavailable = staticVolume
    ? "A conforming tetrahedral volume is static; there is nothing to run."
    : busy ? "The body is initializing or closing; the transport is held until it settles."
    : "";
  return {
    owner: SIMULATION,
    label: "Simulation",
    name: "Live body",
    symbol: running ? "❚❚" : "▶",
    // Start and resume are the same request to the reader: run it.
    action: running ? "hold" : started ? "resume" : "start",
    disabled: !!unavailable,
    title: unavailable || (running
      ? "Hold the simulation."
      : started ? "Resume the simulation." : "Run the simulation."),
    scrub: false,
    frames: 0,
    index: 0,
    speed: false,
    note: unavailable || sessionError ||
      (running ? "The body is stepping." : started
        ? "Held at the last accepted body state. Inputs apply on resume."
        : "The simulation is the app. Press play to run the body."),
  };
}

/** The recordings that can be opened, from what the server actually has.
 *
 * The canonical trajectory is one entry; every completed canonical scenario run
 * is another. A run that failed, or that ran without the canonical body, is not
 * a body recording and is not offered.
 */
export function recordingCatalog({ canonical = null, scenarios = null } = {}) {
  const out = [];
  if (canonical?.available) {
    out.push({
      id: "canonical",
      label: "Canonical trajectory",
      detail: canonical.detail || "The materialized native run for this body.",
      url: "/api/body/trajectory?view=display",
    });
  }
  for (const run of scenarios?.runs || []) {
    if (run.status !== "completed" || !run.canonical_body) continue;
    const frames = run.summary?.canonical_body?.frames;
    out.push({
      id: run.id,
      label: `Scenario ${run.id}`,
      detail: Number.isFinite(frames) ? `${frames} computed frames` : "Completed canonical run",
      url: `/api/body/trajectory?view=display&run=${encodeURIComponent(run.id)}`,
    });
  }
  return out;
}

/** What the header says while a recording is open. Never shown otherwise.
 *
 * A recording that would not load still gets its header -- you opened it, so it
 * is the thing on screen -- and the header says what went wrong instead of
 * reciting a frame count the file does not actually deliver.
 */
export function recordingHeader(recording, { frames = 0, index = 0, time_s = null, error = "" } = {}) {
  if (!recording) return null;
  if (error) return { name: recording.label, detail: error, failed: true };
  const at = Number.isFinite(time_s) ? ` · t ${time_s.toFixed(3)} s` : "";
  return {
    name: recording.label,
    detail: `${recording.detail} · frame ${frames ? index + 1 : 0} / ${frames}${at}`,
    failed: false,
  };
}
