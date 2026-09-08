// Tileable value noise, shared by every procedural texture in the app.
//
// It exists because the first pass used one random number per texel. White
// noise at texel resolution has no mipmap: minified onto a floor that runs
// twenty metres to the wall it does not average to a tone, it crawls, and a
// camera that moves a few millimetres repaints half the screen. Noise built on
// a coarse lattice and interpolated has the same statistics at the scale it was
// authored for and filters correctly at every scale below it.
//
// Everything here is deterministic from its seed, so two runs, two machines and
// two screenshots agree.

function lattice(period, seed) {
  let state = seed >>> 0;
  const values = new Float32Array(period * period);
  for (let i = 0; i < values.length; i++) {
    state = (1664525 * state + 1013904223) >>> 0;
    values[i] = state / 4294967296;
  }
  return values;
}

// One octave: `period` cells across the unit tile, wrapped, so the result tiles
// seamlessly however far the texture repeats.
export function valueNoise(period, seed) {
  const values = lattice(period, seed);
  const at = (i, j) =>
    values[(((j % period) + period) % period) * period + (((i % period) + period) % period)];
  return (u, v) => {
    const x = u * period, y = v * period;
    const x0 = Math.floor(x), y0 = Math.floor(y);
    const fx = x - x0, fy = y - y0;
    // Smoothstep rather than a straight ramp: bilinear interpolation alone
    // leaves a visible crease along every lattice line.
    const sx = fx * fx * (3 - 2 * fx), sy = fy * fy * (3 - 2 * fy);
    const a = at(x0, y0) + (at(x0 + 1, y0) - at(x0, y0)) * sx;
    const b = at(x0, y0 + 1) + (at(x0 + 1, y0 + 1) - at(x0, y0 + 1)) * sx;
    return a + (b - a) * sy;
  };
}

// Octaves at halving weight. `periods` are given rather than derived so a
// caller can say how coarse the coarsest structure is in its own terms.
export function fbm(periods, seed) {
  const octaves = periods.map((period, i) => valueNoise(period, seed + 7919 * (i + 1)));
  const weights = periods.map((_, i) => 1 / 2 ** i);
  const total = weights.reduce((a, b) => a + b, 0);
  return (u, v) => octaves.reduce((sum, f, i) => sum + weights[i] * f(u, v), 0) / total;
}
