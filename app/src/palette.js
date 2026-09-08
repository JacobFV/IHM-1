// Tissue colour palettes. The server publishes an index of named palettes and,
// alongside it and on its own axis, a set of skin-tone options; both are
// read-only here. Nothing in this file decides what a tissue should look like.
// It reads what a palette says, converts a declared skin tone into the hex the
// renderer needs, and reports every structure the palette does not cover rather
// than inventing a colour for it.

const STORE = "ihm.tissue-palette.v1";
// The server's own default, restated so a first paint has something to use
// before the index arrives. The index is still the authority: whatever it names
// as default_palette wins over this the moment it is read.
export const DEFAULT_PALETTE = "didactic";
export const DEFAULT_SKIN_TONE = "not_applied";

export function readPaletteChoice() {
  const fallback = { palette: DEFAULT_PALETTE, skin_tone: DEFAULT_SKIN_TONE };
  try {
    const saved = JSON.parse(localStorage.getItem(STORE));
    if (!saved || typeof saved !== "object") return fallback;
    return {
      palette: typeof saved.palette === "string" ? saved.palette : fallback.palette,
      skin_tone: typeof saved.skin_tone === "string" ? saved.skin_tone : fallback.skin_tone,
    };
  } catch {
    return fallback;
  }
}
export function writePaletteChoice(choice) {
  try { localStorage.setItem(STORE, JSON.stringify(choice)); } catch {}
}

// CIE L*a*b* to sRGB hex. Every skin-tone option in the index declares D65 with
// the two-degree observer, which is the white point sRGB itself is defined on,
// so no chromatic adaptation is involved and none is done here; an option that
// declared any other illuminant is refused rather than converted wrongly.
//
// The constants are lifted from scripts/build_tissue_colour_palettes.py verbatim
// -- the exact CIE kappa and epsilon, the same IEC 61966-2-1 matrix, the same
// rounding -- because the substitution below depends on this reproducing the
// build's own hex to the byte. It does: the neutral option, L*50 a*10 b*17,
// comes back as #91705b, the hex the shipped realistic palette carries for its
// skin roles, and every other option reproduces the hex the build's self-test
// prints for it.
const D65_2 = [95.047, 100.0, 108.883];
const XYZ_TO_RGB = [
  [3.2404542, -1.5371385, -0.4985314],
  [-0.969266, 1.8760108, 0.041556],
  [0.0556434, -0.2040259, 1.0572252],
];
export function labToHex(lab, illuminantObserver = "D65_2") {
  if (illuminantObserver !== "D65_2") return null;
  if (!Array.isArray(lab) || lab.length < 3 || !lab.slice(0, 3).every(Number.isFinite)) return null;
  const [L, a, b] = lab;
  const fy = (L + 16) / 116, fx = fy + a / 500, fz = fy - b / 200;
  const expand = (t) => (t ** 3 > 216 / 24389 ? t ** 3 : (116 * t - 16) / (24389 / 27));
  const xyz = [expand(fx) * D65_2[0], expand(fy) * D65_2[1], expand(fz) * D65_2[2]].map((v) => v / 100);
  const channel = (row) => {
    const linear = Math.min(1, Math.max(0, row[0] * xyz[0] + row[1] * xyz[1] + row[2] * xyz[2]));
    const encoded = linear <= 0.0031308 ? 12.92 * linear : 1.055 * linear ** (1 / 2.4) - 0.055;
    return Math.min(255, Math.max(0, Math.round(encoded * 255))).toString(16).padStart(2, "0");
  };
  return "#" + XYZ_TO_RGB.map(channel).join("");
}

const normalise = (hex) => (typeof hex === "string" && /^#[0-9a-fA-F]{6}$/.test(hex) ? hex.toLowerCase() : null);

// What to paint every structure, given one palette and one skin tone.
//
// A structure the palette does not name keeps the colour it already had. That
// is reported rather than hidden: a palette that silently covers less than the
// scene would otherwise read as a palette that covers all of it.
//
// The skin tone is substituted by the colour the palette itself declares it was
// built with, not by guessing which structures are skin. A palette record names
// its skin_tone and the L*a*b* it was built from; every structure painted in
// exactly that colour is a structure the tone governs, and every other colour --
// palmar and plantar skin among them, which is measured on its own -- is left
// alone. A palette that declares no skin_tone has no such axis and is returned
// untouched.
export function resolvePalette({ palette, structures, base, toneHex }) {
  const colours = new Map();
  const missing = [];
  const from = palette?.colours || {};
  const declared = palette?.skin_tone
    ? normalise(labToHex(palette.skin_tone.colour?.lab, palette.skin_tone.colour?.illuminant_observer))
    : null;
  const tone = normalise(toneHex);
  const substitute = declared && tone && tone !== declared ? tone : null;
  let skinCount = 0;
  for (const s of structures || []) {
    const painted = normalise(from[s.id]);
    if (!painted) {
      missing.push(s.id);
      const kept = base?.get(s.id);
      if (kept) colours.set(s.id, kept);
      continue;
    }
    if (declared && painted === declared) {
      skinCount++;
      colours.set(s.id, substitute || painted);
      continue;
    }
    colours.set(s.id, painted);
  }
  return { colours, missing, skinCount, skinToneApplies: !!declared };
}
