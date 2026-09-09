// The prompt bar: attach an image or an audio file, or type text, and read back
// what the brain retrieved.
//
// **Retrieved.** Every output on this path is a nearest entry in the corpus, not
// a generation: the brain does not draw an image, it finds the one its cortical
// state is closest to. Showing a single returned image without that label would
// be the most misleading thing this app could do, so the panel shows what
// scripts/introspect.py shows -- the top five with their cosine similarities and
// the chance level -- and for the same reason it shows them there.
//
// The top-1 similarity is not a confidence. It falls only 0.706 -> 0.686 under a
// scramble that costs 63 points of accuracy, so a wrong answer looks exactly as
// certain as a right one. The caveat is rendered next to the number, not filed
// behind a disclosure.

const MODALITIES = ["image", "audio", "text"];
const ACCEPT = "image/png,image/jpeg,image/webp,image/gif,audio/wav,audio/mpeg,audio/ogg,audio/flac";
const MAX_BYTES = 8 * 1024 * 1024;

/** The caveats the panel always carries, whatever the server sends back. */
export const STANDING_CAVEATS = [
  "The top-1 cosine similarity is NOT a confidence: under a scramble that costs 63 points of accuracy it falls only 0.706 → 0.686. A wrong answer looks as certain as a right one; the margin over rank 2 is the quantity that separates them.",
  "Every result is a nearest entry retrieved from the corpus. Nothing here is generated.",
];

/** Normalize what /api/brain/prompt returns into what the panel draws.
 *
 * The served envelope is {path, output_kind, output_is_retrieval, retrieval:
 * {chance, n_bank_entries}, results[{rank, cosine_similarity, corpus_entry,
 * text, image_png_base64?, audio_wav_base64?}], top1_cosine,
 * margin_over_rank2, similarity_is_not_confidence, caveats[]}. The flatter
 * {modality, results[{similarity}], chance} shape the design doc sketches is
 * read too, so the panel does not break if the endpoint moves back to it.
 *
 * A missing chance level is reported missing rather than defaulted: a number
 * that was never measured must not appear as though it was.
 */
export function readout(response) {
  if (!response || typeof response !== "object")
    return { error: "The prompt endpoint returned nothing readable." };
  if (response.error) return {
    error: String(response.error),
    why: response.why || null,
    declared_paths: Array.isArray(response.declared_paths) ? response.declared_paths : null,
    candidates: Array.isArray(response.candidates) ? response.candidates
      : Array.isArray(response.example_concepts) ? response.example_concepts : null,
  };
  const bank = response.retrieval && typeof response.retrieval === "object" ? response.retrieval : {};
  const raw = Array.isArray(response.results) ? response.results : [];
  const results = raw.slice(0, 5).map((entry, i) => {
    const similarity = Number(entry?.cosine_similarity ?? entry?.similarity ?? entry?.cosine ?? NaN);
    const corpus = entry?.corpus_entry || {};
    const image = entry?.image_png_base64 ? `data:image/png;base64,${entry.image_png_base64}` : entry?.image_url;
    const audio = entry?.audio_wav_base64 ? `data:audio/wav;base64,${entry.audio_wav_base64}` : entry?.audio_url;
    return {
      rank: Number.isFinite(entry?.rank) ? entry.rank : i + 1,
      similarity: Number.isFinite(similarity) ? similarity : null,
      label: String(corpus.concept ?? entry?.text ?? entry?.label ?? entry?.entry_id ?? `entry ${i + 1}`),
      entry_id: corpus.bank_index ?? entry?.entry_id ?? entry?.id ?? null,
      image_url: image || null,
      audio_url: audio || null,
      text: entry?.text ?? null,
      // What the bank actually holds. On this path it is measured neural
      // responses, never the images, and saying so is the whole point.
      bank_is: corpus.bank_is ?? null,
    };
  });
  const chance = Number(bank.chance ?? response.chance);
  const declaredMargin = Number(response.margin_over_rank2);
  const margin = Number.isFinite(declaredMargin) ? declaredMargin
    : results.length > 1 && results[0].similarity != null && results[1].similarity != null
      ? results[0].similarity - results[1].similarity : null;
  const notConfidence = response.similarity_is_not_confidence;
  const caveats = [...new Set([
    ...(notConfidence?.measured ? [`${notConfidence.headline || "Similarity is not confidence"}: ${notConfidence.measured} ${notConfidence.consequence || ""}`.trim()] : []),
    ...(notConfidence?.use_instead ? [String(notConfidence.use_instead)] : []),
    ...(response.output_is ? [String(response.output_is)] : []),
    ...(Array.isArray(response.caveats) ? response.caveats.map(String) : []),
    ...STANDING_CAVEATS,
  ])];
  const [inModality, outModality] = String(response.path || "").split("->");
  return {
    // The server says so on this build. It stays true whether or not it does:
    // this path is a nearest-entry lookup, and the label is not the server's to
    // withdraw.
    retrieval: true,
    declared_retrieval: response.output_is_retrieval !== false && response.retrieval !== false
      && response.generated !== true,
    modality: MODALITIES.includes(inModality) ? inModality
      : MODALITIES.includes(response.modality) ? response.modality : null,
    output_modality: MODALITIES.includes(response.output_kind) ? response.output_kind
      : MODALITIES.includes(outModality) ? outModality
      : MODALITIES.includes(response.output_modality) ? response.output_modality : null,
    corpus: results[0]?.bank_is ?? bank.ranked_by ?? response.corpus ?? response.retrieval_bank ?? null,
    bank_entries: Number.isFinite(bank.n_bank_entries) ? bank.n_bank_entries : null,
    how: response.how ? String(response.how) : null,
    results,
    chance: Number.isFinite(chance) ? chance : null,
    margin,
    caveats,
    empty: !results.length,
  };
}

/** What a similarity bar is worth, as a fraction of the widest one drawn.
 *
 * Cosine runs -1..1 and the interesting band is narrow, so the bars are scaled
 * to the top result rather than to 1.0 -- and the axis says so, because a bar
 * scaled to its own maximum otherwise reads as a probability.
 */
export function barFractions(results) {
  const values = results.map((r) => r.similarity).filter((v) => v != null);
  if (!values.length) return results.map(() => 0);
  const high = Math.max(...values, 1e-6);
  return results.map((r) => (r.similarity == null ? 0 : Math.max(0, r.similarity / high)));
}

const el = (tag, cls, text) => {
  const node = document.createElement(tag);
  if (cls) node.className = cls;
  if (text !== undefined) node.textContent = text;
  return node;
};

async function base64(file) {
  const buffer = await file.arrayBuffer();
  if (buffer.byteLength > MAX_BYTES) throw Error("Attachment is larger than 8 MB.");
  let binary = "";
  const bytes = new Uint8Array(buffer);
  for (let i = 0; i < bytes.length; i += 0x8000)
    binary += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(binary);
}

//: The six declared paths, and only those. image->audio and audio->image are
//: refused by the endpoint because the two corpora share no entry, and text->text
//: because a lookup returning its own word is not a decode. Offering them here
//: and letting the server refuse would be offering something that does not exist.
export const PATHS = {
  image: ["image", "text"],
  audio: ["audio", "text"],
  text: ["image", "audio"],
};
export function wantsFor(modality) { return PATHS[modality] || []; }

/** Mount the bar. `post(body)` does the request; injected so a test can drive it. */
export function mountPromptBar(host, { post } = {}) {
  const panel = el("div", "prompt-panel");
  panel.id = "prompt-panel";
  panel.hidden = true;
  const bar = el("form", "prompt-bar");
  bar.id = "prompt-bar";

  const attach = el("label", "prompt-attach");
  attach.title = "Attach an image or an audio file";
  const file = document.createElement("input");
  file.type = "file";
  file.accept = ACCEPT;
  file.id = "prompt-file";
  file.setAttribute("aria-label", "Attach an image or audio file");
  attach.append(file, el("span", null, "Attach"));

  const text = document.createElement("input");
  text.type = "text";
  text.id = "prompt-text";
  text.placeholder = "Type a concept, or attach an image or a sound…";
  text.setAttribute("aria-label", "Prompt the brain");

  const want = document.createElement("select");
  want.id = "prompt-want";
  want.className = "prompt-want";
  want.setAttribute("aria-label", "What to retrieve");

  const send = el("button", "prompt-send", "Send");
  send.type = "submit";
  send.id = "prompt-send";

  const chip = el("span", "prompt-chip");
  chip.id = "prompt-chip";
  chip.hidden = true;
  bar.append(attach, chip, text, el("span", "prompt-arrow", "→"), want, send);
  host.append(panel, bar);

  let attached = null;
  const inputModality = () =>
    attached ? (attached.type.startsWith("audio") ? "audio" : "image") : "text";
  function syncWants() {
    const options = wantsFor(inputModality());
    const held = want.value;
    want.replaceChildren();
    for (const value of options) {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      want.append(option);
    }
    want.value = options.includes(held) ? held : options[0] || "";
  }
  file.onchange = () => {
    attached = file.files?.[0] || null;
    chip.hidden = !attached;
    chip.textContent = attached ? `${attached.name} · ${attached.type || "unknown type"}` : "";
    syncWants();
  };
  chip.onclick = () => { attached = null; file.value = ""; chip.hidden = true; syncWants(); };
  syncWants();

  function status(kind, message, detail) {
    panel.hidden = false;
    panel.dataset.state = kind;
    panel.replaceChildren(el("p", "prompt-status", message));
    if (detail) panel.append(el("p", "note", detail));
  }

  // A refusal is a result too: the endpoint says which paths exist and which
  // words the corpus holds, and hiding that behind "request failed" would make
  // a bounded corpus look like a broken server.
  function failure(view) {
    status("error", view.error);
    if (view.why) panel.append(el("p", "note", view.why));
    if (view.declared_paths?.length)
      panel.append(el("p", "note", `Declared paths: ${view.declared_paths.join(", ")}.`));
    if (view.candidates?.length)
      panel.append(el("p", "note", `The corpus holds, for example: ${view.candidates.slice(0, 12).join(", ")}.`));
  }

  function draw(view) {
    panel.hidden = false;
    panel.dataset.state = "results";
    panel.replaceChildren();

    const head = el("div", "prompt-head");
    const flag = el("span", "prompt-flag", "RETRIEVAL · NOT A GENERATION");
    const where = el("span", "note",
      view.corpus ? `Nearest entries, ranked by ${view.corpus}.`
        : "Nearest entries in the corpus. The brain does not draw an image; it finds the one its cortical state is nearest to.");
    head.append(flag, where);
    if (view.bank_entries)
      head.append(el("span", "note", `${view.bank_entries} bank entries.`));
    panel.append(head);

    if (view.empty) {
      panel.append(el("p", "note", "The endpoint returned no results for this prompt."));
      return;
    }

    // The strip: five entries, each labelled with its own cosine, the way
    // out/introspect/readout_still.png lays them out.
    const strip = el("div", "prompt-strip");
    const fractions = barFractions(view.results);
    view.results.forEach((result, i) => {
      const cell = el("div", "prompt-cell");
      if (i === 0) cell.dataset.top = "true";
      cell.append(el("div", "prompt-rank",
        `#${result.rank}  cos ${result.similarity == null ? "unavailable" : (result.similarity >= 0 ? "+" : "") + result.similarity.toFixed(3)}`));
      const frame = el("div", "prompt-thumb");
      if (result.image_url) {
        const image = document.createElement("img");
        image.src = result.image_url;
        image.alt = `Retrieved corpus entry ${result.rank}: ${result.label}`;
        frame.append(image);
      } else if (result.audio_url) {
        const audio = document.createElement("audio");
        audio.controls = true;
        audio.src = result.audio_url;
        frame.append(audio);
      } else {
        frame.append(el("span", "prompt-thumb-text", result.text || result.label));
      }
      cell.append(frame, el("div", "prompt-cell-label", result.label));
      strip.append(cell);
    });
    panel.append(strip);

    // The bars, with chance drawn on them. A retrieval whose chance level is
    // off the page is a retrieval nobody can read.
    const chart = el("div", "prompt-bars");
    view.results.forEach((result, i) => {
      const row = el("div", "prompt-bar-row");
      const fill = el("i");
      fill.style.width = `${(fractions[i] * 100).toFixed(1)}%`;
      const track = el("div", "prompt-bar-track");
      track.append(fill);
      row.append(el("span", "prompt-bar-name", result.label), track,
        el("span", "prompt-bar-value",
          result.similarity == null ? "—" : result.similarity.toFixed(4)));
      chart.append(row);
    });
    panel.append(chart);
    // The bars are cosine scaled to the top result. Saying so is not pedantry:
    // a bar that fills its track otherwise reads as "80% sure", and the number
    // it is drawn from is explicitly not a confidence.
    panel.append(el("p", "note prompt-axis",
      "Bar length is each cosine divided by the top-1 cosine. It is not a probability and it is not on the same axis as the chance level below."));

    const lines = el("dl", "prompt-figures");
    const figure = (term, value) => { lines.append(el("dt", null, term), el("dd", null, value)); };
    figure("top-1 cos", view.results[0].similarity == null ? "unavailable" : view.results[0].similarity.toFixed(4));
    figure("margin over #2", view.margin == null ? "unavailable" : (view.margin >= 0 ? "+" : "") + view.margin.toFixed(4));
    figure("chance", view.chance == null
      ? "not reported by the endpoint"
      : `${(view.chance * 100).toFixed(2)}%  (1 in ${Math.round(1 / view.chance)})`);
    panel.append(lines);

    if (view.how) panel.append(el("p", "note", view.how));
    const caveats = el("ul", "prompt-caveats");
    for (const line of view.caveats) caveats.append(el("li", null, line));
    panel.append(caveats);
  }

  bar.onsubmit = async (event) => {
    event.preventDefault();
    const typed = text.value.trim();
    if (!attached && !typed) { status("empty", "Attach an image or a sound, or type something."); return; }
    send.disabled = true;
    status("busy", "Driving the nerve and retrieving…");
    try {
      const modality = inputModality();
      const body = attached
        ? { modality, want: want.value,
            payload: modality === "audio"
              ? { audio_base64: await base64(attached) }
              : { image_base64: await base64(attached) } }
        : { modality: "text", want: want.value, payload: { text: typed } };
      const response = await post(body);
      const view = readout(response);
      if (view.error) failure(view);
      else draw(view);
    } catch (error) {
      // The endpoint may not exist yet. That is a state to name, not to hide
      // behind a spinner or to paper over with a plausible-looking result.
      if (error?.httpStatus === 404)
        status("absent", "Prompt endpoint not available.",
          "POST /api/brain/prompt is not served by this workbench build. Nothing was retrieved and nothing here is standing in for a result.");
      else status("error", error?.message || "The prompt request failed.");
    } finally {
      send.disabled = false;
    }
  };

  return {
    element: bar,
    panel,
    show(response) { draw(readout(response)); },
    clear() { panel.hidden = true; panel.replaceChildren(); },
  };
}
