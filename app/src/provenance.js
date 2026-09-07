// Every structure answers where it came from, in one shape and at one prominence.
// Nothing here gates, hides, colours or orders a structure.
const esc = (value) =>
  String(value ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
const short = (hash) => (hash ? String(hash).slice(0, 12) : "—");
const verdict = (verified) =>
  verified === true ? "hash re-verified against the working tree at build time"
    : verified === false ? "hash MISMATCH against the working tree at build time"
      : "not verifiable; no file with that content is retained here";

function transformLine(t) {
  const residual = t.residual && t.residual.value !== undefined && t.residual.value !== null
    ? `${(t.residual.value * 1000).toFixed(2)} mm ${esc(t.residual.metric)}${t.residual.method ? ` · ${esc(t.residual.method)}` : ""}${t.residual.note ? `<br>${esc(t.residual.note)}` : ""}`
    : esc(t.residual_reason || "no residual measured");
  const frames = t.from || t.to ? `${esc(t.from || "?")} → ${esc(t.to || "?")}` : "";
  return `<li><strong>${esc(String(t.kind).replaceAll("_", " "))}</strong>${frames ? `<br>${frames}` : ""}<br><span class="note">${residual}</span></li>`;
}

function html(record) {
  const d = record.dataset || {}, s = record.source_file || {}, b = record.build || {}, g = record.geometry || {};
  const licence = d.license
    ? `${esc(d.license)}${d.license_evidence?.path ? `<br><span class="hash">${esc(d.license_evidence.path)} · ${short(d.license_evidence.sha256)}</span>` : ""}`
    : `<span class="note">${esc(d.license_absent_reason || "not recorded in this repository")}</span>`;
  const consumers = record.derived_artifacts || [];
  return `<section class="structure-provenance"><h4>Where this came from</h4>
    <span class="tag">${esc(record.tier)} evidence</span>
    <dl>
      <dt>Originating dataset</dt><dd>${esc(d.label || d.id)}${d.version ? ` · version ${esc(d.version)}` : ""}${d.revision ? `<br><span class="hash">${esc(d.revision)}</span>` : ""}${d.specimen ? `<br><span class="note">${esc(d.specimen)}</span>` : ""}</dd>
      <dt>Licence</dt><dd>${licence}</dd>
      <dt>Source file</dt><dd>${s.path ? `<span class="hash">${esc(s.path)}</span><br><span class="hash">${short(s.sha256)}</span> · ${esc(verdict(s.sha256_verified))}` : `<span class="note">${esc(s.absent_reason || "no upstream file recorded")}</span>${s.sha256 ? `<br><span class="hash">${short(s.sha256)}</span>` : ""}`}${s.note ? `<br><span class="note">${esc(s.note)}</span>` : ""}</dd>
      <dt>Built by</dt><dd><span class="hash">${esc(b.script || "—")}</span><br>commit <span class="hash">${short(b.commit)}</span>${b.uncommitted_changes ? " · uncommitted changes in the working tree" : ""}<br><span class="hash">script ${short(b.script_sha256)}</span></dd>
      <dt>Current geometry</dt><dd><span class="hash">${esc(g.path || "—")}</span><br><span class="hash">${short(g.sha256)}</span> · ${esc(verdict(g.sha256_verified))}${g.source_face_count ? `<br>${Number(g.source_face_count).toLocaleString()} triangles · ${esc(g.units || "")}` : ""}</dd>
      <dt>Transforms applied</dt><dd>${record.transforms?.length ? `<ul class="provenance-chain">${record.transforms.map(transformLine).join("")}</ul>` : `<span class="note">${esc(record.transform_note || "none recorded")}</span>`}<br><span class="note">Frame relation: ${esc(String(record.frame_relation || "").replaceAll("_", " "))}</span></dd>
      <dt>Evidence tier</dt><dd>${esc(record.tier)} · ${esc(record.tier_basis || "")}</dd>
      <dt>Recorded fields</dt><dd>${record.completeness?.required_present}/${record.completeness?.required_total} present${record.completeness?.missing?.length ? ` · missing ${esc(record.completeness.missing.join(", "))}` : ""}</dd>
    </dl>
    ${consumers.length ? `<details><summary>Retained builds that consumed this geometry (${consumers.length})</summary><ul>${consumers.slice(0, 24).map((c) => `<li><span class="hash">${esc(c.build)}</span> · matched by ${esc(c.matched_by)}</li>`).join("")}</ul></details>` : ""}
    ${record.assumptions?.length ? `<details><summary>Assumptions carried</summary><ul>${record.assumptions.map((a) => `<li><strong>${esc(a.id)}</strong><br>${esc(a.statement || "See the canonical assumption ledger")}</li>`).join("")}</ul></details>` : ""}
    ${record.tier_evidence && Object.keys(record.tier_evidence).length ? `<details><summary>Retained strings the tier rule read</summary><ul>${Object.entries(record.tier_evidence).map(([k, v]) => `<li><span class="hash">${esc(k)}</span><br>${esc(v)}</li>`).join("")}</ul></details>` : ""}
  </section>`;
}

export function mountProvenance(container, api) {
  let token = 0, index = null;
  const cache = new Map();
  api("/api/body/experiments/provenance").then((data) => { index = data; }).catch(() => {});
  return {
    get index() { return index; },
    async show(structure) {
      const current = ++token, target = container();
      if (!target || !structure?.id) return;
      const placeholder = document.createElement("section");
      placeholder.className = "structure-provenance";
      placeholder.innerHTML = '<h4>Where this came from</h4><p class="note">Reading the provenance record…</p>';
      target.append(placeholder);
      try {
        if (!cache.has(structure.id))
          cache.set(structure.id, await api("/api/body/experiments/provenance-" + encodeURIComponent(structure.id)));
        if (current !== token || !placeholder.isConnected) return;
        placeholder.outerHTML = html(cache.get(structure.id));
      } catch (error) {
        cache.delete(structure.id);
        if (current !== token || !placeholder.isConnected) return;
        placeholder.innerHTML = `<h4>Where this came from</h4><p class="note">Provenance unavailable: ${esc(error.message)}</p>`;
      }
    },
  };
}
