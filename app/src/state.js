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
// Row-major homogeneous transform T(c+t) R F T(-c). Sparse frames are full
// snapshots for their listed entities; omitted entities return to reference.
export function bodyTransform(state, centroid) {
  const identity = [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
  if (!state) return identity;
  const finite3 = v => Array.isArray(v) && v.length === 3 && v.every(Number.isFinite);
  if (!finite3(centroid)) throw Error('Body motion needs a finite reference centroid');
  const t = state.translation_m || [0,0,0];
  const R = state.rotation_matrix || [[1,0,0],[0,1,0],[0,0,1]];
  const F = state.deformation_gradient || [[1,0,0],[0,1,0],[0,0,1]];
  if (!finite3(t) || ![R,F].every(m => Array.isArray(m) && m.length === 3 && m.every(finite3)))
    throw Error('Invalid body transform');
  const determinant = m => m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  if (determinant(R) <= 0 || determinant(F) <= 0) throw Error('Folded or collapsed body transform');
  const A = R.map(row => [0,1,2].map(j => row.reduce((v,r,k) => v+r*F[k][j],0)));
  return [...A.flatMap((row,i) => [...row, centroid[i]+t[i]-row.reduce((v,a,j)=>v+a*centroid[j],0)]), 0,0,0,1];
}

export function validBodyTrajectory(data) {
  if (!Array.isArray(data?.frames) || !data.frames.length || !data.centroids_m) return false;
  try {
    for (const [i,frame] of data.frames.entries()) {
      if (!Number.isFinite(frame.time_s) || (i && frame.time_s <= data.frames[i-1].time_s) || !frame.entities || Array.isArray(frame.entities)) return false;
      for (const [id,state] of Object.entries(frame.entities)) bodyTransform(state,data.centroids_m[id]);
      if (frame.respiration?.skin_field) {
        skinFieldParameters(frame.respiration.skin_field);
        if (!Array.isArray(frame.respiration.skin_field.entity_ids) || !frame.respiration.skin_field.entity_ids.every(id => typeof id === "string")) return false;
      }
    }
    return true;
  } catch { return false; }
}

function skinFieldParameters(field) {
  const finite = (v, n) => Array.isArray(v) && v.length === n && v.every(Number.isFinite);
  if (!finite(field.center_m, 3) || !finite(field.bounds_m?.min, 3) || !finite(field.bounds_m?.max, 3) ||
      !finite(field.reference_radii_m, 2) || !finite(field.displacement_m, 3))
    throw Error('Invalid thoracic skin field');
  const c = field.center_m, low = field.bounds_m.min, high = field.bounds_m.max;
  const [a,b] = field.reference_radii_m;
  const y = field.thorax_y_offsets_m
    ? field.thorax_y_offsets_m.map(v => v + c[1]) : field.thorax_y_m;
  const lateralSpan = Math.max(high[0]-c[0], c[0]-low[0])-a;
  if (!finite(y, 2) || !low.every((value,i) => value < high[i]) || a <= 0 || b <= 0 ||
      lateralSpan <= 0 || y[0] <= low[1] || y[1] >= high[1] || y[0] >= y[1])
    throw Error('Invalid thoracic skin support');
  return {c,low,high,a,b,y,lateralSpan};
}

// Matches the regional basis in assembly/respiration.py. Always deform the
// reference vertices, before the entity's centroid transform, without accumulation.
export function deformSkinVertices(reference, field, output = new Float32Array(reference.length)) {
  if (reference.length % 3 || output.length !== reference.length || reference === output)
    throw Error('Skin deformation requires separate reference and output coordinates');
  const params = field ? skinFieldParameters(field) : null;
  const clamp = v => Math.max(0,Math.min(1,v));
  const smooth = v => { const t = clamp(v); return t*t*(3-2*t); };
  for (let i=0;i<reference.length;i+=3) {
    const x=reference[i], y=reference[i+1], z=reference[i+2];
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) throw Error('Invalid skin reference coordinates');
    output[i]=x; output[i+1]=y; output[i+2]=z;
    if (!params) continue;
    const {c,low,high,a,b,y:thorax,lateralSpan}=params;
    const weight = smooth((y-low[1])/(thorax[0]-low[1])) * smooth((high[1]-y)/(high[1]-thorax[1])) *
      (1-smooth((Math.abs(x-c[0])-a)/lateralSpan));
    output[i] += weight*field.displacement_m[0]*Math.max(-1,Math.min(1,(x-c[0])/a));
    output[i+2] += weight*field.displacement_m[1]*clamp((z-c[2]+b)/(2*b));
  }
  return output;
}
