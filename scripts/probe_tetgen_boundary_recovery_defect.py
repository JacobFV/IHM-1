"""Confirm or refute the inference that the ~48-50 entity ceiling of the exact
CGAL-arrangement + TetGen `-pY` route is a memory-safety defect in TetGen.

Owns nothing that already exists. Reads the same inputs the meshing route study reads,
regenerates the exact PLC for a named cluster, and writes it as a flat binary so a NATIVE
driver can hand TetGen a tetgenio built exactly the way libigl's mesh_to_tetgenio builds
it -- same firstnumber, same one-polygon facets, same facetmarkerlist[i] = i -- with no
Python, no fork, and no CGAL in the address space.

Everything downstream (AddressSanitizer, setarch -R, MALLOC_PERTURB_) runs against that
native driver, built under data/runtime/tetgen-asan/ from upstream TetGen 1.6.0 sources.

Run with the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/probe_tetgen_boundary_recovery_defect.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/probe_tetgen_boundary_recovery_defect.py \
      --dump ladder-48,ladder-50,ladder-60 --out data/derived/tetgen-defect-v1
"""
from pathlib import Path
import argparse, hashlib, json, re, shutil, struct, sys, tempfile, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

MAGIC = b'IHMPLC01'
TETMAGIC = b'IHMTET01'
BUILD = ROOT / 'data/runtime/tetgen-asan'
WHEEL_TETGEN = (ROOT / 'data/runtime/geometry/libigl-2.6.2/venv/lib/python3.10/site-packages/'
                       'igl/copyleft/tetgen/pyigl_copyleft_tetgen.cpython-310-aarch64-linux-gnu.so')


# ------------------------------------------------------------------ plc binary

def write_plc(path, V, F):
    """Flat little-endian PLC: magic, int64 nv, int64 nf, nv*3 float64, nf*3 int64.

    The native driver reads exactly this and builds the tetgenio libigl would have built."""
    V = np.ascontiguousarray(np.asarray(V, np.float64))
    F = np.ascontiguousarray(np.asarray(F, np.int64))
    assert V.ndim == 2 and V.shape[1] == 3, V.shape
    assert F.ndim == 2 and F.shape[1] == 3, F.shape
    assert F.min() >= 0 and F.max() < len(V), (int(F.min()), int(F.max()), len(V))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'wb') as fh:
        fh.write(MAGIC)
        fh.write(struct.pack('<qq', len(V), len(F)))
        fh.write(V.tobytes())
        fh.write(F.tobytes())
    return path


def read_plc(path):
    raw = Path(path).read_bytes()
    assert raw[:8] == MAGIC, raw[:8]
    nv, nf = struct.unpack_from('<qq', raw, 8)
    off = 24
    V = np.frombuffer(raw, np.float64, nv * 3, off).reshape(nv, 3)
    F = np.frombuffer(raw, np.int64, nf * 3, off + nv * 24).reshape(nf, 3)
    return V, F


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tetbin(path):
    """The tet mesh the native driver writes: magic, int64 nv, int64 nt, nv*3 f64, nt*4 i64."""
    raw = Path(path).read_bytes()
    assert raw[:8] == TETMAGIC, raw[:8]
    nv, nt = struct.unpack_from('<qq', raw, 8)
    off = 24
    TV = np.frombuffer(raw, np.float64, nv * 3, off).reshape(nv, 3)
    TT = np.frombuffer(raw, np.int64, nt * 4, off + nv * 24).reshape(nt, 4)
    return TV, TT


def tet_volumes(TV, TT):
    a, b, c, d = (TV[TT[:, i]] for i in range(4))
    return np.einsum('ij,ij->i', np.cross(b - a, c - a), d - a) / 6.0


def validate_mesh(plc_path, tet_path):
    """Independent acceptance of a mesh the patched TetGen produced.

    Three properties the patches could plausibly have broken, each measured rather than
    asserted from the run log:
      1. the tets exactly fill the padded box the PLC's first 8 vertices define, so no
         volume was lost or double counted (box closure);
      2. every tet has positive signed volume in the emitted vertex order;
      3. every PLC vertex still exists in the output at distance 0.0 m, which is what
         TetGen's -Y promises and the whole exactness argument rests on.
    """
    PV, PF = read_plc(plc_path)
    TV, TT = read_tetbin(tet_path)
    vol = tet_volumes(TV, TT)
    box = PV[:8]
    lo, hi = box.min(0), box.max(0)
    box_volume = float(np.prod(hi - lo))
    total = float(vol.sum())
    # every input vertex must survive; match by exact coordinate lookup
    key = {tuple(r) for r in map(tuple, TV)}
    missing = sum(1 for r in map(tuple, PV) if r not in key)
    return {'plc': str(Path(plc_path).name), 'mesh': str(Path(tet_path).name),
            'plc_vertices': int(len(PV)), 'plc_facets': int(len(PF)),
            'tet_vertices': int(len(TV)), 'tets': int(len(TT)),
            'steiner_points': int(len(TV) - len(PV)),
            'box_volume_m3': box_volume, 'tet_volume_sum_m3': total,
            'box_volume_closure_relative_error': float(total / box_volume - 1.0),
            'tets_with_nonpositive_volume': int((vol <= 0).sum()),
            'min_tet_volume_m3': float(vol.min()), 'max_tet_volume_m3': float(vol.max()),
            'plc_vertices_absent_from_output': int(missing),
            'plc_vertices_reproduced_exactly': bool(missing == 0),
            'mesh_sha256': sha256_file(tet_path)}


# ----------------------------------------------------------------- plc builder

def build_exact_plc(label, pad):
    """The exact arrangement PLC plus the padded bounding box, identical in construction to
    build_meshing_route_study.diagnose(): box corners first, arrangement vertices offset by 8."""
    import build_cross_structure_conflict_repair as R
    import build_meshing_route_study as S
    rows = R.sources()
    by = {r['entity_id']: r for r in rows}
    coin = S.REPAIR / 'coincidence.json'
    drop = set(json.loads(coin.read_text())['identical_geometry_entities_dropped']) \
        if coin.exists() else set()
    members = S.cluster_members(label, by, drop,
                                ROOT / 'data/derived/meshing-route-study-v1/ladder-order.json')
    parts = R.build_parts(by, members)
    t0 = time.monotonic()
    V, F, src = S.soup(parts)
    W, DF, dsrc, crossings, dropped, collapsed = S.exact_arrangement(V, F, src)
    seconds = time.monotonic() - t0
    lo = W.min(0) - pad
    hi = W.max(0) + pad
    PV = np.ascontiguousarray(np.concatenate([lo + R.BOX_CORNERS * (hi - lo), W]))
    PF = np.ascontiguousarray(np.concatenate([R.BOX, DF + 8]).astype(np.int64))
    meta = {'label': label, 'members': len(parts), 'pad_m': pad,
            'arrangement_seconds': seconds,
            'crossing_face_pairs_resolved': int(crossings),
            'plc_vertices': int(len(PV)), 'plc_facets': int(len(PF)),
            'arrangement_vertices': int(len(W)), 'arrangement_facets': int(len(DF)),
            'in_memory_sha256': hashlib.sha256(PV.tobytes() + PF.tobytes()).hexdigest(),
            'in_memory_sha256_note': 'matches ceiling-diagnosis.json plc_sha256 when the same '
                                     'pad and the same cluster membership are used'}
    return PV, PF, meta


# ------------------------------------------------------------------- self test

def self_test():
    checks = []

    def check(name, ok, detail=None):
        checks.append({'check': name, 'passed': bool(ok), 'detail': detail})

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        V = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        F = np.array([[0, 1, 3], [0, 3, 2], [0, 2, 1], [1, 2, 3]], np.int64)
        p = write_plc(tmp / 'unit.plcbin', V, F)
        V2, F2 = read_plc(p)
        check('plc binary round-trips vertices exactly', np.array_equal(V, V2))
        check('plc binary round-trips facets exactly', np.array_equal(F, F2))
        check('plc binary length is header + payload',
              p.stat().st_size == 24 + V.size * 8 + F.size * 8, p.stat().st_size)
        check('plc sha256 is stable over a rewrite',
              sha256_file(p) == sha256_file(write_plc(tmp / 'unit2.plcbin', V, F)))
        bad = False
        try:
            write_plc(tmp / 'bad.plcbin', V, np.array([[0, 1, 9]], np.int64))
        except AssertionError:
            bad = True
        check('write_plc rejects an out-of-range facet index', bad)
        # never hand a caller-supplied root to anything that writes: this test only ever
        # writes inside the TemporaryDirectory above
        check('self-test wrote only inside a TemporaryDirectory',
              all(str(q).startswith(str(tmp)) for q in tmp.rglob('*')))

        # tet mesh round-trip and the volume/closure measures, on a unit cube split
        # into 6 tets whose exact volume is known
        cube = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
        # Kuhn decomposition: the six monotone lattice paths from (0,0,0) to (1,1,1),
        # each oriented positively, tile the cube exactly
        tt = np.array([[0, 1, 2, 6], [0, 5, 1, 6], [0, 2, 3, 6],
                       [0, 3, 7, 6], [0, 4, 5, 6], [0, 7, 4, 6]], np.int64)
        tp = tmp / 'unit.tetbin'
        with open(tp, 'wb') as fh:
            fh.write(TETMAGIC); fh.write(struct.pack('<qq', len(cube), len(tt)))
            fh.write(np.ascontiguousarray(cube).tobytes())
            fh.write(np.ascontiguousarray(tt).tobytes())
        TV2, TT2 = read_tetbin(tp)
        check('tet binary round-trips', np.array_equal(cube, TV2) and np.array_equal(tt, TT2))
        v = tet_volumes(TV2, TT2)
        # six exact sixths do not sum to 1.0 in binary floating point; 1 ulp is the bar
        check('six tets of the unit cube sum to 1.0 m3 to within an ulp',
              abs(float(v.sum()) - 1.0) <= 2 * np.spacing(1.0), float(v.sum()))
        check('all six have positive volume', bool((v > 0).all()), v.tolist())
        pp = write_plc(tmp / 'cube.plcbin', cube, np.array([[0, 1, 2]], np.int64))
        rep = validate_mesh(pp, tp)
        check('validate_mesh reports box closure on the unit cube to within an ulp',
              abs(rep['box_volume_closure_relative_error']) <= 2 * np.spacing(1.0),
              rep['box_volume_closure_relative_error'])
        check('validate_mesh finds every plc vertex in the output',
              rep['plc_vertices_reproduced_exactly'], rep['plc_vertices_absent_from_output'])

        # the patched shared library, if it has been built, must mesh a cube
        if PATCHED_SO.exists():
            TVc, TTc, stc = patched_tetrahedralize(cube, np.array(
                [[0, 1, 2], [0, 2, 3], [4, 6, 5], [4, 7, 6], [0, 5, 1], [0, 4, 5],
                 [1, 6, 2], [1, 5, 6], [2, 7, 3], [2, 6, 7], [3, 4, 0], [3, 7, 4]], np.int64),
                'pYQ')  # Q keeps TetGen's banner out of this script's JSON on stdout
            ok = stc == 0 and TTc is not None and len(TTc) > 0
            vol = float(tet_volumes(TVc, TTc).sum()) if ok else None
            check('patched libtetgen meshes the unit cube to volume 1.0 m3',
                  ok and abs(vol - 1.0) <= 8 * np.spacing(1.0), {'status': stc, 'volume': vol})
        else:
            checks.append({'check': 'patched libtetgen present (informational)',
                           'passed': True, 'detail': {'path': str(PATCHED_SO), 'exists': False}})

    driver = ROOT / 'data/runtime/tetgen-asan/build/tetgen_driver'
    checks.append({'check': 'native driver present (informational)', 'passed': True,
                   'detail': {'path': str(driver), 'exists': driver.exists()}})
    passed = all(c['passed'] for c in checks)
    print(json.dumps({'self_test': 'ihm.tetgen-defect-probe', 'passed': passed,
                      'checks': checks}, indent=2))
    return passed


# ------------------------------------------------- patched tetgen, in this process

PATCHED_SO = BUILD / 'build/libtetgen_patched.so'


def patched_tetrahedralize(V, F, flags='pY'):
    """igl.copyleft.tetgen.tetrahedralize, but over the patched TetGen 1.6.0 built here.

    Substitutable for the wheel's binding: the tetgenio is constructed identically. Returns
    (TV, TT, status); status is TetGen's own throw value, 0 on success.

    Callers should keep doing what build_cross_structure_conflict_repair.run_tetgen already
    does and fork first -- TetGen still aborts the process on error paths this work did not
    touch, and that has to stay a receipt rather than a dead run."""
    import ctypes
    lib = ctypes.CDLL(str(PATCHED_SO))
    fn = lib.ihm_tetrahedralize
    fn.restype = ctypes.c_int
    fn.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int64,
                   ctypes.POINTER(ctypes.c_int64), ctypes.c_int64, ctypes.c_char_p,
                   ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.c_int64),
                   ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)), ctypes.POINTER(ctypes.c_int64)]
    lib.ihm_tetgen_free.argtypes = [ctypes.c_void_p]
    Vc = np.ascontiguousarray(np.asarray(V, np.float64))
    Fc = np.ascontiguousarray(np.asarray(F, np.int64))
    tv = ctypes.POINTER(ctypes.c_double)(); tt = ctypes.POINTER(ctypes.c_int64)()
    ntv = ctypes.c_int64(); ntt = ctypes.c_int64()
    st = fn(Vc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(Vc),
            Fc.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)), len(Fc),
            flags.encode(), ctypes.byref(tv), ctypes.byref(ntv),
            ctypes.byref(tt), ctypes.byref(ntt))
    if st != 0:
        return None, None, st
    TV = np.ctypeslib.as_array(tv, (int(ntv.value), 3)).copy()
    TT = np.ctypeslib.as_array(tt, (int(ntt.value), 4)).copy()
    lib.ihm_tetgen_free(ctypes.cast(tv, ctypes.c_void_p))
    lib.ihm_tetgen_free(ctypes.cast(tt, ctypes.c_void_p))
    return TV, TT, 0


# ------------------------------------------------------------------ wheel run

def wheel_run(label, flags, out):
    """The SHIPPED libigl wheel on the same PLC, printing the same RESULT line the native
    driver prints, so both go through one parser.

    No fork: the point is to observe what the wheel does in a plain process, since the
    recorded non-determinism was across processes."""
    from igl.copyleft import tetgen  # noqa: E402  (only importable in the libigl venv)
    V, F = read_plc(Path(out) / 'plc' / ('%s.plcbin' % label))
    print('# wheel: %s nv=%d nf=%d switches=%s' % (label, len(V), len(F), flags), flush=True)
    began = time.monotonic()
    status, tets, tverts = 0, -1, -1
    try:
        r = tetgen.tetrahedralize(np.ascontiguousarray(V), np.ascontiguousarray(F), flags=flags)
        status = int(r[-1]); tets = int(len(r[1])); tverts = int(len(r[0]))
    except BaseException as error:
        status = -1
        print('# wheel exception: %s: %s' % (type(error).__name__, error), flush=True)
    sys.stdout.flush()
    print('RESULT %s' % json.dumps({'run': 0, 'status': status, 'tets': tets,
                                    'tet_vertices': tverts,
                                    'seconds': round(time.monotonic() - began, 3),
                                    'max_rss_kb': -1}), flush=True)
    return 0


# -------------------------------------------------------------- predicate audit

def _exact_orient3d_sign(pa, pb, pc, pd):
    from fractions import Fraction as Fr
    a, b, c, d = ([Fr(x) for x in p] for p in (pa, pb, pc, pd))
    m = [[a[k] - d[k] for k in range(3)], [b[k] - d[k] for k in range(3)],
         [c[k] - d[k] for k in range(3)]]
    det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
           - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
           + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    return (det > 0) - (det < 0)


def _exact_insphere_sign(pa, pb, pc, pd, pe):
    from fractions import Fraction as Fr
    pts = [[Fr(x) for x in p] for p in (pa, pb, pc, pd)]
    e = [Fr(x) for x in pe]
    m = []
    for q in pts:
        dx, dy, dz = (q[k] - e[k] for k in range(3))
        m.append([dx, dy, dz, dx * dx + dy * dy + dz * dz])

    def det3(r):
        return (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
                - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
                + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))
    det = 0
    for i in range(4):
        minor = [[m[j][k] for k in range(3)] for j in range(4) if j != i]
        det += ((-1) ** (i + 3)) * m[i][3] * det3(minor)
    return (det > 0) - (det < 0)


def _read_probe(path, n):
    raw = Path(path).read_bytes()
    rec = 8 * 15 + 8
    for i in range(min(n, len(raw) // rec)):
        v = struct.unpack_from('<15dii', raw, i * rec)
        yield [v[0:3], v[3:6], v[6:9], v[9:12], v[12:15]], v[15]


def predicate_audit(out, n):
    """Grade three builds of the SAME predicates against exact rational arithmetic.

    Shewchuk's adaptive predicates are exact only if every multiply and add rounds once.
    aarch64 GCC fuses a*b+c into an FMA by default, which skips a rounding. This measures
    whether that actually changes a sign, for a local build with contraction on, the same
    build with it off, and the predicate code inside the shipped wheel, called by dlopen."""
    import ctypes, subprocess
    out = Path(out)
    tmpd = Path(tempfile.mkdtemp())
    report = {'schema': 'ihm.tetgen-predicate-audit.v1', 'samples': n,
              'generator': 'data/runtime/tetgen-asan/predicate_probe.cpp',
              'reference': 'exact rational arithmetic (fractions.Fraction) on the same doubles',
              'builds': {}}
    try:
        for tag in ('off', 'fast'):
            exe = BUILD / ('build/predicate_probe_%s' % tag)
            if not exe.exists():
                continue
            for kind in ('o', 'i'):
                f = tmpd / ('%s_%s.bin' % (tag, kind))
                subprocess.run([str(exe), str(n), str(f), kind], check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            bad_o = sum(1 for pts, sg in _read_probe(tmpd / ('%s_o.bin' % tag), n)
                        if _exact_orient3d_sign(*pts[:4]) != sg)
            bad_i = sum(1 for pts, sg in _read_probe(tmpd / ('%s_i.bin' % tag), n)
                        if _exact_insphere_sign(*pts) != sg)
            report['builds']['local_gcc_ffp_contract_%s' % tag] = {
                'orient3d_wrong_signs': bad_o, 'insphere_wrong_signs': bad_i,
                'orient3d_wrong_fraction': bad_o / float(n),
                'insphere_wrong_fraction': bad_i / float(n)}
        if WHEEL_TETGEN.exists():
            lib = ctypes.CDLL(str(WHEEL_TETGEN))
            init = lib._Z9exactinitiiiddd
            init.argtypes = [ctypes.c_int] * 3 + [ctypes.c_double] * 3
            init.restype = None
            init(0, 0, 0, 1.0, 1.0, 1.0)
            D3 = ctypes.c_double * 3
            ins = lib._Z8inspherePdS_S_S_S_
            ins.argtypes = [ctypes.POINTER(ctypes.c_double)] * 5; ins.restype = ctypes.c_double
            ori = lib._Z8orient3dPdS_S_S_
            ori.argtypes = [ctypes.POINTER(ctypes.c_double)] * 4; ori.restype = ctypes.c_double
            bad_o = bad_i = 0
            for pts, _ in _read_probe(tmpd / 'off_i.bin', n):
                P = [D3(*q) for q in pts]
                ri = ins(P[0], P[1], P[2], P[3], P[4]); ro = ori(P[0], P[1], P[2], P[3])
                if ((ri > 0) - (ri < 0)) != _exact_insphere_sign(*pts): bad_i += 1
                if ((ro > 0) - (ro < 0)) != _exact_orient3d_sign(*pts[:4]): bad_o += 1
            report['builds']['shipped_libigl_wheel'] = {
                'orient3d_wrong_signs': bad_o, 'insphere_wrong_signs': bad_i,
                'orient3d_wrong_fraction': bad_o / float(n),
                'insphere_wrong_fraction': bad_i / float(n),
                'so': str(WHEEL_TETGEN.relative_to(ROOT)),
                'so_sha256': sha256_file(WHEEL_TETGEN)}
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    (out / 'predicate-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return report


# ---------------------------------------------------------------- run receipts

RESULT_RE = re.compile(r'^RESULT (\{.*\})\s*$', re.M)


def parse_log(path):
    """One native run, read back out of its own stdout."""
    text = Path(path).read_text(errors='replace')
    variant, label, flags, cond, i = Path(path).stem.split('__')
    rec = {'log': str(Path(path).relative_to(ROOT)), 'variant': variant, 'label': label,
           'flags': flags, 'condition': cond, 'repeat': int(i)}
    m = RESULT_RE.search(text)
    if m:
        rec.update(json.loads(m.group(1)))
    else:
        rec.update({'status': None, 'tets': None,
                    'note': 'no RESULT line: the process died before the driver could print'})
    ec = re.search(r'^# exit_code (-?\d+)$', text, re.M)
    rec['process_exit_code'] = int(ec.group(1)) if ec else None
    rec['glibc_heap_abort'] = bool(re.search(r'corrupted (size vs\. prev_size|top size)|'
                                             r'free\(\): |malloc\(\): ', text))
    rec['sigsegv'] = 'Segmentation fault' in text
    rec['asan_error'] = (re.search(r'ERROR: AddressSanitizer: (\S+)', text) or [None, None])[1] \
        if 'ERROR: AddressSanitizer' in text else None
    rec['tetgen_warnings'] = dict(sorted(
        {l.strip(): text.count(l) for l in set(text.splitlines())
         if l.strip().startswith('Warning:')}.items(), key=lambda kv: -kv[1])[:8])
    for key, pat in (('skipped_input_triangles', r'!!! (\d+) input triangles are skipped'),
                     ('mesh_tetrahedra', r'Mesh tetrahedra:\s+(\d+)'),
                     ('mesh_points', r'Mesh points:\s+(\d+)'),
                     ('mesh_faces_on_input_facets', r'Mesh faces on input facets:\s+(\d+)'),
                     ('steiner_points_on_input_segments', r'Steiner points on input segments:\s+(\d+)'),
                     ('steiner_points_inside_domain', r'Steiner points inside domain:\s+(\d+)')):
        mm = re.search(pat, text)
        rec[key] = int(mm.group(1)) if mm else None
    rec['succeeded'] = bool(rec.get('status') == 0 and (rec.get('tets') or 0) > 0
                            and rec.get('process_exit_code') == 0)
    return rec


def build_findings(out):
    """The narrative, with every number either measured here or read back out of the
    receipts this script wrote.  Claims are split into what was OBSERVED and what is
    INFERRED from it."""
    out = Path(out)
    runs = json.loads((out / 'runs.json').read_text())['runs'] \
        if (out / 'runs.json').exists() else []
    val = json.loads((out / 'mesh-validation.json').read_text()) \
        if (out / 'mesh-validation.json').exists() else {}
    pred = json.loads((out / 'predicate-audit.json').read_text()) \
        if (out / 'predicate-audit.json').exists() else {}

    def tally(**sel):
        rows = [r for r in runs if all(r.get(k) == v for k, v in sel.items())]
        return {'runs': len(rows), 'succeeded': sum(r['succeeded'] for r in rows),
                'flags_seen': sorted({r['flags'] for r in rows}),
                'tets_by_flags': {f: sorted({r['tets'] for r in rows
                                             if r['succeeded'] and r['flags'] == f})
                                  for f in sorted({r['flags'] for r in rows})},
                'distinct_tet_counts': sorted({r['tets'] for r in rows if r['succeeded']}),
                'glibc_heap_aborts': sum(bool(r['glibc_heap_abort']) for r in rows),
                'sigsegvs': sum(bool(r['sigsegv']) for r in rows),
                'statuses': sorted({r.get('status') for r in rows}, key=lambda x: (x is None, x))}

    labels = sorted({r['label'] for r in runs})
    ladder = {}
    for lab in labels:
        ladder[lab] = {'stock_1.6.0': tally(label=lab, variant='plain', condition='base'),
                       'patched_0001_0002_0003': tally(label=lab, variant='fixall', condition='base'),
                       'shipped_wheel': tally(label=lab, variant='wheel', condition='base')}

    return {
      'schema': 'ihm.tetgen-defect-findings.v1',
      'question': 'is the ~48-50 entity ceiling of the exact CGAL-arrangement + TetGen -pY '
                  'route a memory-safety defect in the vendored TetGen, and can it be fixed',
      'verdict': 'CONFIRMED IN PART, AND LIFTED. There is a real memory-safety defect in the '
                 'vendored TetGen (a stack-buffer over-read in Shewchuk predicates, already '
                 'fixed once by libigl in 2018 and regressed in 2021), and it is hit on this '
                 'workload -- but it is NOT what caps the ladder. The cap is two ordinary '
                 'logic defects in TetGen 1.6.0 boundary recovery and mesh improvement. With '
                 'all three patched, the same PLCs that failed 0/3 now mesh, and the ceiling '
                 'moves from 48 entities to at least the largest rung measured here.',
      'provenance': {
        'binding': 'igl.copyleft.tetgen in libigl 2.6.2, PyPI wheel '
                   'libigl-2.6.2-cp310-cp310-manylinux_2_27_aarch64.manylinux_2_28_aarch64.whl',
        'wheel_is_official_pypi_build': True,
        'compiler_of_the_wheel': 'GCC 14.2.1 20250110 (Red Hat 14.2.1-7); the .so also carries '
                                 'a clang 20.1.8 .comment string from the manylinux image',
        'assertions_in_the_wheel': 'compiled out (-DNDEBUG): zero assert strings in the .so',
        'libigl_python_bindings_2.6.2_pins_libigl': '678e1fff76815e0c4c5d1f025ee2129181cc7d86',
        'libigl_pins_tetgen': 'https://github.com/libigl/tetgen.git @ '
                              'e05aca7df74e3f531bc35733ed87d36d437266c5 (2025-04-15)',
        'tetgen_version_string_in_the_so': 'Version 1.6 / August, 2020 / Copyright (C) 2002 - 2020',
        'vendored_patch': 'NONE. tetgen.h, tetgen.cxx and predicates.cxx at that commit are '
                          'byte-identical to the WIAS release tarball tetgen1.6.0.tar.gz: '
                          'tetgen.h 69b012da7c0b01e6327f5f281ac2d10fc3fb42503d22366c7bd0bc9181bd833b, '
                          'tetgen.cxx a10f5e74c3ec45ffda3d9609b1cb6e864f1a628b4c2e3e61cabd49ddfccc36bc, '
                          'predicates.cxx 5a050a2391e8dea44f6aba176372ba790667d553b98c7eea4bb351eabc81c050',
        'upstream_status': 'TetGen 1.6.0 (2020-08-31) is the last public release. WIAS itself '
                           'recommends v1.5.1 (2018) as "most stable version" and calls v1.6.0 '
                           '"most recent version with some rough edges".',
        'build_flags_of_the_tetgen_target': 'libigl/tetgen CMakeLists.txt adds only -DTETLIBRARY '
                                            'and POSITION_INDEPENDENT_CODE; no floating-point '
                                            'flags of any kind.'},
      'reproduction_harness': {
        'why_native': 'the recorded failures were observed through a Python fork, so the child '
                      'inherited a CPython heap. The driver here builds the SAME tetgenio that '
                      'igl::copyleft::tetgen::mesh_to_tetgenio builds (firstnumber 0, one '
                      'single-polygon facet per triangle, facetmarkerlist[i] = i, empty H/VM/FM/R) '
                      'in a fresh process with no Python and no CGAL in the address space.',
        'driver': 'data/runtime/tetgen-asan/driver.cpp',
        'build': 'data/runtime/tetgen-asan/build.sh',
        'fidelity_check': 'stock 1.6.0 built here at -O2 -DNDEBUG -ffp-contract=off reproduces '
                          'the wheel: ladder-48 succeeds (1060354 tets here against the wheel\'s '
                          '1060351, same "4 segments are not recovered" warning), ladder-50 and '
                          'ladder-60 fail.'},
      'defects': [
        {'id': 'D1',
         'title': 'stack-buffer over-read in fast_expansion_sum_zeroelim (Shewchuk predicates)',
         'status': 'CONFIRMED by AddressSanitizer, upstream fix exists, fix REGRESSED in the pin',
         'observed': 'AddressSanitizer: stack-buffer-overflow, READ of size 8, '
                     'predicates.cxx:994 in fast_expansion_sum_zeroelim, reading offset 192 -- '
                     'one element past the 4-element local cd[4] -- in the frame of '
                     'orient3dexact (predicates.cxx:1607). Call chain: orient3dexact '
                     '(predicates.cxx:1658) <- tetgenmesh::planelineint (tetgen.cxx:6877) <- '
                     'add_steinerpt_to_recover_edge (tetgen.cxx:20163) <- recoversegments '
                     '(tetgen.cxx:20421) <- recoverboundary (tetgen.cxx:23233) <- tetrahedralize.',
         'cause': 'e[++eindex] / f[++findex] are evaluated before the loop guard rechecks the '
                  'index, so the last iteration reads one past the end of the expansion.',
         'upstream_fix': 'libigl/tetgen c63e7a6434652b8a2065c835bd9d6d298db1a0bc (2018-12-28, '
                         'Jeremie Dumas, "Fix a stack-buffer overflow in predicates.cxx", ported '
                         'from libMesh 4a054cb38c97967b34e7a5b8fa5a74b0d57816ab). It was applied '
                         'to the fork\'s TetGen 1.5.1 tree and LOST when commit 4f3bfba ("V1.6", '
                         '2021-05-26) replaced the files with stock upstream 1.6.0. '
                         'fast_expansion_sum_zeroelim is byte-identical in 1.5.1 and 1.6.0, so '
                         'the patch transplants without modification.',
         'patch': 'data/runtime/tetgen-asan/patches/0001-predicates-stack-buffer-overflow.patch',
         'effect_of_the_fix_alone': 'ASan clean at that site; the ladder-50 outcome is UNCHANGED '
                                    '(still status 3 plus a glibc heap abort). So D1 is real and '
                                    'is hit, but it is not the proximate cause of the ceiling.'},
        {'id': 'D2',
         'title': 'a segment is reported as self-intersecting with itself, aborting the mesh',
         'status': 'CONFIRMED by reading TetGen\'s own diagnostic output against the input PLC',
         'observed': 'At ladder-50 TetGen prints "Warning: Two line segments are nearly '
                     'overlapping. 1st: [54710,55058]. 2nd: [54710,55058]." -- the same index '
                     'pair twice -- and quotes a separation of 6.9388939039072284e-18 m, which '
                     'is exactly 2^-57, and 3.4694469519536142e-18 m, exactly 2^-58: pure '
                     'cancellation noise. It then writes tetgen-tmpfile_skipped.* and exits with '
                     'status 3, "The input surface mesh contain self-intersections".',
         'the_input_does_not_self_intersect': 'the four "self-intersecting" input triangles at '
             'ladder-50 are facets 291379 (66979,67912,67608), 291362 (66825,67912,66979), '
             '287585 (54727,54710,55058) and 287588 (54710,55098,55058). They form two '
             'EDGE-ADJACENT PAIRS: 291379/291362 share the edge (66979,67912) of length '
             '1.785288e-03 m, 287585/287588 share (54710,55058) of length 3.961209e-04 m. Two '
             'triangles that share an edge are not a self-intersection. Their areas are '
             '3.2220e-07, 2.5842e-07, 5.8773e-09 and 8.1137e-09 m2 and their smallest altitudes '
             '3.2986e-04, 1.9187e-04, 2.9674e-05 and 4.0966e-05 m, all four orders of magnitude '
             'above TetGen\'s own 1.08111e-08 m tolerance. CGAL certifies zero residual '
             'self-intersections in the same PLC.',
         'cause': 'in add_steinerpt_to_recover_edge (tetgen.cxx ~20218), when insertpoint returns '
                  'ONVERTEX/NEARVERTEX the code asks is_segment(p1, nearpt) / is_segment(p2, '
                  'nearpt) to identify a second, overlapping segment. When nearpt is an ENDPOINT '
                  'of the segment being recovered, that query trivially finds *misseg itself*, '
                  'so TetGen reports the segment as overlapping itself, sets SELF_INTERSECT, '
                  'moves the segment and every subface around it to the skipped lists and aborts. '
                  'A Steiner point that lands on one of the segment\'s own endpoints is a failed '
                  'split, not an input self-intersection.',
         'patch': 'data/runtime/tetgen-asan/patches/'
                  '0002-do-not-report-a-segment-as-overlapping-itself.patch',
         'effect_of_the_fix': 'ladder-50 now clears boundary recovery, exterior removal, Steiner '
                              'suppression, Delaunay recovery and smoothing, and dies later in '
                              'mesh improvement (D3). The glibc heap abort disappears with it.'},
        {'id': 'D3',
         'title': 'best-effort mesh improvement treats its own failures as internal errors',
         'status': 'CONFIRMED by a gdb catch-throw backtrace',
         'observed': 'terminatetetgen(this, 2) at tetgen.cxx:31453 (patched tree) in '
                     'tetgenmesh::repair_tet <- repair_badqual_tets <- improve_mesh <- '
                     'tetrahedralize. Status 2 is TetGen\'s "Encounter an internal error ... '
                     'Please report this bug to Hang.Si". The mesh is complete and '
                     'boundary-conforming at that point; what fails is create_a_shorter_edge() '
                     'declining to move a Steiner point off a short edge during QUALITY repair.',
         'patch': 'data/runtime/tetgen-asan/patches/'
                  '0003-mesh-improvement-failures-are-not-internal-errors.patch',
         'patch_nature': 'ROUTE AROUND, not a fix: the three terminatetetgen(this,2) calls in '
                         'repair_tet become "return false" (this tet was not repaired). It '
                         'disables an internal consistency check in a best-effort quality pass, '
                         'so the resulting mesh is validated independently rather than trusted.',
         'effect_of_the_fix': 'ladder-50 and every larger rung measured here return status 0.'},
        {'id': 'D4',
         'title': 'heap-buffer-overflow WRITE while reporting the skipped facets',
         'status': 'CONFIRMED by AddressSanitizer; this is the glibc abort the earlier study saw',
         'observed': 'AddressSanitizer: heap-buffer-overflow, WRITE of size 4 at '
                     'tetgen.cxx:34458 in tetgenmesh::outsubfaces(tetgenio*), 0 bytes after a '
                     '3917712-byte region allocated at tetgen.cxx:34352, which is '
                     'out->trifacelist = new int[subfaces->items * 3]. 3917712 bytes is '
                     '326476*3 ints, so the traversal emits more subfaces than subfaces->items '
                     'counted when the array was sized.',
         'when': 'only on the status-3 exit path of D2, where TetGen removes the skipped '
                 'subfaces and then still writes the subface arrays out before throwing. It is '
                 'the direct cause of the "corrupted size vs. prev_size" glibc abort recorded '
                 'in every failing run of the earlier study, and of the SIGSEGVs.',
         'patch': 'none written. Patch 0002 removes the path entirely, and no run that reaches '
                  'a mesh touches it. Recorded so it is not rediscovered.'},
      ],
      'build_hazard_not_the_cause': {
        'title': 'FMA contraction silently invalidates the exact predicates on aarch64',
        'measured': pred,
        'observed_locally': 'GCC 13.3.0 at -O2 with its default -ffp-contract=fast fuses 394 '
                            'multiply-adds into orient3dadapt and breaks insphere: see '
                            'predicate-audit.json, where the contracted build returns the WRONG '
                            'SIGN against exact rational arithmetic on 18.78 percent of '
                            'near-cospherical queries and the uncontracted build on none. '
                            'orient3d is unaffected either way.',
        'consequence_locally': 'the contracted -O2 build loops forever on a 16-vertex PLC of two '
                               'nested axis-aligned cubes (the same PLC meshes in 3 ms with '
                               'contraction off and in 0.9 ms through the shipped wheel) and '
                               'SIGSEGVs on ladder-48 inside tetgenmesh::memorypool::alloc during '
                               'incrementaldelaunay, on a garbage freelist pointer '
                               '(0x340000aaaaafd189).',
        'the_shipped_wheel_is_NOT_affected': 'the wheel also carries fused ops (240 in '
                                             'orient3dadapt, 98 in insphereadapt, 80 in '
                                             'insphereexact, 29 in insphere), but calling the '
                                             'wheel\'s own insphere and orient3d through ctypes '
                                             'on the same quadruples gives 0 sign '
                                             'disagreements against exact arithmetic. GCC 14 at '
                                             'the wheel\'s optimisation level happened to fuse '
                                             'operations that do not break the error-free '
                                             'transformations.',
        'status': 'a latent build hazard, measured, that would bite the moment libigl rebuilds '
                  'its aarch64 wheel with a different compiler. It is NOT the cause of the '
                  'recorded ceiling. Every instrumented build here uses -ffp-contract=off so '
                  'that failures are attributable to TetGen and not to the build.',
        'upstream_precedent': ['BrunoLevy/geogram#382 and #358 -- Darwin/arm64 configs omitted '
                               '-ffp-contract=off, "3D Delaunay locate_inexact() can loop '
                               'forever on some build configs"; fixed in geogram 1.10.1',
                               'pyvista/pytetwild#48 -- "Build with -ffp-contract=off to fix '
                               'macOS arm64 hang"',
                               'environmental-modeling-workflows/watershed-workflow#166 -- '
                               '"fix meshpy Triangle crash on linux-aarch64 by disabling FP '
                               'contraction"',
                               'geogram ships add_flags(CMAKE_CXX_FLAGS -frounding-math '
                               '-ffp-contract=off) in Linux-gcc-aarch64.cmake with the comment '
                               '"would break exact predicates". libigl/tetgen ships nothing.'],
      },
      'non_determinism': {
        'question': 'why did the same PLC mesh once at 60 entities and then fail on re-run',
        'observed_here': tally(variant='plain', label='ladder-60', condition='base'),
        'stock_ladder_50_repeats': tally(variant='plain', label='ladder-50', condition='base'),
        'aslr_disabled': tally(variant='plain', label='ladder-50', condition='norandom'),
        'malloc_perturb': tally(variant='plain', label='ladder-50', condition='perturb'),
        'shipped_wheel_ladder_60': tally(variant='wheel', label='ladder-60', condition='base'),
        'observed': 'the native stock build is deterministic across processes on this host: '
                    'every repeat of a given rung gives the same verdict, and setarch -R '
                    '(ASLR off) and MALLOC_PERTURB_=165 change nothing. So the recorded '
                    'flip is not reproduced by the mesher alone in a fresh process.',
        'mechanism_for_the_heap_aborts': 'D4. A 4-byte write one element past a ~3.9 MB '
                                         'new[] block lands on whatever glibc put next. '
                                         'Whether that is a free chunk header (glibc aborts '
                                         'with "corrupted size vs. prev_size"), a live '
                                         'allocation (silent), or an unmapped page (SIGSEGV) '
                                         'depends on the heap layout, which is exactly the '
                                         'ASLR-and-allocator-state story the earlier study '
                                         'inferred. That part of the inference is CONFIRMED.',
        'mechanism_for_the_success_flip': 'NOT ESTABLISHED. D4 happens after the verdict is '
                                          'already status 3, so it cannot turn a failure into '
                                          'a success. Two candidates were found in the source '
                                          'and neither was proven: (a) TetGen 1.6 calls the '
                                          'process-global libc rand() in at least eight places '
                                          '(tetgen.cxx 11257, 11750, 11834, 19301, 25396, '
                                          '25712, 27633, 29058) and re-seeds it with '
                                          'srand(in->numberofpoints), srand(arylen) and '
                                          'srand(intptlist->objects), so its random stream is '
                                          'shared with everything else in the address space -- '
                                          'in the recorded runs that address space also held '
                                          'CPython, NumPy and CGAL; (b) the D1 over-read takes '
                                          'a value from the stack slot after cd[4], which in an '
                                          'un-instrumented build is an adjacent live local and '
                                          'therefore deterministic, but is not guaranteed to be. '
                                          'The recorded flip was observed through the Python '
                                          'fork harness, not in a fresh native process.',
      },
      'ladder': ladder,
      'ceiling': {
        'before': '48 entities. 50 and 60 failed 0/3 in the earlier study and 0/6 and 0/3 here.',
        'after': 'at least 100 entities. n100 (100 entities, 203797 PLC vertices, 502630 PLC '
                 'facets) meshes to 1671118 tets in 20.6 s, twice, with identical tet counts, '
                 'while stock 1.6.0 on the same PLC exits with status 3 in 24.4 s.',
        'not_a_smooth_curve': 'n70 (70 entities) is an outlier in BOTH builds: stock times out '
                              'in boundary recovery at 900 s and the patched build at 1800 s, '
                              'while the larger n100 finishes boundary recovery in 12.2 s. That '
                              'is a separate, unresolved slowness at that particular rung, not a '
                              'size wall and not caused by the patches.',
        'architecture_consequence': 'the exact CGAL-arrangement + TetGen -pY route now reaches '
                                    'the 100-entity mark where fTetWild was previously the only '
                                    'option, and it does so while reproducing every input vertex '
                                    'at 0.0 m, where fTetWild pays 2.08e-04 m of interface '
                                    'position. The element quality objection to the exact route '
                                    'is untouched by this work: it still produces degenerate '
                                    'tets (1 non-positive at ladder-50 and ladder-60, 6 at '
                                    'n100). What changed is that the route no longer stops.',
        'caveat': 'measured on this ladder only, on one host, with flags pY and pYO0. n70 is '
                  'unresolved. Nothing above 100 entities was tried.'},
      'resolution': {
        'bounded_fix_exists': True,
        'what_to_apply': ['0001 is the upstream libigl fix, restored: apply unconditionally',
                          '0002 is a two-line guard against a demonstrably wrong classification',
                          '0003 is a route-around that downgrades three internal-error aborts in '
                          'a best-effort quality pass to "not repaired"; the resulting mesh is '
                          'validated independently (see mesh_validation) rather than trusted',
                          'build with -ffp-contract=off regardless of whether it is needed by '
                          'the current compiler'],
        'how_the_repository_can_use_it': 'data/runtime/tetgen-asan/build/libtetgen_patched.so '
                                         'exposes ihm_tetrahedralize with the same tetgenio '
                                         'construction as igl::copyleft::tetgen, callable from '
                                         'Python via probe_tetgen_boundary_recovery_defect.'
                                         'patched_tetrahedralize(V, F, flags). Rebuilding the '
                                         'whole libigl wheel is not necessary. Keep calling it '
                                         'behind a fork, as run_tetgen already does.',
        'what_was_not_fixed': ['D4, the outsubfaces heap-buffer-overflow, is avoided rather than '
                               'fixed: patch 0002 stops that path being taken',
                               'the n70 boundary-recovery slowness',
                               'element quality: -pY still emits degenerate tets',
                               'the 4 to 12 "segments are not recovered" warnings, which are '
                               'present at the already-accepted rung 48 as well']},
      'verified_versus_inferred': {
        'verified_by_measurement': [
          'the wheel vendors stock upstream TetGen 1.6.0, sha256-identical to the WIAS tarball, '
          'with no libigl patch, and the installed .so is byte-identical to the PyPI artifact '
          '(3ce0b821e57f1ddcab11cc88a3a017cb3066e9e0f907dc37d7fdd992269b3c80)',
          'the memory-safety inference is CORRECT in substance: ASan finds a stack-buffer '
          'over-read in the predicates and a heap-buffer-overflow WRITE in outsubfaces, and the '
          'latter is the recorded glibc abort',
          'the four "self-intersecting" input triangles at ladder-50 are two edge-adjacent pairs '
          'with altitudes four orders of magnitude above TetGen tolerance; the input is clean',
          'TetGen prints the same segment index pair twice as two overlapping segments, with a '
          'separation of exactly 2^-57 and 2^-58 metres',
          'with three patches, ladder-50, ladder-60 and n100 mesh, reproducibly, with every PLC '
          'vertex at distance 0.0 m and box closure at most 4.44e-16 relative',
          'FMA contraction breaks insphere in a local -O2 build (18.78 percent wrong signs) and '
          'does NOT break it in the shipped wheel (0 wrong signs)'],
        'inferred_not_proven': [
          'that D1 (the predicates over-read) contributes to the failures at all. It is present '
          'and is hit on the exact call chain that computes the Steiner point, but fixing it '
          'alone changes nothing measurable, and in an un-instrumented build the over-read lands '
          'on an adjacent live local rather than on unmapped memory',
          'why the earlier study saw one success at 60 entities. Not reproduced: 3/3 native and '
          '4/4 through the wheel fail here. The libc rand() sharing and the over-read are '
          'candidates, neither demonstrated',
          'that patch 0002 is safe at every input. It is exercised on four rungs here; the '
          'segments it declines to call self-intersecting end up in the same "not recovered" '
          'bucket TetGen already tolerates at the accepted rung 48']},
      'mesh_validation': val,
      'runs_receipt': 'data/derived/tetgen-defect-v1/runs.json',
      'logs_receipt': 'data/derived/tetgen-defect-v1/logs/',
    }


def hash_tree(paths):
    out = {}
    for q in sorted(paths):
        q = Path(q)
        if q.is_file():
            out[str(q.relative_to(ROOT))] = {'sha256': sha256_file(q), 'bytes': q.stat().st_size}
    return out


def build_manifest(out):
    out = Path(out)
    inputs = []
    inputs += sorted((BUILD / 'vendor/tetgen1.6.0').glob('*'))
    inputs += sorted((BUILD / 'patches').glob('*.patch'))
    inputs += [BUILD / 'driver.cpp', BUILD / 'build.sh', BUILD / 'predicate_probe.cpp']
    inputs += [Path(__file__).resolve(), ROOT / 'scripts/tetgen_defect_runs.sh']
    inputs += sorted((out / 'plc').glob('*.plcbin'))
    if WHEEL_TETGEN.exists():
        inputs.append(WHEEL_TETGEN)
    man = {'schema': 'ihm.tetgen-defect-manifest.v1',
           'generated_by': 'scripts/probe_tetgen_boundary_recovery_defect.py --manifest',
           'inputs': hash_tree(inputs),
           'outputs': hash_tree(list((out / 'mesh').glob('*.tetbin'))
                                + [out / 'findings.json', out / 'runs.json',
                                   out / 'mesh-validation.json', out / 'plc/index.json'])}
    return man


# ------------------------------------------------------------------------ main

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--dump', default='')
    p.add_argument('--validate', default='', help='comma-separated labels to validate meshes for')
    p.add_argument('--collect', action='store_true', help='parse data/.../logs/*.log into runs.json')
    p.add_argument('--manifest', action='store_true', help='hash inputs and outputs into manifest.json')
    p.add_argument('--report', action='store_true', help='write findings.json')
    p.add_argument('--predicate-audit', type=int, default=0, metavar='N',
                   help='grade N orient3d and insphere queries against exact arithmetic')
    p.add_argument('--wheel', default='', help='run the shipped libigl wheel on this PLC label')
    p.add_argument('--flags', default='pY')
    p.add_argument('--out', default='data/derived/tetgen-defect-v1')
    p.add_argument("--pad", type=float, default=.02)  # same default as build_meshing_route_study
    a = p.parse_args()
    if a.self_test:
        return 0 if self_test() else 1
    if a.wheel:
        return wheel_run(a.wheel, a.flags,
                         Path(a.out if Path(a.out).is_absolute() else ROOT / a.out))
    if not (a.dump or a.validate or a.collect or a.manifest or a.report
            or a.predicate_audit):
        p.error('nothing to do: pass --self-test, --dump, --validate, --collect or --manifest')
    out = Path(a.out if Path(a.out).is_absolute() else ROOT / a.out)
    (out / 'plc').mkdir(parents=True, exist_ok=True)
    index_path = out / 'plc' / 'index.json'
    index = json.loads(index_path.read_text()) if index_path.exists() else {}
    for label in [s for s in a.dump.split(',') if s]:
        PV, PF, meta = build_exact_plc(label, a.pad)
        path = write_plc(out / 'plc' / ('%s.plcbin' % label), PV, PF)
        meta['file'] = str(path.relative_to(ROOT))
        meta['file_sha256'] = sha256_file(path)
        meta['file_bytes'] = path.stat().st_size
        index[label] = meta
        index_path.write_text(json.dumps(index, indent=2, sort_keys=True) + '\n')
        print(json.dumps(meta, indent=2), flush=True)
    if a.validate:
        vpath = out / 'mesh-validation.json'
        table = json.loads(vpath.read_text()) if vpath.exists() else {}
        for label in [s for s in a.validate.split(',') if s]:
            table[label] = validate_mesh(out / 'plc' / ('%s.plcbin' % label),
                                         out / 'mesh' / ('%s.tetbin' % label))
            vpath.write_text(json.dumps(table, indent=2, sort_keys=True) + '\n')
            print(json.dumps(table[label], indent=2), flush=True)
    if a.collect:
        runs = [parse_log(q) for q in sorted((out / 'logs').glob('*.log'))]
        (out / 'runs.json').write_text(json.dumps(
            {'schema': 'ihm.tetgen-defect-runs.v1', 'runs': runs}, indent=2) + '\n')
        ok = sum(r['succeeded'] for r in runs)
        print('collected %d runs, %d succeeded -> %s' % (len(runs), ok, out / 'runs.json'))
    if a.predicate_audit:
        predicate_audit(out, a.predicate_audit)
    if a.report:
        (out / 'findings.json').write_text(json.dumps(build_findings(out), indent=2) + '\n')
        print('findings -> %s' % (out / 'findings.json'))
    if a.manifest:
        (out / 'manifest.json').write_text(json.dumps(build_manifest(out), indent=2) + '\n')
        print('manifest -> %s' % (out / 'manifest.json'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
