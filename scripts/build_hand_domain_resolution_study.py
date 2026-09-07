"""Three tetrahedralisations of the same 109-entity left-hand domain, measured the same way.

  ftetwild  re-runs the tightening that two earlier attempts lost to wall clock under host
            contention: epsr 3e-4 / lr 0.02 first, epsr 1e-4 / lr 0.01 only if that one lands
            inside budget. Same mesher, same soup, same validation as the coarse run.
  exact     the same exact CGAL arrangement PLC through the PATCHED TetGen 1.6.0
            (ihm_tetrahedralize in data/runtime/tetgen-asan/build/libtetgen_patched.so), which
            carries the three patches that lifted the exact route past its recorded 48-entity
            ceiling. -pY preserves input facets, so the interface should land at 0.0 m; the
            standing objection is element quality, so the dihedral distribution is measured over
            EVERY tet, not a sample, and non-positive tets are counted before any rewinding.
  compare   coarse (already built), fine, and exact on accuracy, quality, payload and cost.

Every mesh is validated identically: all tets positive, shared nodes at interfaces, residual
pairwise overlap by exact CGAL boolean, per-entity volume against each input surface's divergence
integral, surface deviation on the undisplaced entities, element quality.

Requires the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_hand_domain_resolution_study.py --self-test
  ... --stage ftetwild --resolution fine
  ... --stage exact
  ... --stage compare --stage manifest

Nothing outside --out is written. scripts/build_hand_reflexive_grip_domain.py,
scripts/probe_tetgen_boundary_recovery_defect.py and everything they import are read, never
modified. TetGen runs from a scratch CWD because it drops tetgen-tmpfile_skipped.* on failure.
"""
from pathlib import Path
import argparse, hashlib, json, os, re, shutil, subprocess, sys, tempfile, threading, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_cross_structure_conflict_repair as R   # noqa: E402  imported, never modified
import build_meshing_route_study as S               # noqa: E402  imported, never modified
import build_hand_reflexive_grip_domain as H        # noqa: E402  imported, never modified
import probe_tetgen_boundary_recovery_defect as P   # noqa: E402  imported, never modified
import igl                                          # noqa: E402
import igl.copyleft.cgal as cgal                    # noqa: E402

HAND = ROOT / 'data/derived/hand-reflexive-grip-domain-v1'
DEFECT = ROOT / 'data/derived/tetgen-defect-v1'
OUT_DEFAULT = ROOT / 'data/derived/hand-domain-resolution-study-v1'
DRIVER = ROOT / 'data/runtime/tetgen-asan/build/tetgen_driver_fixall'
FTETWILD = S.FTETWILD

SETTINGS = {'fine':  {'epsr': 3e-4, 'lr': 0.02},
            'finer': {'epsr': 1e-4, 'lr': 0.01}}
EXACT_PAD_M = 0.005                      # the pad the hand lane's exact control used
VOLUME_ERROR_BUDGET = S.VOLUME_ERROR_BUDGET            # 1e-3, relative
SURFACE_DEVIATION_BUDGET_M = S.SURFACE_DEVIATION_BUDGET_M  # 4.89e-6 m

# Host guard. An unrelated job holds most of the unified pool; pushing this host into swap is a
# worse outcome than an abandoned run, so both bounds are enforced by the sampler and any breach
# kills the mesher and is reported as an abandonment with the bound that was hit.
RSS_CAP_BYTES = 8 * (1 << 30)
MEMAVAIL_FLOOR_BYTES = 6 * (1 << 30)
PAGE = os.sysconf('SC_PAGE_SIZE')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for c in iter(lambda: fh.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def mem_available_bytes():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) * 1024
    return -1


def _proc_table():
    out = {}
    for entry in os.listdir('/proc'):
        if not entry.isdigit():
            continue
        try:
            raw = Path('/proc', entry, 'stat').read_bytes()
        except OSError:
            continue
        close = raw.rfind(b')')
        if close < 0:
            continue
        comm = raw[raw.find(b'(') + 1:close].decode('utf-8', 'replace')
        rest = raw[close + 2:].split()
        try:
            out[int(entry)] = (int(rest[1]), comm, int(rest[21]) * PAGE)   # ppid, comm, rss
        except (IndexError, ValueError):
            continue
    return out


def descendants(root_pid, table):
    kids = {}
    for pid, (ppid, _, _) in table.items():
        kids.setdefault(ppid, []).append(pid)
    seen, stack = [], list(kids.get(root_pid, ()))
    while stack:
        pid = stack.pop()
        seen.append(pid)
        stack.extend(kids.get(pid, ()))
    return seen


class Guard:
    """Samples child RSS and host MemAvailable; kills the tree at either bound."""

    def __init__(self, interval=2.0, rss_cap=RSS_CAP_BYTES, floor=MEMAVAIL_FLOOR_BYTES):
        self.interval, self.rss_cap, self.floor = interval, rss_cap, floor
        self.peak_child_rss = 0
        self.peak_self_rss = 0
        self.min_mem_available = mem_available_bytes()
        self.samples = 0
        self.breach = None
        self.killed = []
        self._stop = threading.Event()
        self._thread = None

    def _tick(self):
        table = _proc_table()
        me = os.getpid()
        kids = descendants(me, table)
        rss = sum(table[p][2] for p in kids if p in table)
        self.peak_child_rss = max(self.peak_child_rss, rss)
        self.peak_self_rss = max(self.peak_self_rss, table.get(me, (0, '', 0))[2])
        avail = mem_available_bytes()
        self.min_mem_available = min(self.min_mem_available, avail)
        self.samples += 1
        if self.breach is None and kids and (rss > self.rss_cap or avail < self.floor):
            self.breach = {'reason': 'child_rss_over_cap' if rss > self.rss_cap
                           else 'host_mem_available_under_floor',
                           'child_rss_bytes': int(rss), 'rss_cap_bytes': int(self.rss_cap),
                           'mem_available_bytes': int(avail), 'floor_bytes': int(self.floor),
                           'at_seconds': time.monotonic() - self.began}
            for pid in kids:
                try:
                    os.kill(pid, 9)
                    self.killed.append({'pid': pid, 'comm': table[pid][1]})
                except OSError:
                    pass

    def _loop(self):
        while not self._stop.wait(self.interval):
            try:
                self._tick()
            except Exception:
                pass

    def __enter__(self):
        self.began = time.monotonic()
        self._tick()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=5)
        try:
            self._tick()
        except Exception:
            pass
        return False

    def report(self):
        return {'samples': self.samples, 'sample_interval_s': self.interval,
                'peak_child_rss_bytes': int(self.peak_child_rss),
                'peak_child_rss_gib': round(self.peak_child_rss / (1 << 30), 3),
                'peak_driver_process_rss_bytes': int(self.peak_self_rss),
                'peak_driver_process_rss_gib': round(self.peak_self_rss / (1 << 30), 3),
                'host_mem_available_min_bytes': int(self.min_mem_available),
                'host_mem_available_min_gib': round(self.min_mem_available / (1 << 30), 3),
                'rss_cap_bytes': int(self.rss_cap), 'mem_available_floor_bytes': int(self.floor),
                'breach': self.breach, 'killed': self.killed,
                'basis': 'RSS summed over every descendant process, sampled from /proc/<pid>/stat; '
                         'host headroom from /proc/meminfo MemAvailable'}


def headroom_or_die(need_gib, note):
    avail = mem_available_bytes()
    if avail < need_gib * (1 << 30):
        raise SystemExit('ABORT before %s: MemAvailable %.1f GiB is under the %.1f GiB this stage '
                         'needs. Nothing was run.' % (note, avail / (1 << 30), need_gib))
    return avail


# ------------------------------------------------------------------ element quality, every tet

def full_quality(TV, TT, chunk=250_000):
    """Dihedral and shape statistics over EVERY tet. The stored coarse numbers came from a
    400 000-tet sample; -pY's known objection is a tail, and a tail is not sampled honestly."""
    faces = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
    pairs = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    n = len(TT)
    mind = np.empty(n)
    maxd = np.empty(n)
    gamma = np.empty(n)
    signed = np.empty(n)
    for s in range(0, n, chunk):
        Q = TV[TT[s:s + chunk]]
        E = Q[:, 1:] - Q[:, :1]
        v = np.linalg.det(np.swapaxes(E, 1, 2)) / 6.0
        signed[s:s + chunk] = v
        N = np.cross(Q[:, faces[:, 1]] - Q[:, faces[:, 0]], Q[:, faces[:, 2]] - Q[:, faces[:, 0]])
        N = N / np.maximum(np.linalg.norm(N, axis=2, keepdims=True), 1e-300)
        ang = np.stack([np.degrees(np.arccos(np.clip(
            -np.einsum('ij,ij->i', N[:, a], N[:, b]), -1.0, 1.0))) for a, b in pairs], 1)
        mind[s:s + chunk] = ang.min(1)
        maxd[s:s + chunk] = ang.max(1)
        L = np.stack([np.linalg.norm(Q[:, b] - Q[:, a], axis=1)
                      for a, b in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]], 1)
        gamma[s:s + chunk] = (12.0 * np.cbrt(3.0 * np.abs(v)) ** 2
                              / np.maximum((L ** 2).sum(1), 1e-300))
    hist, _ = np.histogram(mind, bins=np.arange(0, 181, 5))
    return {'tets': int(n), 'sampled': 'all',
            'min_dihedral_deg': float(np.nanmin(mind)),
            'min_dihedral_p01_deg': float(np.nanpercentile(mind, 1)),
            'min_dihedral_p05_deg': float(np.nanpercentile(mind, 5)),
            'min_dihedral_median_deg': float(np.nanmedian(mind)),
            'max_dihedral_deg': float(np.nanmax(maxd)),
            'max_dihedral_p99_deg': float(np.nanpercentile(maxd, 99)),
            'tets_with_min_dihedral_under_1deg': int((mind < 1).sum()),
            'tets_with_min_dihedral_under_5deg': int((mind < 5).sum()),
            'tets_with_min_dihedral_under_10deg': int((mind < 10).sum()),
            'fraction_under_5deg': float((mind < 5).mean()),
            'fraction_under_10deg': float((mind < 10).mean()),
            'nan_dihedral_tets': int(np.isnan(mind).sum()),
            'shape_gamma_min': float(np.nanmin(gamma)), 'shape_gamma_mean': float(np.nanmean(gamma)),
            'shape_gamma_p01': float(np.nanpercentile(gamma, 1)),
            'signed_volume_min_m3': float(signed.min()), 'signed_volume_max_m3': float(signed.max()),
            'tets_with_nonpositive_signed_volume': int((signed <= 0).sum()),
            'tets_with_negative_signed_volume': int((signed < 0).sum()),
            'tets_with_exactly_zero_signed_volume': int((signed == 0).sum()),
            'min_dihedral_histogram_5deg_bins': [int(x) for x in hist],
            'gamma_note': 'gamma = 12*(3V)^(2/3)/sum(edge^2); 1.0 is the regular tetrahedron',
            'orientation_note': 'signed volumes are in the mesher\'s own emitted vertex order, '
                                'before any rewinding'}


# ------------------------------------------------------------------ shared validation

def validate_domain(TV, TT, parts, priority, workers, box_volume=None, quality_all=True):
    """The coarse run's validation, applied to any tet mesh of the padded box."""
    began = time.monotonic()
    raw_quality = full_quality(TV, TT) if quality_all else None
    mesh, TT, owner, det = S.label_and_measure(TV, TT, parts, priority, box_volume=box_volume)
    if raw_quality is not None:
        mesh['quality_all_tets'] = raw_quality
    out = {'conforming_mesh': mesh,
           'shared_node_audit': S.shared_node_audit(TV, TT, owner, parts)}
    env = {p['entity_id'] for p in parts if p['role'] == H.ENVELOPE_ROLE}
    per = {r['entity_id']: r for r in mesh['per_structure']}
    for r in mesh['per_structure']:
        r['volume_conceded_m3'] = r['claimed_volume_m3'] - r['owned_volume_m3']
        r['fraction_conceded'] = (r['volume_conceded_m3'] / r['claimed_volume_m3']
                                  if r['claimed_volume_m3'] > 0 else None)
    tissue = [r for r in mesh['per_structure'] if r['entity_id'] not in env]
    out['ownership'] = {
        'note': 'the outer body envelope claims every structure tet by construction; the '
                'tissue-level numbers exclude it',
        'tissue_volume_conceded_m3': float(sum(r['volume_conceded_m3'] for r in tissue)),
        'envelope_owned_volume_m3': float(sum(per[e]['owned_volume_m3'] for e in env if e in per)),
        'entities_fully_displaced': [r['entity_id'] for r in tissue if r['owned_tets'] == 0]}
    shells, fid = S.interface_fidelity(TV, TT, owner, parts, det)
    undisplaced = {r['entity_id'] for r in tissue
                   if r['owned_tets'] == r['claimed_tets'] and r['claimed_tets'] > 0}
    clean = [r for r in fid['per_structure'] if r.get('faces') and r['entity_id'] in undisplaced]
    fid['undisplaced_reference'] = {
        'entities': len(clean),
        'basis': 'entities whose owned tet set equals their claimed tet set, so nothing displaced '
                 'them and the deviation is the mesher envelope alone',
        'max_deviation_m': max((r['deviation_max_m'] for r in clean), default=None),
        'median_of_per_entity_max_deviation_m':
            float(np.median([r['deviation_max_m'] for r in clean])) if clean else None,
        'p90_of_per_entity_max_deviation_m':
            float(np.percentile([r['deviation_max_m'] for r in clean], 90)) if clean else None,
        'max_abs_volume_relative_error':
            max((abs(r['volume_relative_error'] or 0.0) for r in clean), default=None),
        'median_abs_volume_relative_error':
            float(np.median([abs(r['volume_relative_error'] or 0.0) for r in clean])) if clean else None,
        'worst': sorted([{k: r[k] for k in ('entity_id', 'name', 'role', 'faces',
                                            'deviation_max_m', 'volume_relative_error')}
                         for r in clean], key=lambda r: -r['deviation_max_m'])[:10]}
    out['interface_fidelity'] = fid
    out['resolved_disjointness'] = H.residual_overlap_parallel(shells, parts, workers)
    e = np.array([abs(r['claimed_relative_volume_error'] or 0.0) for r in mesh['per_structure']])
    out['acceptance'] = {
        'all_tets_positive_volume': bool(mesh['all_positive_volume']),
        'tets_negatively_oriented_by_mesher': mesh['tets_negatively_oriented_by_mesher'],
        'shared_nodes_at_interfaces': bool(out['shared_node_audit']['shared_nodes_at_interfaces']),
        'max_abs_per_structure_volume_error': mesh['max_abs_claimed_volume_error'],
        'volume_error_budget': VOLUME_ERROR_BUDGET,
        'per_structure_volume_error_distribution': {
            'entities': int(len(e)), 'median': float(np.median(e)),
            'p90': float(np.percentile(e, 90)), 'max': float(e.max()),
            'entities_over_1e-3': int((e > 1e-3).sum()), 'entities_over_1e-2': int((e > 1e-2).sum()),
            'entities_over_1e-1': int((e > 1e-1).sum())},
        'volume_budget_met': bool(float(np.median(e)) < VOLUME_ERROR_BUDGET
                                  and int((e > 1e-3).sum()) == 0),
        'surface_deviation_budget_m': SURFACE_DEVIATION_BUDGET_M,
        'max_interface_deviation_m': fid['max_deviation_m'],
        'max_undisplaced_deviation_m': fid['undisplaced_reference']['max_deviation_m'],
        'median_undisplaced_deviation_m': fid['undisplaced_reference']['median_of_per_entity_max_deviation_m'],
        'surface_budget_met': bool((fid['undisplaced_reference']['max_deviation_m'] or 1.0)
                                   < SURFACE_DEVIATION_BUDGET_M),
        'residual_pairwise_overlap_m3': out['resolved_disjointness']['total_residual_overlap_volume_m3'],
        'min_dihedral_deg': (raw_quality or mesh['quality'])['min_dihedral_deg'],
        'tets_with_min_dihedral_under_5deg':
            (raw_quality or mesh['quality'])['tets_with_min_dihedral_under_5deg']}
    out['validation_seconds'] = time.monotonic() - began
    return out, {'TV': TV, 'TT': TT, 'owner': owner, 'det': det}


# ------------------------------------------------------------------ stage: fTetWild

def stage_ftetwild(out, parts, priority, label, args):
    setting = SETTINGS[label]
    workdir = out / 'ftetwild' / label
    headroom_or_die(args.need_gib, 'fTetWild %s (epsr %g / lr %g)'
                    % (label, setting['epsr'], setting['lr']))
    print('ftetwild[%s]: epsr=%g lr=%g timeout=%ds cap=%.1f GiB'
          % (label, setting['epsr'], setting['lr'], args.ftw_timeout,
             RSS_CAP_BYTES / (1 << 30)), flush=True)
    if workdir.exists():
        shutil.rmtree(workdir)
    began = time.monotonic()
    rec = TV = TT = None
    with Guard() as guard:
        try:
            rec, TV, TT = S.run_ftetwild(parts, workdir, setting['epsr'], setting['lr'],
                                         timeout=args.ftw_timeout)
        except BaseException as error:
            rec = {'succeeded': False, 'failed': type(error).__name__,
                   'timeout_s': args.ftw_timeout, 'seconds': time.monotonic() - began,
                   'message': str(error)[:400]}
    report = {'label': label, 'route': 'B-ftetwild', 'mesher': 'fTetWild',
              'mesher_binary': str(FTETWILD.relative_to(ROOT)),
              'epsilon_relative': setting['epsr'], 'edge_length_relative': setting['lr'],
              'member_count': len(parts), 'ftetwild': rec, 'memory': guard.report()}
    if guard.breach:
        report['abandoned'] = {'why': 'memory bound', **guard.breach}
    if not rec or not rec.get('succeeded'):
        report['wall_seconds'] = time.monotonic() - began
        report['acceptance'] = {'meshed': False}
        if not report.get('abandoned') and rec and rec.get('failed') == 'TimeoutExpired':
            report['abandoned'] = {'why': 'wall-clock bound',
                                   'timeout_s': args.ftw_timeout,
                                   'seconds': rec.get('seconds')}
        return report, None
    V, _, _ = S.soup(parts)
    diag = float(np.linalg.norm(V.max(0) - V.min(0)))
    report['envelope'] = {'input_bbox_diagonal_m': diag,
                          'epsilon_absolute_m': setting['epsr'] * diag,
                          'target_edge_length_m': setting['lr'] * diag,
                          'note': 'fTetWild normalises both by the bounding box diagonal'}
    lo, hi = TV.min(0), TV.max(0)
    val, mesh = validate_domain(TV, TT, parts, priority, args.workers,
                                box_volume=float(np.prod(hi - lo)))
    report.update(val)
    report['conforming_mesh']['mesher_box_volume_m3'] = float(np.prod(hi - lo))
    report['wall_seconds'] = time.monotonic() - began
    return report, mesh


# ------------------------------------------------------------------ stage: exact patched TetGen

def build_plc(out, parts, pad, snap_below=0.0):
    """The exact CGAL arrangement plus the padded box, built exactly as the hand lane's exact
    control built it: box corners first, arrangement vertices offset by 8.

    snap_below > 0 selects the repair module's documented short-edge-snap fallback, which merges
    edges below TetGen's own segment tolerance and reverts itself if the merge opens a hole."""
    stem = 'hand-109' if not snap_below else 'hand-109-snap%g' % snap_below
    path = out / 'plc' / ('%s.plcbin' % stem)
    meta_path = out / 'plc' / ('%s.json' % stem)
    if path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if meta.get('pad_m') == pad and meta.get('file_sha256') == sha(path):
            PV, PF = P.read_plc(path)
            return PV, PF, meta, path
    began = time.monotonic()
    W, DF, _, plc = R.arrange(parts, snap_below=snap_below)
    lo, hi = W.min(0) - pad, W.max(0) + pad
    PV = np.ascontiguousarray(np.concatenate([lo + R.BOX_CORNERS * (hi - lo), W]))
    PF = np.ascontiguousarray(np.concatenate([R.BOX, DF + 8]).astype(np.int64))
    P.write_plc(path, PV, PF)
    meta = {'member_count': len(parts), 'pad_m': pad, 'snap_below_m': snap_below,
            'box_volume_m3': float(np.prod(hi - lo)),
            'plc_vertices': int(len(PV)), 'plc_facets': int(len(PF)),
            'arrangement': plc, 'build_seconds': time.monotonic() - began,
            'file': str(path.relative_to(ROOT)), 'file_sha256': sha(path),
            'file_bytes': path.stat().st_size}
    write_json(meta_path, meta)
    return PV, PF, meta, path


def tetgen_worker(plc_path, tetbin_path, flags):
    """--tetgen-worker: one ihm_tetrahedralize call in its own process, run from a scratch CWD.

    A separate process is what makes a TetGen abort a receipt instead of a dead run, and it is
    what lets getrusage report this run's own peak RSS."""
    import resource
    PV, PF = P.read_plc(plc_path)
    cwd = os.getcwd()
    scratch = tempfile.mkdtemp(prefix='hand-exact-')
    began = time.monotonic()
    status, TV, TT = None, None, None
    err = None
    try:
        os.chdir(scratch)
        TV, TT, status = P.patched_tetrahedralize(PV, PF, flags=flags)
    except BaseException as error:
        err = '%s: %s' % (type(error).__name__, error)
    finally:
        os.chdir(cwd)
        leftovers = sorted(p.name for p in Path(scratch).iterdir())
        shutil.rmtree(scratch, ignore_errors=True)
    rec = {'flags': flags, 'status': status, 'seconds': time.monotonic() - began,
           'max_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
           'scratch_leftovers': leftovers, 'exception': err,
           'plc_vertices': int(len(PV)), 'plc_facets': int(len(PF))}
    if TV is not None and TT is not None and len(TT):
        with open(tetbin_path, 'wb') as fh:
            import struct
            fh.write(P.TETMAGIC)
            fh.write(struct.pack('<qq', len(TV), len(TT)))
            fh.write(np.ascontiguousarray(TV, np.float64).tobytes())
            fh.write(np.ascontiguousarray(TT, np.int64).tobytes())
        rec.update(tet_vertices=int(len(TV)), tets=int(len(TT)), wrote=str(tetbin_path))
    rec['succeeded'] = bool(status == 0 and rec.get('tets', 0) > 0)
    print('RESULT %s' % json.dumps(rec, allow_nan=False), flush=True)
    return 0 if rec['succeeded'] else 1


def run_patched(out, plc_path, flags, timeout, stem='hand-109'):
    tetbin = out / 'mesh' / ('%s-%s.tetbin' % (stem, flags))
    tetbin.parent.mkdir(parents=True, exist_ok=True)
    log = out / 'logs' / ('patched-%s-%s.log' % (stem, flags))
    log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(Path(__file__).resolve()), '--tetgen-worker',
           '--plc', str(plc_path), '--tetbin', str(tetbin), '--flags', flags]
    began = time.monotonic()
    scratch = tempfile.mkdtemp(prefix='hand-exact-cwd-')
    with Guard() as guard:
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                                  cwd=scratch)
            rc, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as error:
            rc, stdout, stderr = None, error.stdout or '', error.stderr or ''
            if isinstance(stdout, bytes):
                stdout = stdout.decode('utf-8', 'replace')
            if isinstance(stderr, bytes):
                stderr = stderr.decode('utf-8', 'replace')
    leftovers = sorted(p.name for p in Path(scratch).iterdir())
    shutil.rmtree(scratch, ignore_errors=True)
    log.write_text(stdout + '\n--- stderr ---\n' + stderr)
    rec = {'flags': flags, 'returncode': rc, 'wall_seconds': time.monotonic() - began,
           'library': str(P.PATCHED_SO.relative_to(ROOT)), 'entry_point': 'ihm_tetrahedralize',
           'log': str(log.relative_to(out)), 'scratch_leftovers_from_tetgen': leftovers,
           'stdout_tail': stdout.strip().splitlines()[-20:],
           'stderr_tail': stderr.strip().splitlines()[-20:],
           'memory': guard.report()}
    for line in stdout.splitlines():
        if line.startswith('RESULT '):
            rec['result'] = json.loads(line[7:])
    rec['succeeded'] = bool((rec.get('result') or {}).get('succeeded'))
    if rec['succeeded']:
        rec['tetbin'] = str(tetbin.relative_to(out))
    if rc is not None and rc < 0:
        rec['signal'] = -rc
        rec['note'] = 'the worker was terminated by signal %d; TetGen aborts the process on some ' \
                      'error paths, which is why it runs out of process' % -rc
    return rec, (tetbin if rec['succeeded'] else None)


def run_native_driver(out, plc_path, flags, timeout, stem='hand-109'):
    """Cross-check: the same PLC through the standalone patched driver, no Python in its address
    space. Independent of the ctypes path, so a disagreement would be visible."""
    if not DRIVER.exists():
        return {'skipped': 'no native driver at %s' % DRIVER}
    tetbin = out / 'mesh' / ('%s-native-%s.tetbin' % (stem, flags))
    log = out / 'logs' / ('native-%s-%s.log' % (stem, flags))
    log.parent.mkdir(parents=True, exist_ok=True)
    tetbin.parent.mkdir(parents=True, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix='hand-exact-native-')
    began = time.monotonic()
    with Guard() as guard:
        try:
            proc = subprocess.run([str(DRIVER), str(plc_path), flags, '1', str(tetbin)],
                                  capture_output=True, text=True, timeout=timeout, cwd=scratch)
            rc, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired:
            rc, stdout, stderr = None, '', ''
    leftovers = sorted(p.name for p in Path(scratch).iterdir())
    shutil.rmtree(scratch, ignore_errors=True)
    log.write_text(stdout + '\n--- stderr ---\n' + stderr)
    rec = {'driver': str(DRIVER.relative_to(ROOT)), 'driver_sha256': sha(DRIVER), 'flags': flags,
           'returncode': rc, 'wall_seconds': time.monotonic() - began,
           'scratch_leftovers_from_tetgen': leftovers, 'memory': guard.report(),
           'stdout_tail': stdout.strip().splitlines()[-12:]}
    for line in stdout.splitlines():
        if line.startswith('RESULT '):
            rec['result'] = json.loads(line[7:])
    rec['succeeded'] = bool((rec.get('result') or {}).get('status') == 0
                            and (rec.get('result') or {}).get('tets', -1) > 0)
    if rec['succeeded']:
        rec['tetbin'] = str(tetbin.relative_to(out))
    return rec


def stage_exact(out, parts, priority, args):
    headroom_or_die(args.need_gib, 'the exact patched-TetGen route')
    began = time.monotonic()
    PV, PF, plc_meta, plc_path = build_plc(out, parts, EXACT_PAD_M, args.snap_below)
    print('exact: PLC %d vertices %d facets (valid=%s)'
          % (plc_meta['plc_vertices'], plc_meta['plc_facets'],
             plc_meta['arrangement']['plc_valid']), flush=True)
    report = {'label': args.exact_label, 'route': 'A-exact-arrangement + patched TetGen',
              'snap_below_m': args.snap_below,
              'mesher': 'TetGen 1.6.0 + patches 0001/0002/0003',
              'member_count': len(parts), 'plc': plc_meta,
              'patches': {'source': 'data/runtime/tetgen-asan/patches',
                          'findings': 'data/derived/tetgen-defect-v1/findings.json'},
              'attempts': []}
    chosen = None
    for flags in [f for f in args.flags.split(',') if f]:
        print('exact: ihm_tetrahedralize flags=%s' % flags, flush=True)
        rec, tetbin = run_patched(out, plc_path, flags, args.tetgen_timeout, plc_path.stem)
        report['attempts'].append(rec)
        print('  status=%s tets=%s %.1fs' % ((rec.get('result') or {}).get('status'),
                                             (rec.get('result') or {}).get('tets'),
                                             rec['wall_seconds']), flush=True)
        if rec['succeeded'] and chosen is None:
            chosen = (flags, tetbin, rec)
    if chosen is None:
        report['acceptance'] = {'meshed': False}
        report['wall_seconds'] = time.monotonic() - began
        return report, None
    flags, tetbin, rec = chosen
    report['accepted_flags'] = flags
    report['native_driver_crosscheck'] = run_native_driver(out, plc_path, flags,
                                                           args.tetgen_timeout, plc_path.stem)
    report['exactness'] = P.validate_mesh(plc_path, tetbin)
    TV, TT = P.read_tetbin(tetbin)
    TV = np.ascontiguousarray(TV); TT = np.ascontiguousarray(TT)
    val, mesh = validate_domain(TV, TT, parts, priority, args.workers,
                                box_volume=plc_meta['box_volume_m3'])
    report.update(val)
    report['acceptance']['plc_vertices_reproduced_exactly'] = \
        report['exactness']['plc_vertices_reproduced_exactly']
    report['acceptance']['plc_vertices_absent_from_output'] = \
        report['exactness']['plc_vertices_absent_from_output']
    report['acceptance']['steiner_points'] = report['exactness']['steiner_points']
    report['acceptance']['tets_with_nonpositive_volume_as_emitted'] = \
        report['exactness']['tets_with_nonpositive_volume']
    report['wall_seconds'] = time.monotonic() - began
    return report, mesh


# ------------------------------------------------------------------ stage: diagnose

PATCHED_SRC = ROOT / 'data/runtime/tetgen-asan/build/src-p123'
TERM_RE = 'terminatetetgen'


def stage_diagnose(out, args):
    """Locate the exact source site of a status-2 abort, with a backtrace as the receipt.

    Builds an -O0 -g driver from the SAME patched tree the shared library is built from
    (data/runtime/tetgen-asan/build/src-p123) into this study's own directory, then runs it under
    gdb with `catch throw`. Nothing is written outside --out."""
    diag = out / 'diagnosis'
    diag.mkdir(parents=True, exist_ok=True)
    exe = diag / 'tetgen_driver_p123_dbg'
    rec = {'patched_source_tree': str(PATCHED_SRC.relative_to(ROOT)),
           'patched_source_sha256': {n: sha(PATCHED_SRC / n)
                                     for n in ('tetgen.cxx', 'tetgen.h', 'predicates.cxx')
                                     if (PATCHED_SRC / n).exists()},
           'debug_driver': str(exe.relative_to(out))}
    if not exe.exists():
        cmd = ['g++', '-std=c++11', '-DTETLIBRARY', '-I', str(PATCHED_SRC),
               '-ffp-contract=off', '-O0', '-g', '-o', str(exe),
               str(ROOT / 'data/runtime/tetgen-asan/driver.cpp'),
               str(PATCHED_SRC / 'tetgen.cxx'), str(PATCHED_SRC / 'predicates.cxx')]
        began = time.monotonic()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        rec['build'] = {'returncode': proc.returncode, 'seconds': time.monotonic() - began,
                        'stderr_tail': proc.stderr.strip().splitlines()[-10:]}
        if proc.returncode != 0:
            rec['failed'] = 'the debug driver did not build'
            write_json(out / 'exact-blocker.json', rec)
            return rec
    rec['debug_driver_sha256'] = sha(exe)
    log = diag / ('gdb-%s.log' % args.flags.split(',')[0])
    if not log.exists():
        plc = out / 'plc' / 'hand-109.plcbin'
        scratch = tempfile.mkdtemp(prefix='hand-exact-gdb-')
        began = time.monotonic()
        proc = subprocess.run(['gdb', '-batch', '-ex', 'catch throw', '-ex', 'run',
                               '-ex', 'bt 25', '--args', str(exe), str(plc),
                               args.flags.split(',')[0], '1'],
                              capture_output=True, text=True, timeout=args.tetgen_timeout,
                              cwd=scratch)
        shutil.rmtree(scratch, ignore_errors=True)
        log.write_text(proc.stdout + '\n--- stderr ---\n' + proc.stderr)
        rec['gdb_seconds'] = time.monotonic() - began
    text = log.read_text()
    frames = [l for l in text.splitlines() if re.match(r'^#\d+\s', l)]
    rec['backtrace'] = frames[:12]
    rec['backtrace_log'] = str(log.relative_to(out))
    site = func = None
    for f in frames:
        if 'tetgen.cxx:' in f and TERM_RE not in f:
            site = f.rsplit('tetgen.cxx:', 1)[1].split()[0]
            func = f.split(' in ', 1)[1].split(' (')[0] if ' in ' in f else None
            break
    if site and site.isdigit():
        n = int(site)
        src = (PATCHED_SRC / 'tetgen.cxx').read_text().splitlines()
        rec['abort_site'] = {'file': 'tetgen.cxx', 'line': n, 'function': func,
                             'source': [{'line': i + 1, 'text': src[i]}
                                        for i in range(max(0, n - 14), min(len(src), n + 4))]}
    rec['flag_sweep'] = {'basis': 'every switch string tried on this PLC in this study',
                         'attempts': []}
    doc = json.loads((out / 'study.json').read_text()) if (out / 'study.json').exists() else {}
    for label, run in sorted((doc.get('runs') or {}).items()):
        for att in run.get('attempts', ()):
            rec['flag_sweep']['attempts'].append(
                {'run': label, 'plc_snap_below_m': run.get('snap_below_m'), 'flags': att['flags'],
                 'status': (att.get('result') or {}).get('status'),
                 'tets': (att.get('result') or {}).get('tets'),
                 'seconds': (att.get('result') or {}).get('seconds')})
    write_json(out / 'exact-blocker.json', rec)
    return rec


# ------------------------------------------------------------------ stage: compare

def coarse_reference():
    doc = json.loads((HAND / 'domain.json').read_text())
    return doc['resolutions']['coarse']


def row_for(label, rep, npz):
    if not rep or not rep.get('conforming_mesh'):
        return {'label': label, 'meshed': False,
                'abandoned': rep.get('abandoned') if rep else None,
                'wall_seconds': (rep or {}).get('wall_seconds')}
    m = rep['conforming_mesh']
    a = rep.get('acceptance', {})
    q = m.get('quality_all_tets') or m['quality']
    fid = (rep.get('interface_fidelity') or {}).get('undisplaced_reference') or {}
    dist = a.get('per_structure_volume_error_distribution', {})
    mem = ((rep.get('memory') or {}) if 'memory' in rep else
           ((rep.get('attempts') or [{}])[-1].get('memory') or {}))
    return {
        'label': label,
        'route': rep.get('route'),
        'meshed': True,
        'tets': m['tets'], 'tet_vertices': m['tet_vertices'],
        'payload_npz_bytes': int(npz.stat().st_size) if npz and npz.exists() else None,
        'wall_seconds': rep.get('wall_seconds'),
        'mesher_seconds': ((rep.get('ftetwild') or {}).get('seconds')
                           or ((rep.get('attempts') or [{}])[-1].get('result') or {}).get('seconds')),
        'peak_child_rss_gib': mem.get('peak_child_rss_gib'),
        'all_tets_positive_volume': a.get('all_tets_positive_volume'),
        'tets_nonpositive_as_emitted': a.get('tets_with_nonpositive_volume_as_emitted',
                                             q.get('tets_with_nonpositive_signed_volume')),
        'shared_nodes_at_interfaces': a.get('shared_nodes_at_interfaces'),
        'residual_pairwise_overlap_m3': a.get('residual_pairwise_overlap_m3'),
        'volume_error_median': dist.get('median'), 'volume_error_p90': dist.get('p90'),
        'volume_error_max': dist.get('max'),
        'entities_over_1e-3': dist.get('entities_over_1e-3'),
        'volume_budget_met': a.get('volume_budget_met'),
        'undisplaced_entities': fid.get('entities'),
        'undisplaced_deviation_median_m': fid.get('median_of_per_entity_max_deviation_m'),
        'undisplaced_deviation_max_m': fid.get('max_deviation_m'),
        'max_interface_deviation_m': a.get('max_interface_deviation_m'),
        'surface_budget_met': a.get('surface_budget_met'),
        'plc_vertices_reproduced_exactly': a.get('plc_vertices_reproduced_exactly'),
        'min_dihedral_deg': q.get('min_dihedral_deg'),
        'min_dihedral_median_deg': q.get('min_dihedral_median_deg'),
        'tets_under_5deg': q.get('tets_with_min_dihedral_under_5deg'),
        'fraction_under_5deg': q.get('fraction_under_5deg'),
        'shape_gamma_min': q.get('shape_gamma_min'),
        'shape_gamma_mean': q.get('shape_gamma_mean')}


def localise_defects(out, label):
    """When a hard criterion fails, say WHERE. A distributed quality collapse and a single bad
    sliver in one entity are different objections and must not be reported as the same number."""
    path = out / ('tetmesh-%s.npz' % label)
    if not path.exists():
        return None
    names = {}
    for i, line in enumerate((HAND / 'members.jsonl').read_text().splitlines()):
        r = json.loads(line)
        names[i] = {'entity_id': r['entity_id'], 'name': r['name'], 'role': r['role']}
    d = np.load(path, allow_pickle=False)
    TV, TT, owner = d['TV'], d['TT'], d['owner']
    det = np.linalg.det(np.swapaxes(TV[TT[:, 1:]] - TV[TT[:, 0, None]], 1, 2)) / 6.0
    bad = np.where(det <= 0)[0]
    u, inv, c = np.unique(TV, axis=0, return_inverse=True, return_counts=True)
    dup = np.where(c > 1)[0]
    dup_rows = []
    for k in dup:
        ids = np.where(inv == k)[0]
        owners = set()
        incident = []
        for vid in ids:
            t = np.where((TT == vid).any(1))[0]
            incident.append(int(len(t)))
            owners |= set(owner[t].tolist())
        dup_rows.append({'position_m': [float(x) for x in u[k]],
                         'vertex_ids': [int(x) for x in ids],
                         'incident_tets_per_vertex': incident,
                         'owners': sorted(names.get(o, {'entity_id': 'complement'})['entity_id']
                                          for o in owners if o >= 0) or ['complement']})
    return {
        'mesh': path.name,
        'nonpositive_tets': {
            'count': int(len(bad)),
            'signed_volumes_m3': [float(x) for x in det[bad][:20]],
            'owners': sorted({names.get(int(o), {'entity_id': 'complement'})['entity_id']
                              for o in owner[bad]}),
            'owner_names': sorted({names.get(int(o), {'name': 'interstitial complement'})['name']
                                   for o in owner[bad]})},
        'coincident_vertex_positions': {
            'count': int(len(dup)), 'detail': dup_rows[:10],
            'note': 'two vertex ids at one position split the material there into two node sets; '
                    'that is a crack, and it is why shared_nodes_at_interfaces is false'},
        'basis': 'recomputed from the saved tet mesh, not from the run log'}


def stage_compare(out):
    doc = json.loads((out / 'study.json').read_text()) if (out / 'study.json').exists() else {}
    reps = doc.get('runs', {})
    rows = [row_for('coarse-ftetwild', coarse_reference(), HAND / 'tetmesh-coarse.npz')]
    for label in ('fine', 'finer', 'exact', 'exact-snapped', 'exact-flagsweep'):
        if label in reps:
            rows.append(row_for(label, reps[label], out / ('tetmesh-%s.npz' % label)))
    meshed = [r for r in rows if r.get('meshed')]
    for r in meshed:
        if r['label'] != 'coarse-ftetwild' and not (
                r.get('all_tets_positive_volume') and r.get('shared_nodes_at_interfaces')):
            r['defect_localisation'] = localise_defects(out, r['label'])
    clean = [r for r in meshed
             if r.get('all_tets_positive_volume') and r.get('shared_nodes_at_interfaces')
             and r.get('residual_pairwise_overlap_m3') == 0.0]
    by_dev = sorted([r for r in meshed if r.get('undisplaced_deviation_median_m') is not None],
                    key=lambda r: r['undisplaced_deviation_median_m'])
    pick = None
    if clean:
        pick = min(clean, key=lambda r: (r.get('volume_error_median') or 1.0))
    verdict = {
        'meshes_built': [r['label'] for r in meshed],
        'meshes_passing_every_hard_criterion': [r['label'] for r in clean],
        'hard_criteria': 'all tets positive volume, shared nodes at interfaces (one vertex array, '
                         'every interface face carried by exactly two tets, no coincident '
                         'duplicate vertex), residual pairwise overlap exactly 0.0 m3',
        'no_mesh_meets_the_accuracy_budgets': all(not r.get('volume_budget_met')
                                                  and not r.get('surface_budget_met')
                                                  for r in meshed),
        'best_accuracy': by_dev[0]['label'] if by_dev else None,
        'recommended_canonical_domain': pick['label'] if pick else None,
        'recommendation_rule': 'among the meshes that pass every hard criterion, the one with the '
                               'lowest median per-entity volume error. Accuracy alone does not '
                               'decide it: a domain with a zero-volume element or two coincident '
                               'nodes cannot carry a soft-body solve at all, so the hard criteria '
                               'gate and accuracy ranks within the gate.'}
    cmp_doc = {'schema': 'ihm.hand-domain-resolution-comparison.v1',
               'budgets': {'volume_relative': VOLUME_ERROR_BUDGET,
                           'surface_deviation_m': SURFACE_DEVIATION_BUDGET_M,
                           'budget_source': 'build_meshing_route_study.py'},
               'rows': rows, 'verdict': verdict,
               'exact_route_blocker': (json.loads((out / 'exact-blocker.json').read_text())
                                       if (out / 'exact-blocker.json').exists() else None)}
    if cmp_doc['exact_route_blocker']:
        b = cmp_doc['exact_route_blocker']
        cmp_doc['exact_route_blocker'] = {
            'abort_site': b.get('abort_site'), 'backtrace': b.get('backtrace'),
            'flag_sweep': b.get('flag_sweep'),
            'patched_source_sha256': b.get('patched_source_sha256')}
    write_json(out / 'comparison.json', cmp_doc)
    return cmp_doc


# ------------------------------------------------------------------ manifest

def manifest(out):
    inputs = {}
    for rel in ('data/derived/hand-reflexive-grip-domain-v1/manifest.json',
                'data/derived/hand-reflexive-grip-domain-v1/domain.json',
                'data/derived/hand-reflexive-grip-domain-v1/members.jsonl',
                'data/derived/hand-reflexive-grip-domain-v1/clipped-members.npz',
                'data/derived/hand-reflexive-grip-domain-v1/tetmesh-coarse.npz',
                'data/derived/tetgen-defect-v1/findings.json',
                'data/derived/tetgen-defect-v1/manifest.json',
                'data/runtime/tetgen-asan/build/libtetgen_patched.so',
                'data/runtime/tetgen-asan/build/tetgen_driver_fixall',
                'data/runtime/tetgen-asan/receipt.json',
                'data/runtime/tolerant-mesher/build/FloatTetwild_bin',
                'scripts/build_hand_reflexive_grip_domain.py',
                'scripts/probe_tetgen_boundary_recovery_defect.py',
                'scripts/build_meshing_route_study.py',
                'scripts/build_cross_structure_conflict_repair.py',
                'scripts/build_hand_domain_resolution_study.py'):
        p = ROOT / rel
        inputs[rel] = sha(p) if p.exists() else None
    arts = {}
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name != 'manifest.json':
            arts[str(p.relative_to(out))] = {'sha256': sha(p), 'bytes': p.stat().st_size}
    doc = {'schema': 'ihm.hand-domain-resolution-study-manifest.v1',
           'inputs_sha256': inputs, 'artifacts': arts,
           'existing_repo_files_modified': False,
           'writes_confined_to': str(out.relative_to(ROOT)),
           'python': sys.version,
           'host': {'cpus': os.cpu_count(),
                    'mem_available_at_manifest_bytes': mem_available_bytes()}}
    write_json(out / 'manifest.json', doc)
    return doc


# ------------------------------------------------------------------ self-test

def self_test():
    ok = True

    def check(name, cond, detail=''):
        nonlocal ok
        ok = ok and bool(cond)
        print('%-46s %s %s' % (name, 'ok' if cond else 'FAIL', detail), flush=True)

    check('patched .so present', P.PATCHED_SO.exists(), str(P.PATCHED_SO))
    check('fTetWild present', FTETWILD.exists(), str(FTETWILD))
    check('hand domain present', (HAND / 'members.jsonl').exists(), str(HAND))
    check('budgets from route study',
          VOLUME_ERROR_BUDGET == 1e-3 and abs(SURFACE_DEVIATION_BUDGET_M - 4.89e-6) < 1e-12,
          '%g / %g m' % (VOLUME_ERROR_BUDGET, SURFACE_DEVIATION_BUDGET_M))

    scratch = Path(tempfile.mkdtemp(prefix='hand-study-selftest-'))
    try:
        # PLC round trip
        lo, hi = np.zeros(3), np.ones(3)
        BV = lo + R.BOX_CORNERS * (hi - lo)
        inner = np.array([0.3, 0.3, 0.3]) + R.BOX_CORNERS * 0.4
        PV = np.ascontiguousarray(np.concatenate([BV, inner]))
        PF = np.ascontiguousarray(np.concatenate([R.BOX, R.BOX + 8]).astype(np.int64))
        path = P.write_plc(scratch / 'cubes.plcbin', PV, PF)
        RV, RF = P.read_plc(path)
        check('plc round trip', np.array_equal(RV, PV) and np.array_equal(RF, PF))

        # ihm_tetrahedralize through the worker, out of process, from a scratch cwd
        tb = scratch / 'cubes.tetbin'
        args = argparse.Namespace()
        rec, produced = run_patched(scratch, path, 'pY', 120)
        check('ihm_tetrahedralize meshes nested cubes', rec['succeeded'],
              'status=%s tets=%s' % ((rec.get('result') or {}).get('status'),
                                     (rec.get('result') or {}).get('tets')))
        check('tetgen left no file in the run cwd', rec['scratch_leftovers_from_tetgen'] == [],
              str(rec['scratch_leftovers_from_tetgen']))
        if rec['succeeded']:
            v = P.validate_mesh(path, scratch / 'mesh' / 'hand-109-pY.tetbin')
            check('every plc vertex reproduced at 0.0 m', v['plc_vertices_reproduced_exactly'])
            check('box closure exact', abs(v['box_volume_closure_relative_error']) < 1e-12,
                  '%.3e' % v['box_volume_closure_relative_error'])
            TV, TT = P.read_tetbin(scratch / 'mesh' / 'hand-109-pY.tetbin')
            q = full_quality(np.ascontiguousarray(TV), np.ascontiguousarray(TT), chunk=64)
            check('full_quality counts every tet', q['tets'] == len(TT) and q['sampled'] == 'all')
            sq = S.tet_quality(np.ascontiguousarray(TV), np.ascontiguousarray(TT))
            check('full_quality agrees with route study on min dihedral',
                  abs(q['min_dihedral_deg'] - sq['min_dihedral_deg']) < 1e-9,
                  '%.6f vs %.6f' % (q['min_dihedral_deg'], sq['min_dihedral_deg']))
            # a deliberately inverted tet must be counted, not hidden
            TTb = TT.copy(); TTb[0] = TTb[0][[0, 2, 1, 3]]
            qb = full_quality(np.ascontiguousarray(TV), np.ascontiguousarray(TTb), chunk=64)
            check('inverted tet is counted as emitted',
                  qb['tets_with_negative_signed_volume'] >= 1,
                  '%d' % qb['tets_with_negative_signed_volume'])

        # validate_domain on a one-part cube domain: exact volume, no overlap, closed shell
        cube = {'entity_id': 'selftest-cube', 'name': 'cube', 'role': 'rigid_bone',
                'system': 'test', 'frame': 'test', 'source_path': '', 'source_sha256': '',
                'V': np.ascontiguousarray(inner), 'F': np.ascontiguousarray(R.BOX.copy()),
                'surface_volume_m3': abs(R.divergence(inner, R.BOX))}
        TV, TT = P.read_tetbin(scratch / 'mesh' / 'hand-109-pY.tetbin')
        val, _ = validate_domain(np.ascontiguousarray(TV), np.ascontiguousarray(TT),
                                 [cube], {'selftest-cube': 0}, 2, box_volume=1.0)
        check('validate_domain recovers the cube volume',
              abs(val['conforming_mesh']['max_abs_claimed_volume_error']) < 1e-9,
              '%.3e' % val['conforming_mesh']['max_abs_claimed_volume_error'])
        check('validate_domain reports zero residual overlap',
              val['resolved_disjointness']['total_residual_overlap_volume_m3'] == 0.0)
        check('validate_domain deviation is 0 m on a preserved facet set',
              val['interface_fidelity']['max_deviation_m'] < 1e-15,
              '%.3e m' % val['interface_fidelity']['max_deviation_m'])

        # guard: reports and does not fire without children
        with Guard(interval=0.05) as g:
            time.sleep(0.2)
        r = g.report()
        check('memory guard samples', r['samples'] >= 2 and r['breach'] is None,
              '%d samples, %.1f GiB free' % (r['samples'], r['host_mem_available_min_gib']))
        # guard: fires and kills at an impossible cap
        with Guard(interval=0.05, rss_cap=1) as g2:
            subprocess.run([sys.executable, '-c', 'import time; time.sleep(5)'])
        check('memory guard kills at the cap', g2.report()['breach'] is not None,
              str((g2.report().get('breach') or {}).get('reason')))

        check('headroom_or_die passes with a tiny need', headroom_or_die(0.001, 'selftest') > 0)
        try:
            headroom_or_die(1e6, 'selftest')
            check('headroom_or_die aborts on an impossible need', False)
        except SystemExit:
            check('headroom_or_die aborts on an impossible need', True)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    leaked = sorted(p.name for p in ROOT.glob('tetgen-tmpfile*'))
    check('no tetgen tmpfile at the repo root', leaked == [], str(leaked))
    print('self-test %s' % ('PASSED' if ok else 'FAILED'), flush=True)
    return ok


# ------------------------------------------------------------------ driver

def save_mesh(out, label, mesh):
    if mesh is None:
        return None
    path = out / ('tetmesh-%s.npz' % label)
    np.savez_compressed(path, TV=mesh['TV'], TT=mesh['TT'], owner=mesh['owner'])
    return path


def run(args):
    out = Path(args.out if Path(args.out).is_absolute() else ROOT / args.out)
    out.mkdir(parents=True, exist_ok=True)
    doc = json.loads((out / 'study.json').read_text()) if (out / 'study.json').exists() else {
        'schema': 'ihm.hand-domain-resolution-study.v1', 'runs': {}}
    doc.setdefault('runs', {})

    def flush():
        # Stages run concurrently in separate processes, so the file on disk is re-read and this
        # process's runs are merged into it. Writing back the copy loaded at startup would
        # silently drop whatever another stage recorded in the meantime.
        cur = json.loads((out / 'study.json').read_text()) if (out / 'study.json').exists() else {}
        merged = dict(cur)
        merged.update({k: v for k, v in doc.items() if k != 'runs'})
        runs = dict(cur.get('runs') or {})
        runs.update(doc.get('runs') or {})
        merged['runs'] = runs
        merged['budgets'] = {'volume_relative': VOLUME_ERROR_BUDGET,
                             'surface_deviation_m': SURFACE_DEVIATION_BUDGET_M}
        write_json(out / 'study.json', merged)

    stages = [s for s in args.stage.split(',') if s]
    parts = priority = None
    if any(s in ('ftetwild', 'exact') for s in stages):
        parts = H.load_parts(HAND)
        priority = H.priority_table(parts)
        print('loaded %d members from %s' % (len(parts), HAND.relative_to(ROOT)), flush=True)
    for stage in stages:
        if stage == 'ftetwild':
            for label in [l for l in args.resolution.split(',') if l]:
                rep, mesh = stage_ftetwild(out, parts, priority, label, args)
                doc['runs'][label] = rep
                p = save_mesh(out, label, mesh)
                if p:
                    rep['tetmesh_npz'] = {'file': p.name, 'bytes': p.stat().st_size}
                flush()
                a = rep.get('acceptance', {})
                print('  %s: meshed=%s tets=%s median_vol_err=%s median_dev=%s m'
                      % (label, a.get('all_tets_positive_volume') is not None,
                         (rep.get('conforming_mesh') or {}).get('tets'),
                         (a.get('per_structure_volume_error_distribution') or {}).get('median'),
                         a.get('median_undisplaced_deviation_m')), flush=True)
                if label == 'fine' and not a.get('volume_budget_met'):
                    print('  fine did not meet the volume budget; finer is only worth running if '
                          'asked for explicitly', flush=True)
        elif stage == 'exact':
            rep, mesh = stage_exact(out, parts, priority, args)
            doc['runs'][args.exact_label] = rep
            p = save_mesh(out, args.exact_label, mesh)
            if p:
                rep['tetmesh_npz'] = {'file': p.name, 'bytes': p.stat().st_size}
            flush()
        elif stage == 'diagnose':
            doc['exact_blocker'] = stage_diagnose(out, args)
            flush()
        elif stage == 'compare':
            doc['comparison'] = stage_compare(out)
            flush()
        elif stage == 'manifest':
            manifest(out)
        else:
            raise SystemExit('unknown stage %r' % stage)
    leaked = sorted(p.name for p in ROOT.glob('tetgen-tmpfile*'))
    if leaked:
        for name in leaked:
            (ROOT / name).unlink()
        print('cleaned %d tetgen tmpfile(s) from the repo root: %s' % (len(leaked), leaked),
              flush=True)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--out', default=str(OUT_DEFAULT.relative_to(ROOT)))
    p.add_argument('--stage', default='')
    p.add_argument('--resolution', default='fine', help='fine, finer, or fine,finer')
    p.add_argument('--flags', default='pY,pYO0', help='TetGen switch strings, tried in order')
    p.add_argument('--snap-below', type=float, default=0.0,
                   help='short-edge snap threshold for the arrangement; 0 is the exact PLC')
    p.add_argument('--exact-label', default='exact')
    p.add_argument('--workers', type=int, default=10)
    p.add_argument('--ftw-timeout', type=int, default=10800)
    p.add_argument('--tetgen-timeout', type=int, default=3600)
    p.add_argument('--need-gib', type=float, default=10.0,
                   help='MemAvailable required before a heavy stage starts')
    p.add_argument('--tetgen-worker', action='store_true', help=argparse.SUPPRESS)
    p.add_argument('--plc', default='')
    p.add_argument('--tetbin', default='')
    a = p.parse_args()
    if a.tetgen_worker:
        return tetgen_worker(a.plc, a.tetbin, a.flags)
    if a.self_test:
        return 0 if self_test() else 1
    if not a.stage:
        p.error('nothing to do: pass --self-test or --stage ftetwild,exact,compare,manifest')
    return run(a)


if __name__ == '__main__':
    sys.exit(main())
