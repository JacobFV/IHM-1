"""Persistent native physiology with explicit units, ownership, and action receipts.

Missing native quantities remain None. Compartments may overlap and must never
be summed as independent stores. A session is single-controller/single-thread.
"""
from dataclasses import dataclass, asdict
import json
import math
import os
from pathlib import Path
import queue
import re
import subprocess
import threading
import shutil
import xml.etree.ElementTree as ET
from . import BASE, RUNTIME, SOURCE_REVISION, available_patients, number, _sha

VARIANTS = ('upstream', 'saturation_bounds', 'saturation_bounds_heatflux',
            'saturation_bounds_heatflux_thermal_units', 'whole_body_integrity', 'whole_body_integrity_renal',
            'whole_body_integrity_gi_water', 'whole_body_integrity_energy', 'whole_body_integrity_depletion')


def runtime_patient(patient):
    return RUNTIME/'biogears-build/runtime/patients'/f'{patient}.xml'


def _ticks(seconds, maximum, name):
    number(seconds, .02, maximum, name)
    ticks = round(seconds*50)
    if abs(ticks-seconds*50) > 1e-7:
        raise ValueError(f'{name} must align to the native 0.02 s step')
    return ticks


@dataclass(frozen=True)
class SessionConfig:
    patient: str = 'StandardMale'
    state_path: str | Path | None = None
    engine_variant: str = 'upstream'
    horizon_s: float = 86400

    def __post_init__(self):
        _ticks(self.horizon_s, 86400, 'horizon_s')
        if not isinstance(self.patient, str) or self.patient not in available_patients():
            raise ValueError('Unknown native patient')
        if self.engine_variant not in VARIANTS:
            raise ValueError('Unknown engine variant')
        if self.state_path is not None:
            if not isinstance(self.state_path, (str, Path)) or not Path(self.state_path).is_file():
                raise ValueError('Native state file unavailable')
            object.__setattr__(self, 'state_path', Path(self.state_path).resolve())


@dataclass(frozen=True)
class Meal:
    name: str = 'mixed_meal'
    carbohydrate_g: float = 0
    protein_g: float = 0
    fat_g: float = 0
    sodium_g: float = 0
    calcium_mg: float = 0
    water_ml: float = 0

    def __post_init__(self):
        if not isinstance(self.name, str) or not re.fullmatch(r'[A-Za-z0-9_]{1,64}', self.name):
            raise ValueError('Meal name must be 1–64 letters, digits, or underscores')
        for key, value in asdict(self).items():
            if key != 'name':
                number(value, 0, 10000, key)
        if not any(v > 0 for k, v in asdict(self).items() if k != 'name'):
            raise ValueError('Meal must contain a positive amount')


class NativeSession:
    """Own one native process; commands act on its continuing state, without replay."""

    def __init__(self, config, output_dir, timeout_s=1800):
        if not isinstance(config, SessionConfig):
            raise ValueError('Expected SessionConfig')
        number(timeout_s, 1, 3600, 'timeout_s')
        self.config, self.timeout_s = config, timeout_s
        self.out = Path(output_dir).resolve()
        executable = RUNTIME/'native_biogears_stream'
        if not executable.is_file():
            raise RuntimeError('Build scripts/build_native_adapter.py first')
        compiled=json.loads(executable.with_suffix('.manifest.json').read_text())
        if _sha(executable)!=compiled['executable_sha256'] or any(_sha(BASE/p)!=sha for p,sha in compiled['source_sha256'].items()):
            raise ValueError('Native stream source/build identity changed; rebuild adapter')
        environment = os.environ.copy()
        variant_manifest = None
        if config.engine_variant != 'upstream':
            variant = RUNTIME/'variants'/config.engine_variant
            variant_manifest = json.loads((variant/'manifest.json').read_text())
            if _sha(variant/'libbiogears.so.8.0.0') != variant_manifest['library_sha256']:
                raise ValueError('Variant library integrity mismatch')
            environment['LD_LIBRARY_PATH'] = str(variant)+os.pathsep+environment.get('LD_LIBRARY_PATH', '')
        library = (RUNTIME/'biogears-build/outputs/Release/lib' if config.engine_variant == 'upstream'
                   else RUNTIME/'variants'/config.engine_variant)/'libbiogears.so.8.0.0'
        linkage=subprocess.check_output(['ldd',str(executable)],env=environment,text=True)
        dependencies={}
        for line in linkage.splitlines():
            if 'not found' in line:
                raise RuntimeError('Native dependency unresolved: '+line)
            if '=>' in line:
                resolved=Path(line.split('=>',1)[1].split(' (',1)[0].strip())
                if resolved.is_file():dependencies[str(resolved.resolve())]=_sha(resolved)
        if str(library.resolve()) not in dependencies:
            raise RuntimeError('Native stream resolves the wrong physiology library')
        patient_input=Path(config.state_path) if config.state_path else runtime_patient(config.patient)
        patient_input_sha256=_sha(patient_input)
        patient_tree=ET.parse(patient_input).getroot()
        if config.state_path:
            patient_tree=next((e for e in patient_tree if e.tag.split('}')[-1]=='Patient'),None)
            if patient_tree is None:raise ValueError('Saved state has no patient identity')
        patient_data={e.tag.split('}')[-1]:(dict(e.attrib) if e.attrib else e.text) for e in patient_tree}
        if self.out.exists() and any(self.out.iterdir()):
            raise ValueError('Native session requires a fresh output directory')
        self.out.mkdir(parents=True, exist_ok=True)
        runtime = RUNTIME/'biogears-build/runtime'
        resource_names=('patients', 'substances', 'environments', 'nutrition', 'config',
                        'ecg', 'xsd', 'UCEDefs.conf', 'BioGearsConfiguration.xml')
        for name in resource_names:
            (self.out/name).symlink_to(runtime/name, target_is_directory=(runtime/name).is_dir())
        (self.out/'states').mkdir()
        state_input=None
        if config.state_path:
            state_input=self.out/'input-state.xml'
            shutil.copyfile(patient_input,state_input)
            if _sha(state_input)!=patient_input_sha256:
                raise ValueError('Initial state changed during input capture')
        self._sequence = 0
        self._elapsed_ticks = 0
        self._pending_meal = False
        self._owner = threading.get_ident()
        self._closed = False
        self._queue = queue.Queue()
        self._journal = (self.out/'receipts.jsonl').open('x')
        command = [str(executable), config.patient, str(state_input or '-'), str(round(config.horizon_s*50))]
        manifest = dict(schema='ihm.native-session.v1', configuration=asdict(config),
                        command=command, source_revision=SOURCE_REVISION,
                        executable_sha256=_sha(executable), variant_manifest=variant_manifest,
                        adapter_sha256={p:_sha(BASE/'scripts'/p) for p in ('native_biogears_stream.cpp','native_body_ports.h')},
                        state_sha256=patient_input_sha256 if config.state_path else None,
                        patient_identity=patient_data, patient_identity_input=str(patient_input),
                        patient_identity_input_sha256=patient_input_sha256,
                        initialization_patient_used=config.state_path is None,
                        build_manifest=compiled, dependency_sha256=dependencies,
                        native_step_s=.02, owner='BioGears', experimental_validation=False,
                        checkpoint_exactness='not established; native serializer may omit controller latches')
        manifest['library_sha256'] = _sha(library)
        manifest['library_path'] = str(library)
        (self.out/'manifest.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')
        try:
            from ihm.assembly.native_environment_evidence import freeze_native_environment,materialize_native_resources
            freeze_native_environment(BASE,self.out,before_start=True)
            detached=materialize_native_resources(BASE,self.out)
            for name in resource_names:
                link=self.out/name
                if not link.is_symlink() or link.resolve()!=(runtime/name).resolve():
                    raise ValueError('Native input root changed before resource selection')
                link.unlink()
                link.symlink_to(detached/name,target_is_directory=(detached/name).is_dir())
            selected_patient=state_input if state_input else detached/'patients'/f'{config.patient}.xml'
            if _sha(selected_patient)!=patient_input_sha256:
                raise ValueError('Detached patient identity differs from the startup manifest')
            (self.out/'execution-inputs.json').write_text(json.dumps(dict(
                native_manifest_sha256=_sha(self.out/'manifest.json'),
                environment_manifest_sha256=_sha(self.out/'environment-inputs/manifest.json'),
                resource_materialization_sha256=_sha(detached/'materialization.json'),
                selected_resource_tree=str(detached),selected_patient_input=str(selected_patient),
                selected_patient_input_sha256=patient_input_sha256,
                selection='Detached resource tree and copied initial state selected before Popen; libraries remain at manifest-resolved paths'),indent=2)+'\n')
            self.process = subprocess.Popen(command, cwd=self.out, env=environment, stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
            self._reader = threading.Thread(target=self._read, daemon=True)
            self._reader.start()
            self._last = self._receive('ready', 0)
            self._pending_meal=self._last['pending_meal']
        except BaseException:
            self.close(graceful=False)
            raise

    def _read(self):
        try:
            with (self.out/'runner_stdout.log').open('x') as log:
                for line in self.process.stdout:
                    log.write(line)
                    log.flush()
                    if line.startswith('IHM\t'):
                        self._queue.put(json.loads(line[4:]))
                    elif re.search(r'FATAL|\[ERROR\]|<ERROR>|entered irreversible state',line,re.I):
                        self._queue.put(RuntimeError('Native engine reported an error; inspect preserved logs'))
        except BaseException as exc:
            self._queue.put(exc)
        finally:
            self._queue.put(EOFError('Native session exited; inspect preserved logs'))

    def _receive(self, status, expected_ticks):
        try:
            value = self._queue.get(timeout=self.timeout_s)
        except queue.Empty:
            raise TimeoutError('Native command timed out; inspect preserved logs') from None
        if isinstance(value, BaseException):
            raise value
        if value.get('sequence') != self._sequence or value.get('status') != status:
            raise RuntimeError(f'Native command rejected or out of sequence: {value}')
        for field in ('time_s', 'elapsed_s', 'origin_s'):
            if not isinstance(value.get(field), (int, float)) or not math.isfinite(value[field]):
                raise RuntimeError('Malformed native clock')
        if abs(value['elapsed_s']-expected_ticks*.02) > 1e-8 or abs(value['time_s']-value['origin_s']-expected_ticks*.02) > 1e-4:
            raise RuntimeError('Native clock differs from acknowledged solver steps')
        if not isinstance(value.get('values'), dict) or not value['values']:
            raise RuntimeError('Missing native ports')
        if not isinstance(value.get('pending_meal'),bool):
            raise RuntimeError('Missing native action state')
        for v in value['values'].values():
            if v is not None and (isinstance(v, bool) or not isinstance(v,(int,float)) or not math.isfinite(v)):
                raise RuntimeError('Invalid native port quantity')
        return value

    def _command(self, op, *args, advance_ticks=0, status='ok'):
        if threading.get_ident() != self._owner:
            raise RuntimeError('NativeSession requires a single controller thread')
        if self._closed:
            raise RuntimeError('Native session is closed')
        self._sequence += 1
        line = ' '.join(map(str, (self._sequence, op, *args)))
        # Preserve the attempted command before sending it, including failed commands.
        self._journal.write(json.dumps({'command': line, 'before_ticks': self._elapsed_ticks})+'\n')
        self._journal.flush()
        try:
            self.process.stdin.write(line+'\n')
            self.process.stdin.flush()
            result = self._receive(status, self._elapsed_ticks+advance_ticks)
        except BaseException:
            self.close(graceful=False)
            raise
        self._elapsed_ticks += advance_ticks
        self._last = result
        self._pending_meal=result['pending_meal']
        self._journal.write(json.dumps({'acknowledgment': result}, allow_nan=False)+'\n')
        self._journal.flush()
        return result

    def snapshot(self):
        return self._command('snapshot')

    def step(self, seconds):
        ticks = _ticks(seconds, self.config.horizon_s, 'step')
        if ticks>round(self.config.horizon_s*50)-self._elapsed_ticks:
            raise ValueError('Step exceeds the remaining native horizon')
        result=self._command('step', ticks, advance_ticks=ticks)
        return result

    def apnea(self, severity):
        number(severity, 0, 1, 'apnea severity')
        return self._command('apnea', severity)

    def exercise(self, intensity):
        number(intensity, 0, .5, 'exercise intensity')
        return self._command('exercise', intensity)

    def meal(self, meal):
        if not isinstance(meal, Meal):
            raise ValueError('Expected Meal')
        if self._pending_meal:
            raise ValueError('Native meal is pending digestion; step before delivering another meal')
        result=self._command('meal', *asdict(meal).values())
        return result

    def save_state(self):
        if self._pending_meal:
            raise ValueError('Consume the pending meal with step() before saving; native pending-meal XML is not reloadable')
        self._command('save')
        path = self.out/'states/native_session.xml'
        if not path.is_file():
            raise RuntimeError('Native save did not produce a state file')
        retained = path.with_name(f'native_session_{self._sequence:08d}.xml')
        shutil.copyfile(path, retained)
        self._journal.write(json.dumps({'saved_state_sha256': _sha(retained), 'path': str(retained), 'sequence': self._sequence})+'\n')
        self._journal.flush()
        return retained

    def close(self, graceful=True):
        if self._closed:
            return
        pending_error=graceful and self._pending_meal
        if pending_error:graceful=False
        try:
            process = getattr(self, 'process', None)
            if process is not None:
                if graceful and process.poll() is None:
                    self._command('quit', status='closed')
                if process.poll() is None:
                    if not graceful:
                        process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait(timeout=5)
                self._reader.join(timeout=5)
                process.stdin.close()
                process.stdout.close()
                (self.out/'execution.json').write_text(json.dumps({'exit_code': process.returncode,
                    'elapsed_s': self._elapsed_ticks*.02, 'graceful_close': graceful,
                    'pending_meal_at_close':self._pending_meal,
                    'unconsumed_action_error':pending_error})+'\n')
        finally:
            self._closed = True
            self._journal.close()
        if pending_error:
            raise ValueError('Closed with an unconsumed meal; native process terminated and action receipts preserved. Step before closing to consume it.')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close(graceful=exc_type is None)
