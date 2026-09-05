"""Traceable extraction from upstream physiology repositories (not a calibration).

Only literal numbers are numeric values. Expressions, missing units, and uncertain
provenance remain explicit; no inferred value is presented as a measurement.
"""
from __future__ import annotations
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from functools import lru_cache

SOURCES = {
    'biogears': ('https://github.com/BioGearsEngine/core.git', 'Apache-2.0; third-party components retain their licenses'),
    'betse': ('https://github.com/betsee/betse.git', 'BSD-2-Clause'),
    'physiomodel': ('https://github.com/physiology/Physiomodel.git', 'Physiomodel License 1.0: GPL-3.0 public route; registered BSD route not assumed'),
    'physiolibrary': ('https://github.com/MarekMatejak/Physiolibrary.git', 'BSD-3-Clause (current checkout LICENSE); historical dependency version not matched'),
    'hummod': ('https://github.com/HumMod/hummod-standalone.git', 'XML descriptions GPL-2.0; executable and other material restricted, see README.md'),
}
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

@lru_cache(maxsize=None)
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()

def numeric(value):
    return float(value) if re.fullmatch(NUMBER, value.strip()) else None

def record(source, path, root, text, pos, name, expression, kind, unit=None, **extra):
    line = text.count('\n', 0, pos) + 1
    lines = text.splitlines()
    return dict(source=source, source_file=str(path.relative_to(root)), source_sha256=sha(path),
                line=line, name=name, expression=expression, value=numeric(expression),
                units=unit, units_status='explicit' if unit else 'unresolved', kind=kind,
                evidence_class='upstream_model_default', experimental_calibration_verified=False,
                equation_context='\n'.join(lines[max(0,line-4):line+5]), **extra)

def hummod(root):
    for path in sorted((root/'hummod'/'Structure').rglob('*.DES')):
        text = path.read_text(errors='replace')
        for match in re.finditer(r'<(parm|constant|def)>(.*?)</\1>', text, re.S):
            tag, body = match.groups()
            name = re.search(r'<name>(.*?)</name>', body, re.S)
            val = re.search(r'<val>(.*?)</val>', body, re.S)
            if name and val and name[1].strip() and val[1].strip():
                yield record('hummod', path, root, text, match.start(), name[1].strip(), val[1].strip(),
                             'parameter' if tag=='parm' else 'constant' if tag=='constant' else 'equation_definition',
                             reference_files=[str(p.relative_to(root)) for p in sorted(path.parent.glob('*.REF'))])
        for match in re.finditer(r'<curve>(.*?)</curve>', text, re.S):
            name = re.search(r'<name>(.*?)</name>', match[1], re.S)
            if not name: continue
            for j, point in enumerate(re.finditer(r'<point>(.*?)</point>', match[1], re.S)):
                for field in ('x','y','slope'):
                    val = re.search(fr'<{field}>(.*?)</{field}>', point[1], re.S)
                    if val:
                        yield record('hummod', path, root, text, match.start(), f'{name[1].strip()}[{j}].{field}', val[1].strip(), 'response_curve_coordinate')

def modelica(root):
    for path in sorted((root/'physiomodel').rglob('*.mo')):
        text = path.read_text(errors='replace')
        # Mask line comments with equal-length blanks, preserving source positions.
        clean = re.sub(r'//[^\n]*', lambda m:' '*len(m[0]), text)
        for match in re.finditer(r'\b(parameter|constant)\s+([\w.]+)\s+(\w+)\s*(\([^;]*?\))?\s*=\s*([^;]+);', clean):
            kind, typ, name, modifiers, rhs = match.groups()
            expression = re.split(r'\s+annotation\b|\s*"', rhs)[0].strip()
            if not expression: continue
            unit = re.search(r'\bunit\s*=\s*"([^"]+)"', modifiers or '')
            item = record('physiomodel', path, root, text, match.start(), name, expression, kind,
                          unit=unit[1] if unit else None, modelica_type=typ, declaration=text[match.start():match.end()])
            if not unit and typ.startswith(('Physiolibrary.Types.','Modelica.Units.SI.','Modelica.SIunits.')):
                item['units_status']='typed_requires_library_resolution'
            yield item

def biogears(root):
    for path in sorted((root/'biogears'/'share'/'data').rglob('*.xml')):
        if 'states' in path.parts: continue # dynamic snapshots are not coefficients
        text=path.read_text(errors='replace')
        for m in re.finditer(r'<([\w]+)\b([^<>]*\bvalue="[^"]+"[^<>]*)/?>',text):
            attrs=dict(re.findall(r'(\w+)="([^"]*)"',m[2]))
            if numeric(attrs['value']) is not None:
                yield record('biogears',path,root,text,m.start(),m[1],attrs['value'],'configuration_value',attrs.get('unit'))
    engine=root/'biogears/projects/biogears/libBiogears/src/engine'
    for path in sorted(engine.rglob('*.cpp')):
        text=path.read_text(errors='replace')
        # These are coefficient CANDIDATES, not a claim all literals are fit parameters.
        for m in re.finditer(fr'^\s*(?:(?:const\s+)?double\s+)?([\w]+)\s*=\s*({NUMBER})\s*;',text,re.M):
            yield record('biogears',path,root,text,m.start(),m[1],m[2],'numeric_assignment_candidate')

def betse(root):
    for path in sorted((root/'betse').rglob('*.yaml')):
        if not any(p in path.parts for p in ('data','yaml','_data')): continue
        text=path.read_text(errors='replace'); stack=[]; offset=0
        for line in text.splitlines(keepends=True):
            m=re.match(r'^(\s*)([^#\s][^:]*):\s*(.*?)\s*(?:#(.*))?$',line.rstrip())
            if m:
                depth=len(m[1]); key=m[2].strip(); val=m[3].strip()
                while stack and stack[-1][0]>=depth: stack.pop()
                name='.'.join([s[1] for s in stack]+[key])
                if not val: stack.append((depth,key))
                elif numeric(val) is not None:
                    item=record('betse',path,root,text,offset,name,val,'yaml_numeric_setting',comment=(m[4] or '').strip())
                    item['published_example']='paper' in path.parts
                    yield item
            offset+=len(line)

def validation(root):
    directory=root/'biogears/share/data/validation'
    bib=root/'biogears/share/doc/methodology/BioGears.bib'
    bibtext=bib.read_text(errors='replace')
    for path in sorted(directory.glob('*Validation.csv')):
        with path.open(newline='') as f:
            reader=csv.reader(f)
            for row in reader:
                if not row or len(row)<4: continue
                refs=[x.strip() for x in re.split('[,\n]',row[3]) if x.strip()]
                yield dict(source='biogears',source_file=str(path.relative_to(root)),source_sha256=sha(path),
                           line=reader.line_num, name=row[0], units=row[1] or None, reference_value=row[2],
                           citation_keys=refs, citation_keys_found_in_bibliography={r:bool(re.search(r'@\w+\s*\{\s*'+re.escape(r)+r'\s*,',bibtext,re.I)) for r in refs},
                           bibliography_file=str(bib.relative_to(root)),notes=row[4] if len(row)>4 else '',
                           system=row[5] if len(row)>5 else path.stem.replace('Validation',''),
                           evidence_class='literature_reference_target_as_curated_by_upstream',
                           raw_participant_measurement=False, locally_reproduced=False)

def calibration_evidence(root):
    """Locate upstream fit claims without promoting them to independently verified fits."""
    for path in sorted((root/'biogears/share/doc/methodology').glob('*Methodology.md')):
        text=path.read_text(errors='replace')
        offset=0
        for line in text.splitlines(keepends=True):
            if re.search(r'\b(fit|fitted|calibrated|tuned|tuning)\b',line,re.I):
                yield dict(source='biogears', source_file=str(path.relative_to(root)),
                           source_sha256=sha(path), line=text.count('\n',0,offset)+1,
                           upstream_statement=line.strip(),
                           citation_keys=re.findall(r'@cite\s+(\w+)',line),
                           evidence_class='upstream_calibration_statement_requires_parameter_mapping',
                           independently_reproduced=False)
            offset+=len(line)

def curated_fits(root):
    path=root/'biogears/projects/biogears/libBiogears/src/engine/Systems/Endocrine.cpp'
    text=path.read_text()
    for name,unit,unit_status in [('e50_W','W','explicit_variable_suffix'),
                                 ('eta','1/W','dimensional_inference_from_logistic_equation'),
                                 ('maxMultiplier','1','dimensionless_release_multiplier')]:
        m=re.search(fr'double {name} = ({NUMBER});',text)
        if not m: raise ValueError(f'Upstream fitted coefficient changed: {name}')
        row=record('biogears',path,root,text,m.start(),name,m[1],'upstream_fitted_parameter',unit)
        row.update(units_status=unit_status,evidence_class='upstream_reported_fit_to_published_data',
                   equation='releaseMultiplier = 1 + maxMultiplier / (1 + exp(-eta*(exercise_W-e50_W)))',
                   equation_source='biogears/projects/biogears/libBiogears/src/cdm/utils/GeneralMath.cpp:390',
                   calibration_context=text[text.index('// If we have exercise'):text.index('// If we have a stress/anxiety')],
                   citation_keys=['tidgren1991renal','stratton1985hemodynamic'],
                   methodology_source='biogears/share/doc/methodology/EndocrineMethodology.md:75',
                   experimental_calibration_verified=False,
                   limitation='Source authors report fit; original observations and fitting residuals not recovered or independently refitted.')
        yield row

def collect(base, download=False):
    base=Path(base); root=base/'data/raw/physiology'; out=base/'data/derived/physiology'
    root.mkdir(parents=True,exist_ok=True); out.mkdir(parents=True,exist_ok=True)
    for name,(url,_) in SOURCES.items():
        if not (root/name).exists() and download:
            subprocess.run(['git','clone','--depth','1',url,str(root/name)],check=True)
        if not (root/name/'.git').exists(): raise FileNotFoundError(root/name)
    inventories=[]; sources=[]
    for name,(url,license_) in SOURCES.items():
        folder=root/name
        rev=subprocess.check_output(['git','-C',str(folder),'rev-parse','HEAD'],text=True).strip()
        # Only tracked source assets: excludes .git and files generated by local installs/runs.
        files=subprocess.check_output(['git','-C',str(folder),'ls-files','-z']).decode().split('\0')
        total=0; count=0
        for filename in filter(None,files):
            p=folder/filename
            if not p.is_file(): continue
            size=p.stat().st_size; total+=size; count+=1
            inventories.append(dict(source=name,path=str(p.relative_to(base)),size_bytes=size,sha256=sha(p),revision=rev))
        sources.append(dict(source=name,url=url,revision=rev,license=license_,tracked_files=count,bytes=total))
    if sum(s['bytes'] for s in sources)>3_000_000_000: raise RuntimeError('Source assets exceed 3 GB cap')
    counts=Counter(); valuecounts=Counter(); n=0
    with (out/'coefficients.jsonl').open('w') as f:
        for parser in (hummod,modelica,biogears,betse):
            for row in parser(root):
                f.write(json.dumps(row,ensure_ascii=False)+'\n'); n+=1; counts[row['source']]+=1; valuecounts[row['kind']]+=1
    targets=list(validation(root))
    fits=list(curated_fits(root))
    evidence=list(calibration_evidence(root))
    for name,rows in [('assets.jsonl',inventories),('validation_targets.jsonl',targets),('reported_fitted_coefficients.jsonl',fits),('calibration_evidence.jsonl',evidence)]:
        with (out/name).open('w') as f:
            for row in rows:f.write(json.dumps(row,ensure_ascii=False)+'\n')
    summary=dict(collected_at=datetime.now(timezone.utc).isoformat(),sources=sources,coefficient_records=n,
                 by_source=dict(counts),by_kind=dict(valuecounts),validation_targets=len(targets),upstream_reported_fits=len(fits),calibration_statements=len(evidence),
                 calibration_claim='Extraction only. No coefficients are certified experimentally fitted by this collector.',
                 limitations=['Numeric assignment candidates include algorithm constants and initial values.',
                              'Modelica parameter regex is an index, not a full language parser.',
                              'Units are unresolved unless explicit; typed units require library resolution.',
                              'Upstream validation output is simulated data, not participant measurement.',
                              'No engine coupling into IHM runtime is performed by collection.'])
    (out/'collection.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

def verify(base):
    """Check downloaded bytes, derived provenance and independent source anchors."""
    base=Path(base); raw=base/'data/raw/physiology'; out=base/'data/derived/physiology'
    assets=[json.loads(x) for x in (out/'assets.jsonl').read_text().splitlines()]
    sha.cache_clear()
    for row in assets:
        path=base/row['path']
        if path.stat().st_size!=row['size_bytes'] or sha(path)!=row['sha256']:
            raise ValueError(f'Asset integrity failure: {path}')
    hashes={str((base/x['path']).relative_to(raw)):x['sha256'] for x in assets}
    rows=[json.loads(x) for x in (out/'coefficients.jsonl').read_text().splitlines()]
    for row in rows:
        if row['source_sha256']!=hashes.get(row['source_file']) or not row['name'] or not row['expression']:
            raise ValueError(f'Invalid coefficient provenance: {row}')
        if row['experimental_calibration_verified']:
            raise ValueError('Generic coefficient extractor must not certify calibration')
    anchors={
        'BioGears e50_W=190':any(r['source']=='biogears' and r['name']=='e50_W' and r['value']==190 for r in rows),
        'Physiomodel ventricular EDV=70e-6':any(r['source']=='physiomodel' and r['name']=='NormalEndDiastolicVolume' and r['value']==70e-6 for r in rows),
        'HumMod ADH/Kidney/Skin source coverage':all(any(k in r['source_file'] for r in rows) for k in ('ADH/','Kidney','Skin')),
    }
    if not all(anchors.values()): raise ValueError(anchors)
    report=dict(status='passed',tracked_files_hash_verified=len(assets),coefficient_provenance_rows_verified=len(rows),
                source_anchor_checks=anchors,claim='Artifact integrity and parser anchors, not physiological validation')
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
