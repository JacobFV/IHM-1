"""Source-target deviations with exact native quantity bindings and no pass criteria."""
from pathlib import Path
import csv,hashlib,json,re
import numpy as np

NUMBER=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
UNITS={
    'mmHg':('pressure',133.322387415,0), 'Pa':('pressure',1,0),
    'mL':('volume',1e-6,0),'L':('volume',1e-3,0),'m^3':('volume',1,0),
    'mL/min':('flow',1e-6/60,0),'L/min':('flow',1e-3/60,0),'L/day':('flow',1e-3/86400,0),
    '1/min':('frequency',1/60,0),'pmol/min':('amount_rate',1e-12/60,0),
    'W':('power',1,0),'kcal/day':('power',4184/86400,0),
    'degC':('temperature',1,273.15),'K':('temperature',1,0),
    'g/s':('mass_rate',.001,0),'dimensionless':('dimensionless',1,0)}
LIMITATIONS=[
    'Targets are literature reference values curated by upstream BioGears, not local participant measurements or an independently assembled validation cohort.',
    'No tolerance, clinical pass/fail classification, aggregate success percentage or coefficient calibration is inferred from a source point or range.',
    'The final 60 seconds are a descriptive native output window; serial samples are not independent subjects and no confidence interval is manufactured.',
    'One-hour cool-room rest is a nonstationary source experiment; thermal conditions and gastrointestinal water redistribution affect the comparison.',
    'Mixed point/range or multiple-range strings remain unparsed; citations do not establish that different reference populations match the simulated subject.',
    'Native instantaneous lung volume and arterial pulse outputs at 1 Hz can alias; undersampled oscillatory quantities are excluded from comparisons.'
]


def parse_target(text):
    text=str(text).strip()
    if re.fullmatch(NUMBER,text):
        value=float(text)
        if np.isfinite(value):return dict(kind='point',value=value)
    match=re.fullmatch(r'\[\s*('+NUMBER+r')\s*,\s*('+NUMBER+r')\s*\]',text)
    if match:
        low,high=map(float,match.groups())
        if np.isfinite([low,high]).all() and low<=high:return dict(kind='range',low=low,high=high)
    return None


def convert(values,source_unit,target_unit):
    source_unit=(source_unit or 'dimensionless').replace(' ','');target_unit=(target_unit or 'dimensionless').replace(' ','')
    if source_unit not in UNITS or target_unit not in UNITS:raise ValueError('unsupported unit')
    ds,ss,os=UNITS[source_unit];dt,st,ot=UNITS[target_unit]
    if ds!=dt:raise ValueError('incompatible quantity dimensions')
    return (np.asarray(values,float)*ss+os-ot)/st


def compare_run(path,targets,root):
    path=Path(path);root=Path(root)
    with path.open() as f:
        reader=csv.reader(f);headers=next(reader);a=np.array([[float(v) for v in row] for row in reader if row])
    if headers[0]!='#Time(s)' or a.ndim!=2 or a.shape[1]!=len(headers) or not np.isfinite(a).all() or np.any(np.diff(a[:,0])<=0):raise ValueError('finite native CSV with explicit increasing seconds required')
    t=a[:,0];dt=float(np.median(np.diff(t)));selected=t>t[-1]-60
    columns={}
    for i,h in enumerate(headers[1:],1):
        match=re.fullmatch(r'(.+)\(([^()]*)\)',h);name,unit=match.groups() if match else (h,None)
        if name in columns:raise ValueError('ambiguous duplicate native quantity')
        columns[name]=(i,unit)
    rows=[]
    for target in targets:
        name=target['name'];parsed=parse_target(target['reference_value']);row=dict(target=target,parsed_target=parsed,binding='exact quantity name',status='unmapped_quantity')
        if name in columns:
            i,source_unit=columns[name];row.update(native_column=headers[i],native_unit=source_unit)
            if dt>=1-1e-6 and name in ('ArterialPressure','TotalLungVolume'):row['status']='undersampled_oscillatory_quantity'
            elif parsed is None:row['status']='unparsed_target'
            else:
                try:values=convert(a[selected,i],source_unit,target['units'])
                except ValueError as error:row.update(status='unsupported_unit_binding',reason=str(error))
                else:
                    mean=float(values.mean());low=float(values.min());high=float(values.max())
                    row.update(status='compared',comparison_unit=target['units'] or 'dimensionless',observed_mean=mean,observed_min=low,observed_max=high)
                    if parsed['kind']=='point':
                        difference=mean-parsed['value'];relative=None
                        if name not in ('BloodPH','ArterialBloodPH') and (target['units'] or '').replace(' ','') not in ('degC','K') and parsed['value']!=0:relative=difference/abs(parsed['value'])
                        row.update(signed_difference=difference,absolute_difference=abs(difference),relative_difference=relative)
                    else:
                        distance=mean-parsed['low'] if mean<parsed['low'] else mean-parsed['high'] if mean>parsed['high'] else 0.
                        row.update(signed_distance_to_reference_range=distance,observed_mean_in_reference_range=parsed['low']<=mean<=parsed['high'])
        rows.append(row)
    config=path.parent/'configuration.json'
    return dict(id=path.parent.name,source=str(path.relative_to(root)),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        configuration=json.loads(config.read_text()) if config.exists() else None,
        configuration_sha256=hashlib.sha256(config.read_bytes()).hexdigest() if config.exists() else None,source_samples=len(t),source_end_s=float(t[-1]),
        sample_interval_s=dt,window=dict(rule='source time > final time minus 60 seconds',start_s=float(t[selected][0]),end_s=float(t[-1]),samples=int(selected.sum())),
        coverage={status:sum(row['status']==status for row in rows) for status in sorted({row['status'] for row in rows})},rows=rows)


def audit_native_targets(root):
    root=Path(root).resolve();source=root/'data/derived/physiology/validation_targets.jsonl';targets=[json.loads(line) for line in source.read_text().splitlines() if line.strip()]
    paths=[root/'data/derived/physiology/biogears_native_run/native_multisystem.csv',root/'data/derived/physiology/native_baseline_v2/native_multisystem.csv',*sorted((root/'data/derived/physiology').glob('native_hour_*/native_multisystem.csv'))]
    runs=[];skipped=[]
    for path in paths:
        if not path.exists():continue
        run=compare_run(path,targets,root);config=run['configuration'] or {}
        if config.get('seconds') and run['source_end_s']<config['seconds']-run['sample_interval_s']/2:
            skipped.append(dict(id=run['id'],reason='run has not reached configured horizon',source_end_s=run['source_end_s'],configured_seconds=config['seconds']));continue
        runs.append(run)
    result=dict(schema_version=1,kind='upstream_reference_target_deviation_audit',targets_source=str(source.relative_to(root)),targets_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),target_rows=len(targets),runs=runs,skipped_runs=skipped,limitations=LIMITATIONS)
    out=root/'data/derived/calibration/native-target-audit.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,separators=(',',':'),allow_nan=False));return result
