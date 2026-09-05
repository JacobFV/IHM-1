"""Acquire pinned authors' JOS-3 and execute source-qualified one-hour cases."""
from pathlib import Path
import sys,json,hashlib,subprocess,urllib.request
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.thermal import REVISION,SOURCE_URL,BEDDING_DOI,PROFILES,run_thermal
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,separators=(',',':'),allow_nan=False))
def main():
    raw=ROOT/'data/raw/thermal';raw.mkdir(parents=True,exist_ok=True);repo=raw/'JOS-3'
    if not repo.exists():subprocess.run(['git','clone',SOURCE_URL+'.git',str(repo)],check=True)
    actual=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
    if actual!=REVISION:raise ValueError('Checkout must be pinned to '+REVISION)
    if subprocess.check_output(['git','-C',str(repo),'status','--porcelain'],text=True).strip():raise ValueError('JOS-3 source checkout is modified')
    papers=[('bedding.pdf','https://researchmap.jp/shinichitanabe/published_papers/50138511/attachment_file.pdf',BEDDING_DOI),('jos3-paper.pdf','https://researchmap.jp/shinichitanabe/published_papers/30582046/attachment_file.pdf','10.1016/j.enbuild.2020.110575')]
    paper_metadata=[]
    for name,url,doi in papers:
        path=raw/name
        if not path.exists():urllib.request.urlretrieve(url,path)
        if not path.read_bytes().startswith(b'%PDF-'):raise ValueError('not a PDF: '+name)
        paper_metadata.append(dict(path=str(path.relative_to(ROOT)),url=url,doi=doi,sha256=sha(path)))
    sources={str(p.relative_to(repo)):sha(p) for p in [*sorted((repo/'src/jos3').glob('*.py')),repo/'LICENSE',repo/'README.md']}
    provenance=dict(url=SOURCE_URL,revision=REVISION,license='MIT',version='0.5.0',source_sha256=sources,papers=paper_metadata,adapter_path='ihm/native/thermal.py',adapter_sha256=sha(ROOT/'ihm/native/thermal.py'))
    write(raw/'source.json',provenance)
    out=ROOT/'data/derived/thermal';out.mkdir(parents=True,exist_ok=True);runs=[]
    for profile in PROFILES:
        cases={dt:run_thermal(ROOT,profile,dt=dt) for dt in [60,30,15,7.5,3.75,1.875,.9375]}
        fine=cases[.9375]
        refinement={}
        for coarse,finer in [(60,30),(30,15),(15,7.5),(7.5,3.75),(3.75,1.875),(1.875,.9375)]:
            delta=np.asarray(cases[coarse]['node_temperature_C'])-np.asarray(cases[finer]['node_temperature_C'])[::int(coarse/finer)]
            refinement[f'dt{coarse}_vs_dt{finer}_max_node_temperature_C']=float(abs(delta).max())
            refinement[f'dt{coarse}_vs_dt{finer}_final_max_node_temperature_C']=float(abs(delta[-1]).max())
            refinement[f'dt{coarse}_vs_dt{finer}_after60s_max_node_temperature_C']=float(abs(delta[int(60/coarse):]).max())
        refinement['interpretation']='Common-clock step refinement. Early maxima at moving first-step times are not uniformly monotonic; endpoint and post60s errors reported separately. No clinical accuracy inference.'
        fine['refinement']=refinement
        coefficient=fine.pop('coefficients');write(out/f'{profile}-coefficients.json',coefficient)
        write(out/f'{profile}.json',fine)
        runs.append(dict(id=profile,label=fine['label'],trajectory_path=str((out/f'{profile}.json').relative_to(ROOT)),coefficient_path=str((out/f'{profile}-coefficients.json').relative_to(ROOT)),configuration=fine['configuration'],refinement=refinement,audit={k:v for k,v in fine['audit'].items() if k!='steps'},final={c['id']:c['values'][-1] for c in fine['channels'] if c['id'] in ['MeanSkinTemperature','CentralBloodTemperature','Met','RES']},limitations=fine['limitations']))
        print(profile,refinement,runs[-1]['audit'],runs[-1]['final'],flush=True)
    index=dict(schema_version=1,id='jos3_thermal',source_kind='source_model_simulation',source=provenance,node_count=coefficient['node_count'],body_regions=coefficient['body_regions'],runs=runs)
    write(out/'index.json',index)
    # Standalone scientific comparison of independently sourced cases.
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for run in runs:
        data=json.loads((ROOT/run['trajectory_path']).read_text());values={x['id']:x['values'] for x in data['channels']}
        for ax,key in zip(axes,['CentralBloodTemperature','MeanSkinTemperature']):ax.plot(np.array(data['time_s'])/60,values[key],label=run['id']);ax.set(xlabel='Time after lying/boundary change (min)',ylabel=key+' (°C)');ax.grid(alpha=.2)
    axes[1].legend(fontsize=7);fig.suptitle('Pinned JOS-3 source simulation · distinct defaults and measured total-insulation boundaries')
    fig.savefig(out/'comparison.png',dpi=160);fig.savefig(out/'comparison.svg')
if __name__=='__main__':main()
