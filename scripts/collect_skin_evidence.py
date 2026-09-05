"""Import published human wound-field observations without fitting channel parameters."""
from html import unescape
import json
from pathlib import Path
import re
from ihm.forge.acquisition import fetch

URL='https://pmc.ncbi.nlm.nih.gov/articles/PMC3228273/'

def main():
    source=Path('data/raw/integumentary/human-wound-field-2011.html')
    receipt=fetch(URL,source,max_bytes=5*1024**2)
    html=source.read_text()
    if 'Human clinical trial data' not in html: raise ValueError('Expected scientific table, not access challenge')
    start=html.index('id="T1"');section=html[start:html.index('</section>',start)]
    tables=re.findall(r'<table\b[^>]*>(.*?)</table>',section,re.S)
    assert len(tables)==2
    observations=[];targets=[]
    for table,sites in zip(tables,[['control','LA-1'],['LA-2','LL-1']],strict=True):
        rows=[]
        for row in re.findall(r'<tr\b[^>]*>(.*?)</tr>',table,re.S):
            rows.append([re.sub(r'\s+',' ',unescape(re.sub('<[^>]+>',' ',cell))).strip() for cell in re.findall(r'<t[dh]\b[^>]*>(.*?)</t[dh]>',row,re.S)])
        means=next(r for r in rows if r[0]=='Mean');sems=next(r for r in rows if r[0]=='SEM')
        for col in range(8):
            context={'site':sites[col//4],'age_group_table':'18-29' if col%4<2 else '65-80','sex':'male' if col%2==0 else 'female',
                'age_discrepancy':'Methods and abstract give 18-25; table heading gives 18-29. Retain both; unresolved.',
                'units':'mV/mm','si_units':'V/m','si_scale':1,'species':'human','measurement':'lateral field above epidermis, below stratum corneum',
                'source_file':str(source),'source_sha256':receipt['sha256'],'reference':URL+'#T1','dynamic_parameter_fit':False}
            targets.append({**context,'name':'wound_field/'+sites[col//4]+'/'+context['age_group_table']+'/'+context['sex'],
                'value':float(means[col+1]),'sem':float(sems[col+1]),'basis':'published_human_experiment_summary',
                'calibration_context':'Dermacorder human clinical trial NCT00355823; published outlier-excluded mean; not a membrane voltage or channel conductance'})
            for row in rows:
                if not row[0].isdigit(): continue
                raw=row[col+1];match=re.fullmatch(r'\*?\s*(\d+)±(\d+)',raw)
                observations.append({**context,'participant_within_group':int(row[0]),'raw_value':raw,
                    'value':float(match[1]) if match else None,'sem':float(match[2]) if match else None,
                    'excluded_by_authors':'*' in raw,'missing':match is None,'replicates':'usually three scans; exact count not provided per cell'})
    assert len(observations)==160 and len(targets)==16
    out=Path('data/derived/integumentary');out.mkdir(parents=True,exist_ok=True)
    for name,rows in [('human-wound-observations',observations),('validation-targets',targets)]:
        (out/(name+'.jsonl')).write_text(''.join(json.dumps(r,allow_nan=False)+'\n' for r in rows))
    summary={'source':receipt,'table_cells':len(observations),'numeric_cells':sum(r['value'] is not None for r in observations),
        'author_excluded_cells':sum(r['excluded_by_authors'] for r in observations),'group_summary_targets':len(targets),
        'participants':40,'raw_scans_available':False,'independent_validation':False,'age_heading_inconsistent_with_methods':True}
    (out/'index.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))

if __name__=='__main__':main()
