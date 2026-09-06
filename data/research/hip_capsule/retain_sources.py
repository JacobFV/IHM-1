"""Receipt existing bounded source bytes; no network, native or model mutation."""
import hashlib,json
from pathlib import Path
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent
ENTRIES=[
 ('van_Arkel_2015_envelope.xml','PMC4655836','CC BY4.0',
  '10 fresh-frozen pelvises,6male; one morphology exclusion and one capsule rupture leave8 analyzed rotation specimens; age76±9years; left hips.',
  'Skeletonized with capsule preserved;110N compression20degmedial/proximal relative to femur;roomtemperature,water spray;10s sinusoidal cycle,second iteration;±5Nm at30 flexion/abduction positions.',
  'Numeric morphology and significance tables, plot means/CIs, mid-slack regression; no raw specimen torque-angle files or covariance found. Aggregate includes labrum; not capsule-only.'),
 ('van_Arkel_2015_resection.xml','PMC4491667','CC BY (version not specified in retained permissions)',
  'Nine fresh-frozen cadaveric left hips,6male,age76.4years(range61–89).',
  'Sequential resection with return to original intact±5Nm rotational angles across18 flexion/abduction positions;110N functional compression; superposition estimates contributions.',
  'Capsular ligaments dominate restraint but labrum/teres are present. Resection contributions do not identify a global elastic potential or independent constitutive laws.'),
 ('Karunaseelan_2021_stabilization.xml','PMC8479567','CC BY-NC-ND4.0',
  'Nine cadaveric hips from earlier resection study; not an independent new cohort.',
  'MATLAB/Twente2.0 enthesis centroid scaling plus ArtiSynth wrapping; effective moment arms; stiffness inferred from4–5Nm endpoint region.',
  'Modeled forces and attachment sensitivity; no executable MATLAB/ArtiSynth model or input deck recovered. Unchanged XML retained for research attribution; do not propagate this restrictive license to original code.'),
 ('Anantha_Krishnan_2024_model.xml','PMC10813259','CC BY4.0',
  'Base model one implanted cadaver hip; property priors from10 prior calibrated specimens; independent validation5THA hips.',
  'Fiber-reinforced membrane with six sectors and tension-only springs; synthetic500-case dataset; experimental0–5Nm I/E at0,30,60,90degflexion.',
  'Article explicitly says data available upon request; no downloadable calibrated model or regression weights recovered. THA geometry and sector properties cannot be transplanted into native hip.'),
 ('Pieroh_2016_material.xml','PMC5042535','CC BY (version not specified in retained permissions)',
  'Reported17 cadavers,age83.65±10.54years; authors report9female/7male and inconsistent specimen counts in places; preserve rather than correct source.',
  'Ethanol-glycerin embalming,excised direction-selected ligament strips; saline storage,roomtemperature uniaxial tests.',
  'Material data include preparation-sensitive stress/strain values. Do not combine with fresh-frozen whole-joint torsional stiffness without geometry/slack-length and preparation model.')]


def generate():
    cards=[];files=[]
    for filename,pmcid,license,specimens,conditions,scope in ENTRIES:
        path=HERE/filename;raw=path.read_bytes();root=ET.fromstring(raw)
        url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
        doi=root.find('.//article-id[@pub-id-type="doi"]').text
        cards.append({'file':filename,'title':''.join(root.find('.//article-title').itertext()),'doi':doi,
            'primary_article_url':f'https://doi.org/{doi}','retrieval_url':url,'license':license,
            'license_verbatim':''.join(root.find('.//permissions').itertext()),
            'specimens':specimens,'conditions':conditions,'usable_scope_and_limitations':scope})
        files.append({'file':filename,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'retrieval_url':url})
    for filename,url in [
        ('van_Arkel_2015_envelope.pdf','https://europepmc.org/articles/PMC4655836?pdf=render'),
        ('Myers_2020_metadata.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:32432892%20AND%20SRC:MED&format=json&resultType=core'),
        ('envelope_page4_1800.png','Local rendering of retained CC BY4.0 PDF page4: pdftoppm -f4 -singlefile -scale-to1800 -png'),
        ('envelope_page5_1800.png','Local rendering of retained CC BY4.0 PDF page5: pdftoppm -f5 -singlefile -scale-to1800 -png'),
        ('neutral_slice.json','Local original research interpolation with explicit figure-derived calibration')]:
        raw=(HERE/filename).read_bytes();files.append({'file':filename,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'retrieval_url':url})
    meta=json.loads((HERE/'Myers_2020_metadata.json').read_bytes())['resultList']['result'][0]
    cards.append({'file':'Myers_2020_metadata.json','title':meta['title'],'doi':'10.1080/10255842.2020.1764543',
        'primary_article_url':'https://doi.org/10.1080/10255842.2020.1764543',
        'license':'Metadata retained; no fulltext/model reuse license established',
        'usable_scope_and_limitations':'Abstract identifies calibrated fiber-reinforced FE capsule and4.7degmean/5.1degprobabilistic RMS. EuropePMC lists subscription-required fulltext; no model files or calibrated parameter table acquired.'})
    return {'schema':'ihm.hip-capsule-source-cards.v1','native_activation_allowed':False,'acquisition_date':'2026-09-06',
        'source_cards':cards,'retained_files':files,
        'decision':'Retain one provisional conservative scalar research slice only. Full3DOF passive potential and subject-specific activation are not identified by available data.',
        'prohibited_inferences':['Do not treat±40degmodelROM as capsule slack or force onset.',
            'Do not set production stiffness to remove a support residual.',
            'Do not relax ROM bounds or apply fixed-pose torque at another flexion/abduction.',
            'Do not add aggregate periarticular torque on top of separately modeled labrum/ligament forces without overlap accounting.']}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    result=generate();path=HERE/'sources.json'
    if a.check:
        if json.loads(path.read_bytes())!=result:raise SystemExit('Source receipt mismatch')
        print('Verified',len(result['retained_files']),'retained files and',len(result['source_cards']),'source cards')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
