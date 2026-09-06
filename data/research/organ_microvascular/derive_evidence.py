"""Offline extraction of held primary tables and metadata; no physiology fit."""
import collections, hashlib, json, pathlib, xml.etree.ElementTree as ET, zipfile
from html.parser import HTMLParser
ROOT=pathlib.Path(__file__).resolve().parents[3]
HERE=pathlib.Path(__file__).resolve().parent
class Tables(HTMLParser):
 def __init__(self): super().__init__();self.tables=[];self.table=None;self.row=None;self.cell=None
 def handle_starttag(self,tag,attrs):
  if tag=='table':self.table=[]
  elif tag=='tr' and self.table is not None:self.row=[]
  elif tag in ('th','td') and self.row is not None:self.cell=[]
 def handle_data(self,s):
  if self.cell is not None:self.cell.append(s)
 def handle_endtag(self,tag):
  if tag in ('td','th') and self.cell is not None:self.row.append(' '.join(''.join(self.cell).split()));self.cell=None
  elif tag=='tr' and self.row is not None:self.table.append(self.row);self.row=None
  elif tag=='table' and self.table is not None:self.tables.append(self.table);self.table=None

def main():
 receipts=json.loads((HERE/'acquisition.json').read_text())
 for n,r in receipts.items():
  if 'sha256' in r:
   b=(HERE/n).read_bytes();assert len(b)==r['bytes'] and hashlib.sha256(b).hexdigest()==r['sha256'],n
 ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 with zipfile.ZipFile(HERE/'DERMA-OCTA.xlsx') as z:
  strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
  rows=[]
  for r in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row',ns):
   row={}
   for c in r:
    v=c.find('s:v',ns);v=v.text if v is not None else '';row[''.join(x for x in c.attrib['r'] if x.isalpha())]=strings[int(v)] if c.get('t')=='s' and v else v
   if row.get('A','').isdigit():rows.append(row)
 assert len(rows)==330
 healthy=[r for r in rows if r['D']=='HEALTHY']
 derma={'scans':len(rows),'subjects':len({r['B'] for r in rows}),'healthy_scans':len(healthy),'healthy_subjects':len({r['B'] for r in healthy}),'healthy_locations':dict(collections.Counter(r['C'] for r in healthy)),'disease_scan_counts':dict(collections.Counter(r['D'] for r in rows)),'metadata_only':True,'images_acquired':False}
 record=json.loads((HERE/'derma_zenodo_record.json').read_text())
 derma['listed_files']=[{'name':f['key'],'bytes':f['size'],'publisher_checksum':f['checksum'],'held':f['key']=='DERMA-OCTA.xlsx'} for f in record['files']]
 derma['listed_total_bytes']=sum(f['size'] for f in record['files'])
 muscle=ET.parse(HERE/'muscle_article.xml').getroot();table=next(t for t in muscle.findall('.//table-wrap') if t.find('label').text=='Table 3')
 muscle_rows=[[' '.join(''.join(c.itertext()).split()) for c in row] for row in table.findall('.//tr')]
 lung=Tables();lung.feed((HERE/'lung_article.html').read_text());lung_table=next(t for t in lung.tables if any('Present Study' in row for row in t))
 held=[]
 for name in ['ihm/assembly/body_microstructure.py','ihm/assembly/body_exchange.py','scripts/build_body_details.py','data/derived/canonical/microvascular.json','data/derived/microstructure/kidney/example_graph.npz','data/derived/microstructure/kidney/statistics.json','data/raw/microstructure/kidney/Test.am','data/raw/microstructure/kidney/Test2.am','data/raw/microstructure/kidney/segmentation_19685382/slab_manifest.json']:
  p=ROOT/name;b=p.read_bytes();held.append({'path':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
 out={'schema':'organ_microvascular_primary_evidence_v1','derma':derma,'muscle_table3_rows':muscle_rows,'muscle_statistic':'mean +/- SD, human vastus lateralis cross sections; Table 3 n totals 46, recruited cohort 47','lung_table2_present_study':[[r[0],r[1]] for r in lung_table if len(r)>1],'lung_scope':'13 single donor lungs, mean +/- SD; alveolar/septal statistics, not capillary diameters','prior_held_files':held,'new_downloaded_bytes':sum(r.get('bytes',0) for r in receipts.values())}
 (HERE/'derived_statistics.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
 print(json.dumps({'verified_downloads':sum('sha256' in r for r in receipts.values()),'new_downloaded_bytes':out['new_downloaded_bytes'],'derma':derma},indent=2))
if __name__=='__main__':main()
