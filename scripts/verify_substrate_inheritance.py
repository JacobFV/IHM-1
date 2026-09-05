"""Run retained inheritance harness with explicit substrate variant substitution."""
import difflib,json
from pathlib import Path
from audit_native_substrate_integrity import ROOT,RUNTIME,sha

def main():
 original=ROOT/'scripts/verify_final_thermal_inheritance.py';before=original.read_text();variant='whole_body_integrity_substrate_availability';library=RUNTIME/'variants'/variant/'libbiogears.so.8.0.0'
 after=before.replace("VARIANT='whole_body_integrity_evaporation_humidity'","VARIANT='"+variant+"'").replace('ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a',sha(library)).replace("prefix='thermal-final-inheritance-'","prefix='substrate-inheritance-'")
 namespace={'__name__':'retained_substrate_inheritance','__file__':str(original)};exec(compile(after,str(original),'exec'),namespace);out=namespace['run']()
 (out/'inherited_harness.py').write_text(after);(out/'harness.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True))))
 (out/'harness_receipt.json').write_text(json.dumps({'original_sha256':sha(original),'adapted_sha256':sha(out/'inherited_harness.py'),'wrapper_sha256':sha(Path(__file__)),'change':'Only variant, corresponding required hash, and fresh output prefix. Existing harness retains native tests and documented historical parity omission.'},indent=2)+'\n')
if __name__=='__main__':main()
