"""Prepare/resume bounded donor-object Simbody instance-mass diagnostic variant.
Default preparation is source-only. Every native chunk needs the heavy slot.
No installation or latest pointer is promoted; completed library remains opt-in.
"""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,time
from patch_simbody_instance_mass import patch,ROOT
RAW=ROOT/'data/raw/mechanics/simbody/Simbody/src'
BUILD=ROOT/'data/runtime/opensim/simbody-build/Simbody'
BASE=ROOT/'data/runtime/opensim/variants'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def prepare(out):
    if out.exists():raise ValueError('Fresh owned variant directory required')
    out.mkdir(parents=True);source=out/'src';shutil.copytree(RAW,source)
    changed=patch(source);objects=out/'objects';objects.mkdir()
    objectbase=BUILD/'CMakeFiles/SimTKsimbody.dir/src'
    link=shlex.split((BUILD/'CMakeFiles/SimTKsimbody.dir/link.txt').read_text())
    deps={};rebuild=[];donors={};inputs={}
    for token in link:
        if not token.endswith('.o'):continue
        donor=(BUILD/token).resolve();name=donor.name
        dependency=donor.with_suffix(donor.suffix+'.d');raw=dependency.read_text()
        needed=any('/'+n in raw for n in changed)
        # Pin retained object and dependency inventory even for rebuilt units.
        inputs[str(donor.relative_to(ROOT))]=sha(donor);inputs[str(dependency.relative_to(ROOT))]=sha(dependency)
        for part in raw.replace('\\\n',' ').split():
            p=Path(part)
            if p.is_file() and p.is_relative_to(ROOT):inputs[str(p.relative_to(ROOT))]=sha(p)
        if needed:rebuild.append(name)
        else:shutil.copy2(donor,objects/name);donors[name]=sha(objects/name)
        deps[name]=str(donor.relative_to(ROOT))
    for token in link:
        candidate=(BUILD/token).resolve()
        if not token.startswith('-') and candidate.is_file() and '.so' in candidate.name:
            inputs[str(candidate.relative_to(ROOT))]=sha(candidate)
    flags={}
    for line in (BUILD/'CMakeFiles/SimTKsimbody.dir/flags.make').read_text().splitlines():
        if line.startswith(('CXX_DEFINES =','CXX_INCLUDES =','CXX_FLAGS =')):
            k,v=line.split('=',1);flags[k.strip()]=shlex.split(v.strip())
    # O1 is deliberately isolated from the retained O2 donor object inventory.
    cflags=[x for x in flags['CXX_FLAGS'] if not x.startswith('-O')]+['-O1']
    inputs[str((ROOT/'scripts/patch_simbody_instance_mass.py').relative_to(ROOT))]=sha(ROOT/'scripts/patch_simbody_instance_mass.py')
    inputs[str((ROOT/'scripts/simbody_instance_mass_source_identity.json').relative_to(ROOT))]=sha(ROOT/'scripts/simbody_instance_mass_source_identity.json')
    manifest={'schema':'ihm.simbody-instance-mass-variant.v1','complete':False,'source_revision':json.loads((ROOT/'scripts/simbody_instance_mass_source_identity.json').read_text())['revision'],'changed_files':changed,'rebuild':sorted(rebuild),'donor_objects':donors,'original_objects':deps,'built_objects':{},'inputs':inputs,'source_overlay_sha256':{p.name:sha(p) for p in source.iterdir() if p.is_file()},'compiler':'/usr/bin/c++','compile_flags':flags['CXX_DEFINES']+flags['CXX_INCLUDES']+cflags,'link_command':link,'chunks':[],'scope':'Opt-in native source variant, no installed-library replacement; O1 replacement objects with untouched O2 donors'}
    dump(out/'manifest.json',manifest);return manifest

def validate(out,m):
    for name,digest in m['inputs'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen donor/source changed: '+name)
    for name,digest in m['source_overlay_sha256'].items():
        if sha(out/'src'/name)!=digest:raise ValueError('Frozen overlay changed: '+name)
    for name,digest in {**m['donor_objects'],**m['built_objects']}.items():
        if sha(out/'objects'/name)!=digest:raise ValueError('Object changed: '+name)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--build-native',action='store_true');parser.add_argument('--max-objects',type=int,default=3);parser.add_argument('--link',action='store_true');parser.add_argument('--only-object');a=parser.parse_args()
    out=(ROOT/a.output).resolve()
    if not out.is_relative_to(BASE):raise ValueError('Variant must be under owned opensim/variants')
    if not 1<=a.max_objects<=3:raise ValueError('Native chunk is one to three translation units')
    if a.link and not a.build_native:raise ValueError('Link requires coordinated --build-native')
    m=json.loads((out/'manifest.json').read_text()) if out.exists() else prepare(out)
    validate(out,m)
    remaining=[x for x in m['rebuild'] if x not in m['built_objects']]
    if a.build_native and not m['complete']:
        limiter=shutil.which('prlimit');assert limiter,'prlimit required'
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        chunk={'objects':[],'started_unix_s':time.time()}
        if a.only_object and a.only_object not in remaining:raise ValueError('Requested object is not pending')
        selected=[a.only_object] if a.only_object else remaining[:a.max_objects]
        for name in selected:
            command=[limiter,'--as=4294967296','--','nice','-n','10',m['compiler'],*m['compile_flags'],'-I',str(out/'src'),'-c',str(out/'src'/name.removesuffix('.o')),'-o',str(out/'objects'/name)]
            start=time.monotonic();dump(out/(name+'.command.json'),command)
            with (out/(name+'.log')).open('w') as log:result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode:raise RuntimeError('Native object failed: '+str(out/(name+'.log')))
            validate(out,m);m['built_objects'][name]=sha(out/'objects'/name)
            chunk['objects'].append({'name':name,'wall_s':time.monotonic()-start,'max_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss})
            dump(out/'manifest.json',m)
        m['chunks'].append(chunk);dump(out/'manifest.json',m)
        remaining=[x for x in m['rebuild'] if x not in m['built_objects']]
        if a.link:
            if remaining:raise ValueError('Cannot link incomplete object inventory')
            original=m['link_command'];command=[];skip=False
            for i,token in enumerate(original):
                if skip:skip=False;continue
                if token=='-o':command+=['-o',str(out/'libSimTKsimbody.so.3.9')];skip=True
                elif token.endswith('.o'):command.append(str(out/'objects'/Path(token).name))
                elif not token.startswith('-') and (BUILD/token).is_file():command.append(str((BUILD/token).resolve()))
                else:command.append(token)
            command=[limiter,'--as=4294967296','--','nice','-n','10',*command]
            dump(out/'link.command.json',command)
            with (out/'link.log').open('w') as log:result=subprocess.run(command,cwd=BUILD,env=env,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode:raise RuntimeError('Variant link failed: '+str(out/'link.log'))
            validate(out,m);m['complete']=True;m['library_sha256']=sha(out/'libSimTKsimbody.so.3.9');dump(out/'manifest.json',m)
    print(json.dumps({'output':str(out),'complete':m['complete'],'rebuild_total':len(m['rebuild']),'built':len(m['built_objects']),'remaining':remaining,'native_build_requested':a.build_native},indent=2))
if __name__=='__main__':main()
