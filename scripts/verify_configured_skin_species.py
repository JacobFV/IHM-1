#!/usr/bin/env python3
"""Prepare/run bounded actual native liquid transport and passive Albumin tests."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile
from prepare_configured_regional_skin import ROOT,PARENT,PARENT_MANIFEST_SHA,sha

SOURCE=ROOT/'data/raw/physiology/biogears'
RUNTIME=ROOT/'data/runtime/physiology'
PREPARED=ROOT/'data/research/configured_regional_skin/prepared_v1'
SNAPSHOT=ROOT/'data/derived/audits/regional-python-vvdedt24/initial.json'
SNAPSHOT_SHA='cb37d4e4f5460d48e2938a374ab5a92148890559f7070482926732b2903e0b4d'
STAGE=ROOT/'data/research/configured_skin_species/prepared_v2'


def prepare():
    if sha(PARENT/'manifest.json')!=PARENT_MANIFEST_SHA or sha(SNAPSHOT)!=SNAPSHOT_SHA:
        raise ValueError('Changed accepted parent or actual species seed snapshot')
    prior=json.loads((PREPARED/'preparation.json').read_text())
    inputs={'native_configured_regional_skin.h':PREPARED/'native_configured_regional_skin.h',
            'native_configured_regional_skin_fixture.cpp':PREPARED/'native_configured_regional_skin_fixture.cpp',
            'native_configured_skin_species_fixture.cpp':ROOT/'scripts/native_configured_skin_species_fixture.cpp'}
    if sha(inputs['native_configured_regional_skin.h'])!=prior['configured_header_sha256'] or sha(inputs['native_configured_regional_skin_fixture.cpp'])!=prior['fixture_sha256']:
        raise ValueError('Changed accepted configured header/fixture')
    values=json.loads(SNAPSHOT.read_text())['values'];prefix='tissue.regional_skin.aggregate.species.'
    species={k[len(prefix):-len('.mass_ug')]:v for k,v in values.items() if k.startswith(prefix) and k.endswith('.mass_ug')}
    volume=values['tissue.regional_skin.aggregate.volume_ml']
    seed='static const std::array<SpeciesSeed,'+str(len(species))+'> species_seeds{{\n'
    seed+=''.join('  {'+json.dumps(name)+','+format(mass/volume,'.17g')+'},\n' for name,mass in sorted(species.items()))+'}};\n'
    blobs={name:path.read_bytes() for name,path in inputs.items()};blobs['configured_species_seed.inc']=seed.encode()
    operators=[SOURCE/'projects/biogears/libBiogears/include/biogears/cdm/substance/SESubstanceTransport.inl',
               ROOT/'data/runtime/physiology/variants/whole_body_integrity_regional_skin_gi_v1/Diffusion.cpp']
    receipt={'schema':'configured_skin_species_preparation_v1','parent_variant':PARENT.name,
             'parent_manifest_sha256':PARENT_MANIFEST_SHA,'parent_library_sha256':prior['parent_library_sha256'],
             'configuration_sha256':prior['configuration_sha256'],'species_seed_snapshot_sha256':SNAPSHOT_SHA,
             'species_count':len(species),'species_names':sorted(species),
             'source_receipts':{str(p.relative_to(ROOT)):sha(p) for p in [*inputs.values(),*operators]},
             'prepared_sha256':{name:hashlib.sha256(raw).hexdigest() for name,raw in blobs.items()},
             'operator_scope':'Separate SELiquidTransporter and native passive Albumin phases; no duplicate operator application, sweat or filtration links',
             'native_executed':False,'production_activation_allowed':False}
    if STAGE.exists():
        for name,raw in blobs.items():
            if (STAGE/name).read_bytes()!=raw:raise ValueError('Frozen preparation differs: '+name)
        if json.loads((STAGE/'preparation.json').read_text())!=receipt:raise ValueError('Frozen species preparation receipt differs')
    else:
        STAGE.mkdir(parents=True)
        for name,raw in blobs.items():(STAGE/name).write_bytes(raw)
        (STAGE/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--compile-and-run',action='store_true');args=parser.parse_args()
    receipt=prepare()
    if not args.compile_and_run:
        print(json.dumps({'prepared':True,'native_executed':False,'species_count':receipt['species_count']}));return
    if sha(PARENT/'libbiogears.so.8.0.0')!=receipt['parent_library_sha256']:raise ValueError('Changed parent library')
    out=Path(tempfile.mkdtemp(prefix='configured-skin-species-',dir=ROOT/'data/derived/audits'))
    for p in STAGE.iterdir():shutil.copyfile(p,out/p.name)
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'fixture'
    cmd=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_configured_skin_species_fixture.cpp')]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        cmd+=['-I',str(p)]
    cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(PARENT)+':'+str(lib)}
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    def run(command,label,timeout):
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            result=subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+command,cwd=out,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits,timeout=timeout)
        if result.returncode:raise RuntimeError(label+' failed; retained '+str(out))
    print(str(out),flush=True);run(cmd,'compile',75)
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        (out/name).symlink_to(build/'runtime'/name)
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
    if str(PARENT/'libbiogears.so.8.0.0') not in linkage:raise ValueError('Wrong native lineage')
    run([str(binary)],'fixture',10)
    line=next(line for line in (out/'fixture.stdout').read_text().splitlines() if line.startswith('RESULT '));result=json.loads(line[7:])
    if not result['passed'] or result['species_count']!=receipt['species_count']:raise ValueError('Species fixture validation failed')
    result.update(preparation=receipt,compile_command=cmd,fixture_sha256=sha(binary),
                  scope='Native species operator/circuit fixture only; no full-patient composition or physiological validation of regional chemical independence')
    prepare()
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
