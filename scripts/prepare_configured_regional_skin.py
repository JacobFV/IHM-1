#!/usr/bin/env python3
"""Stage a new immutable configured header/fixture only; no native build or activation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.regional_skin_configuration import build_configuration, validate_configuration
PARENT=ROOT/'data/runtime/physiology/variants/whole_body_integrity_regional_skin_graph_v2'
PARENT_MANIFEST_SHA='52de403e3dab2682f774755b7bf1370711bc5ab218244c581460b84a82aad130'
HEADER_SHA='1948643357d0d0acdfab82582af886e6c727c7938d56a65965d8674cc157dee6'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def patch_header(source,config):
    if hashlib.sha256(source.encode()).hexdigest()!=HEADER_SHA:raise ValueError('Changed accepted native header')
    fractions=validate_configuration(config)
    names=[r['name'] for r in config['regions']]
    changes=[('NativeRegionalSkin','NativeConfiguredRegionalSkin'),
             ('static constexpr std::array<double,3> fractions{.2,.3,.5};',
              'static constexpr std::array<double,3> fractions{'+','.join(format(x,'.17g') for x in fractions)+'};'),
             ('static inline const std::array<std::string,3> regions{"region_a","region_b","residual"};',
              'static inline const std::array<std::string,3> regions{'+','.join(json.dumps(x) for x in names)+'};')]
    result=source
    for old,new in changes:
        expected=2 if old=='NativeRegionalSkin' else 1
        if result.count(old)!=expected:raise ValueError('Native patch match count changed')
        result=result.replace(old,new)
    receipt='  static constexpr const char* configuration_sha256="'+config['configuration_sha256']+'";\n'
    result=result.replace('public:\n','public:\n'+receipt,1)
    restored=result.replace(receipt,'',1)
    for old,new in reversed(changes):restored=restored.replace(new,old)
    if restored!=source:raise ValueError('Native patch changed source laws')
    return result


def stage(config,output):
    validate_configuration(config)
    if sha(PARENT/'manifest.json')!=PARENT_MANIFEST_SHA:raise ValueError('Changed accepted parent manifest')
    source=(PARENT/'native_regional_skin.h').read_text()
    header=patch_header(source,config)
    fixture=(ROOT/'scripts/native_regional_skin_fixture.cpp').read_text()
    fixture=fixture.replace('native_regional_skin.h','native_configured_regional_skin.h').replace('NativeRegionalSkin','NativeConfiguredRegionalSkin')
    for old,new in [('mL_Per_s)/.2;','mL_Per_s)/split.fractions[0];'),('mL_Per_s)/.3;','mL_Per_s)/split.fractions[1];')]:
        if fixture.count(old)!=1:raise ValueError('Fixture normalization match changed')
        fixture=fixture.replace(old,new)
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    (out/'configuration.json').write_text(json.dumps(config,sort_keys=True,separators=(',',':'))+'\n')
    (out/'native_configured_regional_skin.h').write_text(header)
    (out/'native_configured_regional_skin_fixture.cpp').write_text(fixture)
    parent=json.loads((PARENT/'manifest.json').read_text())
    receipt={'schema':'configured_regional_skin_preparation_v1','parent_variant':PARENT.name,
             'parent_manifest_sha256':PARENT_MANIFEST_SHA,'parent_header_sha256':HEADER_SHA,
             'parent_library_sha256':parent['library_sha256'],'configuration_sha256':config['configuration_sha256'],
             'configured_header_sha256':sha(out/'native_configured_regional_skin.h'),
             'fixture_sha256':sha(out/'native_configured_regional_skin_fixture.cpp'),
             'parent_manifest':parent,'source_algorithm_preserved_except_class_constants_and_receipt':True,
             'native_compiled':False,'native_parity_verified':False,'native_activation_allowed':False,
             'builder_sha256':sha(Path(__file__))}
    (out/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection',type=Path,required=True,help='Explicit named unions and three engineering priors JSON')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();selection=json.loads(args.selection.read_text())
    config=build_configuration(ROOT/'data/research/engineered_skin_territories/materialization.json',**selection)
    receipt=stage(config,args.output)
    print(json.dumps({k:receipt[k] for k in ['configuration_sha256','configured_header_sha256','native_compiled']}))
