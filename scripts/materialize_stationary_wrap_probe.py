"""Private BIClong geodesic clones layered on the accepted BRA-private variant."""
from pathlib import Path
import re,json,hashlib
from materialize_arm26_wrap_probe import materialize as cylinder_materialize,ROOT,SOURCE

def materialize(output):
    output=Path(output);parent=cylinder_materialize(output);before=(output/'candidate.osim').read_bytes();changed=before;clones=[];refs=[]
    import xml.etree.ElementTree as ET
    root=ET.fromstring(before)
    for muscle in root.findall('.//ForceSet/objects/*'):
        if not muscle.get('name','').startswith('arm26_BIClong_'):continue
        name=muscle.findtext('.//wrap_object');new=name+'__biclong_geodesic'
        pattern=rb'<ObservedArm26Ellipsoid name="'+re.escape(name.encode())+rb'">.*?</ObservedArm26Ellipsoid>'
        matches=list(re.finditer(pattern,changed,re.S));assert len(matches)==1;hit=matches[0]
        clone=hit.group().replace(b'ObservedArm26Ellipsoid',b'StationaryArm26Ellipsoid').replace(('name="'+name+'"').encode(),('name="'+new+'"').encode(),1);clones.append(clone)
        changed=changed[:hit.end()]+clone+changed[hit.end():]
        mp=rb'<Thelen2003Muscle name="'+muscle.get('name').encode()+rb'">.*?</Thelen2003Muscle>';hit=next(re.finditer(mp,changed,re.S))
        oldref=('<wrap_object>'+name+'</wrap_object>').encode();newref=('<wrap_object>'+new+'</wrap_object>').encode();assert hit.group().count(oldref)==1
        changed=changed[:hit.start()]+hit.group().replace(oldref,newref)+changed[hit.end():];refs.append((oldref,newref))
    assert len(clones)==2
    restored=changed
    for clone in clones:assert restored.count(clone)==1;restored=restored.replace(clone,b'',1)
    for old,new in refs:restored=restored.replace(new,old)
    assert restored==before
    (output/'candidate.osim').write_bytes(changed)
    parent['parent_cylinder_candidate_sha256']=parent['model_sha256']['candidate'];parent['model_sha256']['candidate']=hashlib.sha256(changed).hexdigest();parent['scope']='Four private corrected wrap clones, BRA cylinder and BIClong stationary ellipsoid; unchanged originals retained for all other94 muscle paths; exact reverse transformations preserve source bytes.'
    (output/'manifest.json').write_text(json.dumps(parent,indent=2)+'\n');return parent
