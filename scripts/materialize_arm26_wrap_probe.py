"""Copy frozen model; change only four wrap class tags for isolated observation."""
from pathlib import Path
import hashlib,json,re,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/derived/effective-potential-build-3a9juno_/inputs/subject_walk_scaled.osim'

def materialize(output):
    output=Path(output).resolve()
    if not output.is_relative_to(ROOT/'data/derived') or (output.exists() and any(output.iterdir())):raise ValueError('Fresh owned output required')
    raw=(ROOT/SOURCE).read_bytes();root=ET.fromstring(raw);changes=[]
    for muscle in root.findall('.//ForceSet/objects/*'):
        if muscle.get('name','').startswith(('arm26_BRA_','arm26_BIClong_')):
            for w in muscle.findall('.//PathWrap'):changes.append((w.findtext('wrap_object'),'Cylinder' if '_BRA_' in muscle.get('name') else 'Ellipsoid'))
    if len(changes)!=4:raise ValueError('Expected exact four wrap objects')
    output.mkdir(parents=True,exist_ok=True);models={}
    for mode in ('observed','candidate'):
        changed=raw;clones=[];references=[]
        for name,kind in changes:
            old='Wrap'+kind;new='ObservedArm26'+kind
            pattern=rb'<'+old.encode()+rb' name="'+re.escape(name.encode())+rb'">.*?</'+old.encode()+rb'>'
            hits=list(re.finditer(pattern,changed,re.S));assert len(hits)==1
            hit=hits[0];chunk=hit.group();replacement=chunk.replace(b'<'+old.encode()+b' ',b'<'+new.encode()+b' ',1).replace(b'</'+old.encode()+b'>',b'</'+new.encode()+b'>')
            if mode=='candidate' and kind=='Cylinder':
                clone_name=name+'__bra_exact'
                clone=replacement.replace(b'ObservedArm26Cylinder',b'ExactArm26Cylinder').replace(('name="'+name+'"').encode(),('name="'+clone_name+'"').encode(),1)
                clones.append(clone);replacement+=clone
                side=name[-1];muscle_name='arm26_BRA_'+side
                muscle_pattern=rb'<Thelen2003Muscle name="'+muscle_name.encode()+rb'">.*?</Thelen2003Muscle>'
                references.append((muscle_pattern,name,clone_name))
            changed=changed[:hit.start()]+replacement+changed[hit.end():]
        for pattern,name,clone_name in references:
            hits=list(re.finditer(pattern,changed,re.S));assert len(hits)==1
            hit=hits[0];chunk=hit.group();old_ref=('<wrap_object>'+name+'</wrap_object>').encode();new_ref=('<wrap_object>'+clone_name+'</wrap_object>').encode();assert chunk.count(old_ref)==1
            changed=changed[:hit.start()]+chunk.replace(old_ref,new_ref)+changed[hit.end():]
        restored=changed
        for clone in clones:
            assert restored.count(clone)==1;restored=restored.replace(clone,b'',1)
        for _,name,clone_name in references:restored=restored.replace(('<wrap_object>'+clone_name+'</wrap_object>').encode(),('<wrap_object>'+name+'</wrap_object>').encode())
        for old,new in [('ObservedArm26Cylinder','WrapCylinder'),('ExactArm26Cylinder','WrapCylinder'),('ObservedArm26Ellipsoid','WrapEllipsoid')]:restored=restored.replace(('<'+old+' ').encode(),('<'+new+' ').encode()).replace(('</'+old+'>').encode(),('</'+new+'>').encode())
        assert restored==raw
        (output/(mode+'.osim')).write_bytes(changed);models[mode]=hashlib.sha256(changed).hexdigest()
    r={'source':SOURCE,'source_sha256':hashlib.sha256(raw).hexdigest(),'model_sha256':models,'changed_wraps':changes,'native_enabled':False,'scope':'Observed class tags plus two BRA-private cylinder clones and two BRA-only reference changes; exact reversal recovers all original bytes; shared triceps wrap remains unchanged.'}
    (output/'manifest.json').write_text(json.dumps(r,indent=2)+'\n');return r
