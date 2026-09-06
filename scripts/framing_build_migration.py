"""Audit one exact protocol-ordering change without weakening physics identity."""
import hashlib,re
from pathlib import Path

OLD='if(command=="evaluate_static_pose"){std::cout<<"@IHM "<<ihm_static_pose::evaluate(model,state,in,environment,!external->loads.empty(),support_plane,surface_foundation)<<std::endl;continue;}'
NEW='if(command=="evaluate_static_pose"){const auto payload=ihm_static_pose::evaluate(model,state,in,environment,!external->loads.empty(),support_plane,surface_foundation);std::cout<<"@IHM "<<payload<<std::endl;continue;}'


def audit(previous,current,root):
    root=Path(root);old=previous['build_files'];new=current['build_files']
    def retained_source(files):
        found=[p for p in files if p.startswith('data/runtime/mechanical-stream/build-') and p.endswith('/native_mechanical_stream.cpp')]
        if len(found)!=1:raise ValueError('Missing unique retained native source')
        path=found[0];content=(root/path).read_bytes()
        if hashlib.sha256(content).hexdigest()!=files[path]:raise ValueError('Retained source hash mismatch')
        return content.decode()
    before=retained_source(old);after=retained_source(new)
    if before.count(OLD)!=1 or before.replace(OLD,NEW)!=after:raise ValueError('Native diff is not exactly audited framing order')
    def physical(files):
        result={}
        for path,expected in files.items():
            if path=='scripts/native_mechanical_stream.cpp':continue
            data=(root/path).read_bytes()
            if hashlib.sha256(data).hexdigest()!=expected:raise ValueError('Retained native dependency hash mismatch')
            if path.endswith('/native_mechanical_stream.cpp') or path.endswith('/native_mechanical_stream'):continue
            normalized=re.sub(r'data/runtime/mechanical-stream/build-[^/]+/', 'retained-build/',path)
            result[normalized]=expected
        return result
    if physical(old)!=physical(new):raise ValueError('Physical native headers/libraries changed')
    return dict(kind='exact-static-dispatch-framing-order-only',old_source_sha256=hashlib.sha256(before.encode()).hexdigest(),
        new_source_sha256=hashlib.sha256(after.encode()).hexdigest(),physical_dependency_count=len(physical(old)))
