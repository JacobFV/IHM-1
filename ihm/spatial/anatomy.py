"""Auditable BP3D system crosswalk using recorded FMA relationships, never substrings.

The app system is a display grouping, not an FMA assertion. Is-a identity takes
precedence over contextual part-of membership; all supported groups are retained.
"""
import csv
import hashlib
import json
from collections import defaultdict, deque
from pathlib import Path

# Exact source concept labels are resolved to recorded FMA IDs at load time.
# This small curated crosswalk does not claim to reproduce all FMA semantics.
ANCHORS = {
    'skeletal': ['bone organ', 'cartilage organ', 'tooth', 'skeletal system'],
    'muscular': ['musculature', 'ascending part of trapezius', 'descending part of trapezius', 'muscle organ', 'zone of muscle organ', 'muscle tissue', 'head of muscle organ', 'set of muscles', 'skeletal muscle tissue'],
    'arterial': ['zone of artery', 'artery', 'set of arteries', 'segment of arterial tree organ', 'segment of artery', 'arterial anastomosis', 'variant artery', 'systemic arterial system'],
    'venous': ['vein', 'set of veins', 'venous network', 'segment of venous tree organ', 'venous anastomosis', 'venous tree organ', 'systemic venous system', 'portal venous system'],
    'cardiac': ['heart', 'leaf of cardiac valve', 'cusp of cardiac valve', 'cavity of cardiac chamber', 'region of wall of heart', 'region of papillary muscle', 'region of papillary muscle of left ventricle', 'papillary muscle of right ventricle', 'papillary muscle of left ventricle'],
    'connective': ['linea alba', 'common tendinous ring', 'tarsal plate of eyelid', 'pterygomandibular raphe', 'pharyngeal raphe', 'thyrohyoid membrane', 'conus elasticus', 'ligament organ', 'ligament organ component', 'interosseous membrane', 'retinaculum', 'raphe', 'tendon', 'deep fascial system', 'zone of investing fascia'],
    'nervous': ['cell part cluster of neuraxis', 'gray matter of neuraxis', 'white matter of neuraxis', 'choroid plexus', 'lamina terminalis', 'septum of telencephalon', 'dura mater', 'region of dura mater', 'interventricular foramen', 'subdivision of subarachnoid space', 'choroid plexus of cerebral hemisphere', 'nervous system', 'nerve', 'ganglion', 'organ component of neuraxis', 'segment of neural tree organ', 'neural tree organ', 'subdivision of nervous system', 'brain', 'spinal cord', 'segment of brain', 'nucleus of brain', 'ventricular system of neuraxis'],
    'respiratory': ['respiratory system'],
    'urinary': ['urinary system'],
    'reproductive': ['genital system', 'testis', 'epididymis', 'deferent duct', 'seminal vesicle', 'penis', 'corpus cavernosum of penis', 'corpus spongiosum of penis', 'glans penis'],
    'endocrine': ['endocrine system', 'adrenal gland', 'pituitary gland', 'pineal body', 'thyroid gland', 'parathyroid gland'],
    'digestive': ['segment of liver', 'hepatovenous subsector', 'segment of biliary tree', 'subdivision of mouth', 'subdivision of pharynx', 'hepatovenous segment', 'duct of caudate lobe of liver', 'alimentary system', 'pancreaticobiliary system', 'salivary gland', 'bile duct'],
    'integumentary': ['eyebrow', 'integumentary system', 'skin', 'nail', 'hair'],
    'sensory': ['right eye', 'left eye', 'external ear', 'region of layer of wall of eyeball', 'layer of wall of eyeball', 'chamber of eyeball', 'segment of lacrimal duct', 'lacrimal lake', 'optic part of retina', 'lacrimal gland'],
    'lymphatic': ['spleen', 'thymus', 'lymph node', 'lymphatic vessel'],
}


class AnatomyOntology:
    def __init__(self, raw_dir='data/raw/anatomy/bodyparts3d'):
        self.raw_dir = Path(raw_dir)
        self.names = {}
        self.parents = {kind:defaultdict(set) for kind in ('isa', 'partof')}
        self.elements = {kind:defaultdict(set) for kind in ('isa', 'partof')}
        self.ontology_sha256 = {}
        for kind in self.parents:
            for suffix in ('parts_list_e', 'inclusion_relation_list', 'element_parts'):
                path = self.raw_dir / f'{kind}_{suffix}.txt'
                self.ontology_sha256[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
                with path.open() as handle:
                    rows = list(csv.reader(handle, delimiter='\t'))[1:]
                for row in rows:
                    if suffix == 'parts_list_e':
                        self.names[row[0]] = row[2]
                    elif suffix == 'inclusion_relation_list':
                        parent, pname, child, cname = row
                        self.names.update({parent:pname, child:cname})
                        self.parents[kind][child].add(parent)
                    else:
                        concept, name, element = row
                        self.names[concept] = name
                        self.elements[kind][element].add(concept)
        self.ids_by_name = {name:cid for cid, name in self.names.items()}
        self.anchors = {self.ids_by_name[name]:system for system,names in ANCHORS.items()
                        for name in names if name in self.ids_by_name}
        self.crosswalk_sha256 = hashlib.sha256(json.dumps(ANCHORS, sort_keys=True).encode()).hexdigest()
        self._cache = {}
        self.concept_frequency = defaultdict(int)
        for ids in self.elements['isa'].values():
            for cid in ids:
                self.concept_frequency[cid] += 1

    def paths(self, concept, kind):
        key = concept, kind
        if key not in self._cache:
            queue = deque([(concept, [concept])]); seen = set(); paths = []
            while queue:
                current, path = queue.popleft()
                if current in seen:
                    continue
                seen.add(current)
                if current in self.anchors:
                    paths.append({'system':self.anchors[current], 'relation':kind,
                                  'concept_id':concept, 'anchor_id':current,
                                  'path':path, 'path_names':[self.names.get(c,c) for c in path]})
                queue.extend((parent, path+[parent]) for parent in sorted(self.parents[kind][current]))
            self._cache[key] = paths
        return self._cache[key]

    def classify(self, concepts, element_id=None):
        supplied = {c['concept_id'] for c in concepts}
        by_kind = {kind: supplied | self.elements[kind].get(element_id,set()) for kind in self.parents}
        # Source tables include broad ancestor memberships. Start at the most
        # specific recorded concepts so the evidence includes meaningful paths.
        evidence = []
        preferred_name = None
        for kind, ids in by_kind.items():
            ancestors = set()
            queue = list(ids)
            while queue:
                cid = queue.pop()
                for parent in self.parents[kind][cid]:
                    if parent not in ancestors:
                        ancestors.add(parent); queue.append(parent)
            leaves = ids - ancestors
            if kind == 'isa' and leaves:
                ranked = sorted(leaves, key=lambda cid:(self.concept_frequency[cid], cid))
                if len(ranked)==1 or self.concept_frequency[ranked[0]] < self.concept_frequency[ranked[1]]:
                    preferred_name = self.names.get(ranked[0])
            for cid in sorted(leaves):
                evidence.extend(self.paths(cid, kind))
        # Intrinsic tissue/vessel identity prevents e.g. thyroid artery from
        # becoming endocrine or an organ's supporting bone becoming muscle.
        intrinsic = [e for e in evidence if e['relation'] == 'isa']
        candidates = intrinsic or evidence
        candidates = sorted(candidates, key=lambda e:(len(e['path']), e['system'], e['concept_id']))
        primary = candidates[0]['system'] if candidates else 'other'
        # Papillary myocardium has a direct source cardiac anchor even when
        # contextual muscle memberships also exist.
        if any(e['system']=='cardiac' for e in intrinsic):
            primary = 'cardiac'
        systems = sorted({e['system'] for e in evidence})
        best = {}
        for e in sorted(evidence, key=lambda e:(len(e['path']), e['relation'], e['concept_id'])):
            best.setdefault((e['system'],e['relation']),e)
        return {'preferred_name':preferred_name, 'system':primary, 'systems':systems or ['other'],
                'status':'classified' if candidates else 'unresolved',
                'reason':'Recorded FMA is-a identity, then part-of context; exact concept anchors.' if candidates else 'No recorded path to the curated system crosswalk.',
                'evidence':list(best.values()),
                'classification_basis':'isa_identity' if intrinsic else 'partof_context' if evidence else 'unresolved',
                'intrinsic_systems':sorted({e['system'] for e in intrinsic}),
                'contextual_systems':sorted({e['system'] for e in evidence if e['relation']=='partof'}),
                'ambiguous_primary':len({e['system'] for e in candidates if len(e['path'])==len(candidates[0]['path'])})>1 if candidates else False,
                'source_concept_ids':sorted(set.union(*by_kind.values())),
                'ontology_sha256':self.ontology_sha256, 'crosswalk_sha256':self.crosswalk_sha256}
