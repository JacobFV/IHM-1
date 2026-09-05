"""Fail-closed regression tests; all mutations use temporary synthetic fixtures."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import build_extended_anatomy as builder


class SourceIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        root=Path(self.temp.name);raw=root/'raw';out=root/'out';raw.mkdir();out.mkdir()
        self.root,self.raw,self.out=root,raw,out
        self.patchers=[patch.object(builder,k,v) for k,v in [('ROOT',root),('RAW',raw),('OUT',out),('PINNED_REVISION','fixture')]]
        for p in self.patchers:p.start();self.addCleanup(p.stop)
        records=[];expected={}
        for name in ['archive.zip','README','LICENSE']:
            path=raw/name;path.write_bytes(name.encode());digest=builder.sha256(path);expected[name]=digest
            records.append({'path':str(path.relative_to(root)),'sha256':digest,'bytes':path.stat().st_size})
        blend=raw/'fixture.blend';blend.write_bytes(b'fixture blend')
        for key,value in [('PINNED_FILES',expected),('PINNED_BLEND_SHA256',builder.sha256(blend))]:
            p=patch.object(builder,key,value);p.start();self.addCleanup(p.stop)
        self.provenance={'revision':'fixture','files':records,'blend_path':'raw/fixture.blend','blend_sha256':builder.sha256(blend)}
        (raw/'provenance.json').write_text(json.dumps(self.provenance))
        for name in ['source','base']:
            np.savez(root/(name+'.npz'),vertices=np.eye(3),faces=np.array([[0,1,2]]))
        self.entry={'id':'fixture','source_geometry_path':'source.npz','source_geometry_sha256':builder.sha256(root/'source.npz'),'base_geometry_path':'base.npz','base_geometry_sha256':builder.sha256(root/'base.npz')}
        self.index={'schema_version':2,'extraction':builder.EXTRACTION_DESCRIPTION,'source_blend_sha256':self.provenance['blend_sha256'],'meshes':[self.entry]}
        self.save_index()

    def save_index(self):
        (self.out/'source_index.json').write_text(json.dumps(self.index))

    def test_legacy_semantic_v1_cache_is_accepted(self):
        self.assertEqual(builder.validate_cached_sources(builder.validate_provenance())['meshes'],[self.entry])

    def test_extraction_rejects_different_loaded_document(self):
        with patch.dict('sys.modules',{'bpy':SimpleNamespace(data=SimpleNamespace(filepath=str(self.raw/'different.blend')))}):
            with self.assertRaisesRegex(ValueError,'Loaded Blender document'):builder.extract()
        self.assertFalse((self.out/'source_geometry').exists())

    def test_changed_source_is_rejected_before_build_outputs(self):
        (self.root/'source.npz').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'source geometry'):builder.build()
        self.assertFalse((self.out/'geometry').exists())

    def test_changed_base_is_rejected(self):
        (self.root/'base.npz').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'base geometry'):builder.validate_cached_sources(builder.validate_provenance())

    def test_missing_blend_is_rejected(self):
        (self.raw/'fixture.blend').unlink()
        with self.assertRaisesRegex(ValueError,'blend'):builder.validate_provenance()

    def test_changed_raw_record_is_rejected(self):
        (self.raw/'README').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'raw source'):builder.validate_provenance()

    def test_incompatible_extractor_is_rejected(self):
        self.index['extractor_semantic_version']=999;self.save_index()
        with self.assertRaisesRegex(ValueError,'extractor'):builder.validate_cached_sources(builder.validate_provenance())

    def test_unknown_legacy_pipeline_is_rejected(self):
        self.index['extraction']='decimated meshes';self.save_index()
        with self.assertRaisesRegex(ValueError,'extractor'):builder.validate_cached_sources(builder.validate_provenance())

    def test_append_checks_all_geometry_before_copy(self):
        geom=self.out/'mesh.json.gz';geom.write_bytes(b'geometry')
        s={'id':'fixture','geometry_path':'out/mesh.json.gz','geometry_sha256':builder.sha256(geom),'source':{'sha256':self.provenance['blend_sha256'],'geometry_sha256':self.entry['source_geometry_sha256'],'base_geometry_sha256':self.entry['base_geometry_sha256']}}
        (self.out/'manifest_fragment.json').write_text(json.dumps({'structures':[s]}))
        geom.write_bytes(b'changed')
        target=self.root/'app';target.mkdir();manifest=target/'manifest.json';original='{"models":[],"structures":[]}'
        manifest.write_text(original)
        with self.assertRaisesRegex(ValueError,'display geometry'):builder.append_manifest(manifest)
        self.assertEqual(manifest.read_text(),original);self.assertFalse((target/'geometry').exists())


if __name__=='__main__':unittest.main()
