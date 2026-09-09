"""Bounded immutable bytes registered by native owners, never by HTTP paths."""
import gzip
import hashlib
import re
import threading

DIGEST = re.compile(r'[0-9a-f]{64}\Z')


class SurfaceAssetRegistry:
    def __init__(self, *, maximum_asset_bytes=64*1024*1024, maximum_total_bytes=256*1024*1024):
        self.maximum_asset_bytes=maximum_asset_bytes
        self.maximum_total_bytes=maximum_total_bytes
        self._records={};self._bytes=0;self._lock=threading.Lock()

    def register(self, records):
        if not isinstance(records,dict):raise ValueError('Surface asset records must be a digest-to-bytes mapping')
        staged={}
        for digest,raw in records.items():
            if not isinstance(digest,str) or not DIGEST.fullmatch(digest):raise ValueError('Invalid surface asset digest')
            if not isinstance(raw,bytes) or not raw or len(raw)>self.maximum_asset_bytes:raise ValueError('Surface asset exceeds immutable byte envelope')
            if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('Surface asset integrity mismatch')
            staged[digest]=(raw,gzip.compress(raw,compresslevel=1,mtime=0))
        with self._lock:
            fresh={k:v for k,v in staged.items() if k not in self._records}
            for digest,value in staged.items():
                if digest in self._records and self._records[digest][0]!=value[0]:raise ValueError('Surface asset identity collision')
            additional=sum(len(raw)+len(compressed) for raw,compressed in fresh.values())
            if self._bytes+additional>self.maximum_total_bytes:raise ValueError('Surface asset registry capacity exceeded')
            self._records.update(fresh);self._bytes+=additional

    def get(self,digest,*,compressed=False):
        if not isinstance(digest,str) or not DIGEST.fullmatch(digest):raise ValueError('Invalid surface asset digest')
        with self._lock:
            value=self._records.get(digest)
            if value is None:raise KeyError('Surface asset is not retained by this server')
            return value[1 if compressed else 0]


def register_body_surface_assets(body,registry):
    binding=getattr(getattr(body,'plant',None),'surface_binding',None)
    if binding is not None:
        records=getattr(binding,'asset_records',None)
        if records is not None:registry.register(records())
