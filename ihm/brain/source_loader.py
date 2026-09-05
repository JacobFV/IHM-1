"""Import IBM exclusively from verified source snapshots, bypassing bytecode."""
import importlib.abc
import importlib.util
import sys


class SnapshotLoader(importlib.abc.MetaPathFinder,importlib.abc.Loader):
    def __init__(self,source,snapshots,identity):
        self.source=source;self.snapshots=snapshots;self.identity=identity

    def find_spec(self,fullname,path=None,target=None):
        if fullname!='ibm' and not fullname.startswith('ibm.'):return None
        base=fullname.replace('.','/')
        candidates=[(base+'/__init__.py',True),(base+'.py',False)]
        for rel,package in candidates:
            if rel in self.snapshots:
                spec=importlib.util.spec_from_loader(fullname,self,is_package=package)
                spec.origin=str(self.source/rel);spec.loader_state={'relative_path':rel}
                if package:spec.submodule_search_locations=[str(self.source/base)]
                return spec
        raise ModuleNotFoundError('Module absent from pinned IBM source: '+fullname)

    def create_module(self,spec):return None

    def exec_module(self,module):
        rel=module.__spec__.loader_state['relative_path']
        module.__file__=str(self.source/rel)
        module.__ihm_source_identity__=self.identity
        # Compiling captured bytes also closes timestamp-valid .pyc and
        # post-verification filesystem races; no archive source is modified.
        exec(compile(self.snapshots[rel],module.__file__,'exec'),module.__dict__)


def install_snapshot(source,snapshots,identity):
    for name,module in tuple(sys.modules.items()):
        if name=='ibm' or name.startswith('ibm.'):
            if getattr(module,'__ihm_source_identity__',None)!=identity:
                raise RuntimeError('Another unverified IBM module is imported; use a fresh process')
    if not any(isinstance(f,SnapshotLoader) and f.identity==identity for f in sys.meta_path):
        sys.meta_path.insert(0,SnapshotLoader(source,snapshots,identity))
