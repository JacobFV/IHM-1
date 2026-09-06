"""Continuing native state with explicit post-preprocess mechanical boundaries."""
from .session import NativeSession
from . import number
from pathlib import Path


class CoupledNativeSession(NativeSession):
    executable_name='native_biogears_coupled'
    adapter_sources=('native_biogears_coupled.cpp','native_coupled_engine.h','native_body_ports.h','native_tissue_ports.h','native_tissue_compression.h')

    def launch_command(self,command):
        # Address-space ceiling bounds this owned process; it is not a host-wide
        # memory budget and does not constrain IBM or the visualization process.
        for path in ('/usr/bin/prlimit','/usr/bin/nice'):
            if not Path(path).is_file():raise RuntimeError('Required bounded native launcher unavailable: '+path)
        return ['/usr/bin/prlimit','--as=4294967296','--','/usr/bin/nice','-n','10',*command]

    def respiratory_load(self,external_pressure_pa):
        number(external_pressure_pa,-5000,5000,'external respiratory pressure Pa')
        return self._command('respiratory_load',external_pressure_pa)

    def skin_compression(self,external_pressure_pa):
        number(external_pressure_pa,0,5000,'whole-skin external compression Pa')
        result=self._command('skin_compression',external_pressure_pa)
        self._compression_installed=True
        return result

    def save_state(self):
        if getattr(self,'_compression_installed',False):
            raise ValueError('Native serialization of the added compression topology is not verified')
        if self._last['values'].get('coupling.external_pressure_pa',0)!=0:
            raise ValueError('Release external respiratory load before native save')
        return super().save_state()
