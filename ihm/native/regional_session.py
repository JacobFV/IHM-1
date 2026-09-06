"""Signed native physiology with exclusively owned regional Skin compartments.

The fixed 0.2/0.3/0.5 allocation is an engineering materialization. Region IDs
have no implied anatomical territory until an explicit registration supplies it.
"""
from . import number
from .coupled_session import CoupledNativeSession,SignedCoupledNativeSession


class RegionalSignedNativeSession(SignedCoupledNativeSession):
    executable_name='native_biogears_regional_signed'
    adapter_sources=('native_biogears_regional_signed.cpp','native_regional_coupled_engine.h',
        'native_body_ports.h','native_tissue_ports.h','native_regional_skin.h',
        'native_regional_species.h','native_signed_muscle_port.h')
    regions=('region_a','region_b','residual')

    def __init__(self,config,output_dir,timeout_s=1800):
        if config.engine_variant!='whole_body_integrity_regional_skin_graph_v2':
            raise ValueError('Regional signed adapter requires its verified regional source variant')
        if getattr(config,'state_path',None) is None:raise ValueError('Regional installation requires a retained native state')
        CoupledNativeSession.__init__(self,config,output_dir,timeout_s)

    def regional_skin_pressure(self,region,pressure_pa):
        if not isinstance(region,str) or region not in self.regions:raise ValueError('Unknown native Skin region')
        number(pressure_pa,0,5000,'regional external skin pressure Pa')
        return self._command('regional_skin_pressure',region,pressure_pa)

    def skin_compression(self,external_pressure_pa):
        raise ValueError('Whole-Skin compression topology is incompatible; specify regional boundaries')

    def save_state(self):
        raise ValueError('Regional topology and signed boundary serialization is not verified')
