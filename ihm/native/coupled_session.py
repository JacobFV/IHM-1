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


class SignedCoupledNativeSession(CoupledNativeSession):
    """Muscle-only chemical/heat/work increments on the audited native variant."""
    executable_name='native_biogears_signed'
    adapter_sources=('native_biogears_signed.cpp','native_coupled_engine.h','native_body_ports.h',
        'native_tissue_ports.h','native_tissue_compression.h','native_signed_muscle_port.h','native_intake_receipts.h')

    def __init__(self,config,output_dir,timeout_s=1800):
        if config.engine_variant not in ('whole_body_integrity_signed_muscle_v2','whole_body_integrity_gi_absorption'):raise ValueError('Signed adapter requires the matching source variant')
        super().__init__(config,output_dir,timeout_s)

    def step(self,seconds):
        if hasattr(self,'_metabolic_reference'):
            raise ValueError('Use signed_step for every interval after metabolic reference binding')
        return super().step(seconds)

    def signed_step(self,reference_id,delta_m_w,delta_h_w,delta_w_w):
        import re
        if not isinstance(reference_id,str) or not re.fullmatch('[a-f0-9]{64}',reference_id):raise ValueError('SHA256 metabolic reference required')
        for value in (delta_m_w,delta_h_w,delta_w_w):number(value,-5000,5000,'signed muscle increment W')
        if abs(delta_m_w-delta_h_w-delta_w_w)>1e-10*(1+abs(delta_m_w)+abs(delta_h_w)+abs(delta_w_w)):raise ValueError('Chemical/heat/work incidence mismatch')
        if self._elapsed_ticks>=round(self.config.horizon_s*50):raise ValueError('Native horizon reached')
        if getattr(self,'_metabolic_reference',reference_id)!=reference_id:raise ValueError('Metabolic reference cannot change during a session')
        try:result=self._command('signed_step',reference_id,delta_m_w,delta_h_w,delta_w_w,advance_ticks=1)
        except EOFError as error:
            log=self.out/'runner_stdout.log'
            if log.is_file():
                with log.open('rb') as stream:
                    stream.seek(max(0,log.stat().st_size-8192));tail=stream.read().decode(errors='replace')
                diagnostics=[line for line in tail.splitlines() if line.startswith('SIGNED_NATIVE_FAILURE:')]
                if diagnostics:raise RuntimeError(diagnostics[-1][:2048]) from error
            raise
        self._metabolic_reference=reference_id
        return result

    def save_state(self):
        if hasattr(self,'_metabolic_reference'):raise ValueError('Signed coupling reference/ledger serialization is not verified')
        return super().save_state()
