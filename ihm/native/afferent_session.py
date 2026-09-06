"""Opt-in observation-only native nervous adapters; receptor state stays native."""
from .coupled_session import SignedCoupledNativeSession
from .regional_session import RegionalSignedNativeSession

class AfferentSignedNativeSession(SignedCoupledNativeSession):
    executable_name='native_biogears_afferent_signed'
    adapter_sources=SignedCoupledNativeSession.adapter_sources+('native_nervous_afferents.h','native_biogears_afferent_signed.cpp')

class AfferentRegionalNativeSession(RegionalSignedNativeSession):
    executable_name='native_biogears_afferent_regional'
    adapter_sources=RegionalSignedNativeSession.adapter_sources+('native_nervous_afferents.h','native_biogears_afferent_regional.cpp')
