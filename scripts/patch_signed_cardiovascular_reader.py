"""Prepare isolated reader/vascular-prior source overlays; never build/promote."""
from pathlib import Path
import argparse,difflib,hashlib,json,tempfile
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/runtime/physiology/variants/whole_body_integrity_signed_muscle_v2/Cardiovascular.cpp'
SOURCE_SHA='ad93e3930f307f487e166dac1795234074f73062fc11b44eead288fe5d305f6d'
ANCHOR='    MetabolicToneResponse();\n  }\n}\n//--------------------------------------------------------------------------------------------------'
READER='''    MetabolicToneResponse();
  } else if (ihm_signed::active() && m_data.GetEnergy().HasTotalMetabolicRate()) {
    // Routine signed-demand observation only. No MAP/resistance/energy write.
    (void)ihm_signed::effective(m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::kcal_Per_day), Convert(1.0, PowerUnit::W, PowerUnit::kcal_Per_day), 0);
  }
'''
PRIOR='''  // EXPERIMENTAL: normalized transfer of the inherited exercise resistance prior.
  // Existing autonomic baseline reset must own every affected path each step;
  // otherwise this multiplier would compound rather than release at zero demand.
  if (auto* signedMuscle = ihm_signed::active()) {
    if (!m_data.GetNervous().HasResistanceScaleMuscle() || !m_data.GetNervous().HasResistanceScaleExtrasplanchnic() || !m_data.GetNervous().HasResistanceScaleSplanchnic())
      throw std::runtime_error("signed vascular prior requires native regional baseline resets");
    const auto factors = ihm_signed_vascular::response(m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::W), m_patient->GetBasalMetabolicRate(PowerUnit::W), signedMuscle->delta_m_W);
    auto* trace = ihm_signed_vascular::ihm_signed_vascular_trace();
    *trace = {double(signedMuscle->sequence), m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::W), m_patient->GetBasalMetabolicRate(PowerUnit::W), signedMuscle->delta_m_W, factors.ratio};
    for (auto* path : m_muscleResistancePaths) trace->muscle_before += path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
    for (auto* path : m_extrasplanchnicResistancePaths) trace->other_before += path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
    for (auto* path : m_splanchnicResistancePaths) trace->other_before += path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
    if (signedMuscle->delta_m_W != 0) {
      auto applyPrior = [&](const auto& paths, double factor) {
        for (SEFluidCircuitPath* path : paths) {
          if (!path->HasNextResistance()) throw std::runtime_error("signed vascular prior missing native resistance");
          const double prior = path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
          const double updated = prior * factor;
          if (!std::isfinite(updated) || prior < 0) throw std::runtime_error("signed vascular prior invalid native resistance");
          path->GetNextResistance().SetValue(std::max(updated, m_minIndividialSystemicResistance__mmHg_s_Per_mL), FlowResistanceUnit::mmHg_s_Per_mL);
        }
      };
      applyPrior(m_muscleResistancePaths, factors.muscle);
      applyPrior(m_extrasplanchnicResistancePaths, factors.other);
      applyPrior(m_splanchnicResistancePaths, factors.other);
    }
    for (auto* path : m_muscleResistancePaths) trace->muscle_after += path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
    for (auto* path : m_extrasplanchnicResistancePaths) trace->other_after += path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
    for (auto* path : m_splanchnicResistancePaths) trace->other_after += path->GetNextResistance(FlowResistanceUnit::mmHg_s_Per_mL);
  }
'''
def patch(source,*,vascular=False):
    if hashlib.sha256(source.encode()).hexdigest()!=SOURCE_SHA or source.count(ANCHOR)!=1:raise ValueError('Changed or ambiguous pinned cardiovascular source')
    replacement=READER+(PRIOR if vascular else '')+'}\n//--------------------------------------------------------------------------------------------------'
    result=source.replace(ANCHOR,replacement,1)
    if vascular:result='#include "native_signed_vascular_prior.h"\nnamespace ihm_signed_vascular { extern "C" Trace* ihm_signed_vascular_trace() { static thread_local Trace trace; return &trace; } }\n'+result
    return result

def prepare(output=None):
    source=SOURCE.read_text();out=Path(output) if output else Path(tempfile.mkdtemp(prefix='signed-cardiovascular-overlay-',dir=ROOT/'data/derived/audits'))
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Fresh empty overlay output required')
    receipts={}
    for name,vascular in [('reader_only',False),('vascular_prior',True)]:
        directory=out/name;directory.mkdir();after=patch(source,vascular=vascular);target=directory/'Cardiovascular.cpp';target.write_text(after)
        (directory/'native_signed_muscle_port.h').write_bytes((SOURCE.parent/'native_signed_muscle_port.h').read_bytes())
        if vascular:(directory/'native_signed_vascular_prior.h').write_bytes((ROOT/'scripts/native_signed_vascular_prior.h').read_bytes())
        (directory/'Cardiovascular.patch').write_text(''.join(difflib.unified_diff(source.splitlines(True),after.splitlines(True),fromfile=str(SOURCE),tofile=str(target))))
        receipts[name]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir()}
    report={'schema':'ihm.signed-cardiovascular-source-overlay.v1','built':False,'source':str(SOURCE.relative_to(ROOT)),'source_sha256':SOURCE_SHA,'files':receipts,'scope':'Separate measurement-only reader correction and uncalibrated transferred exercise resistance prior; no production source/library mutation'}
    (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args();print(prepare(args.output))
