#!/usr/bin/env python3
"""Explicit bounded fake/native tests for observational consumed-intake receipts."""
from pathlib import Path
import argparse,json,subprocess,tempfile,hashlib,time
BASE=Path(__file__).resolve().parents[1]
FAKE=r'''
#include "native_intake_receipts.h"
#include <cassert>
#include <iostream>
using namespace ihm_intake;
struct Nutrition {
 double weight=.013;std::string GetName()const{return "actual_held";}
 bool HasCarbohydrate()const{return true;}bool HasProtein()const{return true;}
 bool HasFat()const{return true;}bool HasSodium()const{return true;}
 bool HasCalcium()const{return true;}bool HasWater()const{return true;}
 double GetWeight(int)const{return weight;}double GetCarbohydrate(int)const{return 1;}
 double GetProtein(int)const{return 1;}double GetFat(int)const{return 1;}
 double GetSodium(int)const{return 0;}double GetCalcium(int)const{return 0;}double GetWater(int)const{return 10;}
};
template<class F>void reject(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}assert(caught);}
void emit(const char* name,const Ledger& ledger){std::cout<<name<<'\t';ledger.write_json(std::cout);std::cout<<'\n';}
int main(){
 Nutrition n;const Payload p=capture(n,0,1,2,3);n.weight=.5;assert(capture(n,0,1,2,3).mass_kg==.5);
 const std::uint64_t meal=9007199254740999ULL;
 Ledger ledger("fixture-owner");emit("initial",ledger);ledger.accepted(meal,p);emit("queued",ledger);
 ledger.before(p,meal+1,0,3600);ledger.commit(false,1,3600.02);emit("consumed",ledger);
 ledger.before(std::nullopt,meal+2,1,3600.02);ledger.commit(false,2,3600.04);emit("later",ledger);
 ledger.accepted(meal+3,p);ledger.before(p,meal+4,2,3600.04);ledger.commit(false,3,3600.06);emit("second",ledger);
 Ledger duplicate("duplicate");duplicate.accepted(1,p);reject([&]{duplicate.accepted(1,p);});
 Ledger missing("missing");reject([&]{missing.before(p,2,0,0);});
 Ledger changed("changed");changed.accepted(1,p);auto q=p;q.nutrients[5]+=1;reject([&]{changed.before(q,2,0,0);});
 Ledger transition("transition");transition.accepted(1,p);transition.before(p,2,0,0);reject([&]{transition.commit(true,1,.02);});
 reject([&]{transition.write_json(std::cout);});
 Ledger failed("failed");failed.accepted(1,p);failed.before(p,2,0,0);failed.terminal_failure();
 reject([&]{failed.commit(false,1,.02);});reject([&]{failed.write_json(std::cout);});
 Ledger foreign("foreign");foreign.before(std::nullopt,1,0,0);reject([&]{foreign.commit(true,1,.02);});
 Ledger repeated("repeated");repeated.accepted(1,p);repeated.before(p,2,0,0);repeated.commit(false,1,.02);reject([&]{repeated.commit(false,1,.02);});
}
'''
def fake():
    from ihm.assembly.intake_mass import consumed_intake_delta
    out=Path(tempfile.mkdtemp(prefix='native-intake-fake-',dir=BASE/'data/derived/audits'))
    source=out/'fixture.cpp';source.write_text(FAKE);binary=out/'fixture'
    command=['c++','-std=c++20','-Wall','-Wextra','-Werror','-I'+str(BASE/'scripts'),str(source),'-o',str(binary)]
    subprocess.run(command,check=True,timeout=30)
    result=subprocess.run([str(binary)],text=True,capture_output=True,check=True,timeout=5)
    (out/'stdout.log').write_text(result.stdout)
    rows={key:json.loads(raw) for key,raw in (line.split('\t',1) for line in result.stdout.splitlines())}
    assert consumed_intake_delta(rows['initial'],rows['queued']) is None
    first=consumed_intake_delta(rows['queued'],rows['consumed']);assert first['mass_kg']==.013
    assert rows['consumed']['last_consumed']['meal_sequence']==9007199254740999
    assert consumed_intake_delta(rows['consumed'],rows['later']) is None
    assert consumed_intake_delta(rows['later'],rows['second'])['mass_kg']==.013
    (out/'verification.json').write_text(json.dumps({'passed':True,'fixture_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'header_sha256':hashlib.sha256((BASE/'scripts/native_intake_receipts.h').read_bytes()).hexdigest(),'command':command,'receipts':rows,
        'scope':'Fake native quantity interface + real ledger state machine/parser; no physiology'},indent=2)+'\n')
    print(json.dumps({'passed':True,'kind':'fake','output':str(out)}))

def native():
    from ihm.native.session import SessionConfig,Meal
    from ihm.native.coupled_session import SignedCoupledNativeSession
    from ihm.native.regional_session import RegionalSignedNativeSession
    from ihm.assembly.intake_mass import consumed_intake_delta
    from build_signed_coupled_native import intake_source_receipt,RUNTIME
    out=Path(tempfile.mkdtemp(prefix='native-intake-boundary-',dir=BASE/'data/derived/audits'));started=time.monotonic()
    state_receipt=BASE/'data/derived/systemic/exertion_v3/exercise/native/manifest.json'
    state=Path(json.loads(state_receipt.read_text())['configuration']['state_path'])
    ref=hashlib.sha256(b'intake adapter boundary fixture only').hexdigest();cases={};session=None;epochs=[]
    try:
        for label,kind,variant in [('signed',SignedCoupledNativeSession,'whole_body_integrity_gi_absorption'),
                                   ('regional',RegionalSignedNativeSession,'whole_body_integrity_regional_skin_graph_v2')]:
            config=SessionConfig(state_path=state,engine_variant=variant,horizon_s=.08)
            session=kind(config,out/label);initial=session.snapshot();epochs.append(initial['intake']['owner_epoch'])
            queued=session.meal(Meal(name='native_intake_one',carbohydrate_g=1,protein_g=2,fat_g=3,sodium_g=.1,calcium_mg=100,water_ml=10))
            assert queued['pending_meal'] is True and queued['elapsed_s']==0
            assert queued['intake']['pending']['meal_sequence']==queued['sequence']
            assert consumed_intake_delta(initial['intake'],queued['intake']) is None
            consumed=session.step(.04) # Two native calls: ingestion once, then internal GI redistribution.
            first=consumed_intake_delta(queued['intake'],consumed['intake']);assert abs(first['mass_kg']-.0162)<1e-14
            assert first['interval_start_tick']==0 and first['interval_end_tick']==1
            assert consumed['pending_meal'] is False and consumed['intake']['consumed_count']==1
            snapshot=session.snapshot();assert snapshot['intake']==consumed['intake']
            queued2=session.meal(Meal(name='native_intake_two',water_ml=20))
            assert consumed_intake_delta(snapshot['intake'],queued2['intake']) is None
            consumed2=session.signed_step(ref,0,0,0)
            second=consumed_intake_delta(queued2['intake'],consumed2['intake']);assert abs(second['mass_kg']-.02)<1e-14
            assert second['interval_start_tick']==2 and second['interval_end_tick']==3
            assert second['native_advance_sequence']==consumed2['sequence']
            continued=session.signed_step(ref,0,0,0);assert consumed_intake_delta(consumed2['intake'],continued['intake']) is None
            assert continued['intake']['consumed_count']==2 and abs(continued['intake']['cumulative_mass_kg']-.0362)<1e-14
            session.close();session=None
            cases[label]={'frames':[initial,queued,consumed,snapshot,queued2,consumed2,continued],'deltas':[first,second],
                'intake_consumer_guard':intake_source_receipt(RUNTIME/'variants'/variant),
                'diagnostic_patient_weight_delta_kg':consumed['values']['patient_weight_kg']-initial['values']['patient_weight_kg']}
            # GI can consume before a later Energy/Tissue failure. No successful
            # receipt is permitted and that uncertain owner must never be retried.
            session=kind(config,out/(label+'-failed'));before=session.snapshot();epochs.append(before['intake']['owner_epoch'])
            pending=session.meal(Meal(name='terminal_intake',water_ml=10))
            try:session.signed_step(ref,-100,-100,0)
            except (RuntimeError,EOFError):pass
            else:raise AssertionError('Expected native signed thermal-domain failure')
            assert session._closed and session.process.poll() is not None
            assert session._last['intake']['consumed_count']==0
            log=(out/(label+'-failed')/'runner_stdout.log').read_text()
            attempts=[json.loads(line.split(' ',1)[1]) for line in log.splitlines() if line.startswith('INTAKE_ATTEMPT ')]
            assert attempts and attempts[-1]['pending']['meal_sequence']==pending['sequence']
            assert attempts[-1]['pending']['native_payload']['water_ml']==10
            cases[label]['failed_attempt']=attempts[-1];cases[label]['failed_process_exit_code']=session.process.returncode
            session=None
        assert len(set(epochs))==len(epochs)
        record={'passed':True,'wall_s':time.monotonic()-started,'cases':cases,'owner_epochs':epochs,
            'scope':'Native queued/consumed distinction, exact enqueue identity, once-only multi-step and signed credit, epoch reset, terminal failure; no mechanical mass application or all-boundary body-mass conservation'}
        (out/'verification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'passed':True,'kind':'native','output':str(out),'wall_s':record['wall_s']}))
    except BaseException as error:
        (out/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-started},indent=2)+'\n');raise
    finally:
        if session:session.close(graceful=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--fake',action='store_true');parser.add_argument('--native',action='store_true');args=parser.parse_args()
    if not (args.fake or args.native):raise SystemExit('Explicit --fake or --native and coordinated resource slot required')
    if args.fake:fake()
    if args.native:native()
