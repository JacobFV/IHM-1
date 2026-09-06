"""Compile the actual dispatch fragment with a logging/throwing stub evaluator."""
from pathlib import Path
import json,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[1]

def main():
    source=(ROOT/'scripts/native_mechanical_stream.cpp').read_text()
    line=next(line.strip() for line in source.splitlines() if 'if(command=="evaluate_static_pose")' in line)
    expression='ihm_static_pose::evaluate(model,state,in,environment,!external->loads.empty(),support_plane,surface_foundation)'
    assert expression in line;line=line.replace(expression,'evaluate(fail)')
    code='''#include <iostream>
#include <stdexcept>
#include <string>
std::string evaluate(bool fail){std::cout<<"[warn] evaluator diagnostic\\n";if(fail)throw std::runtime_error("physical domain failure");return "{\\"kind\\":\\"static_pose_evaluated\\"}";}
int main(int argc,char**){bool fail=argc>1;std::string command="evaluate_static_pose";for(int i=0;i<1;++i){try{'''+line+'''}catch(const std::exception& e){std::cout<<"@IHM {\\"error\\":\\""<<e.what()<<"\\"}\\n";}}}
'''
    with tempfile.TemporaryDirectory() as temporary:
        folder=Path(temporary);(folder/'fixture.cpp').write_text(code)
        subprocess.run(['c++','-std=c++20',str(folder/'fixture.cpp'),'-o',str(folder/'fixture')],check=True,capture_output=True,timeout=15)
        for args in [[],['fail']]:
            result=subprocess.run([str(folder/'fixture'),*args],check=True,capture_output=True,text=True,timeout=5)
            lines=result.stdout.splitlines();frames=[line[5:] for line in lines if line.startswith('@IHM ')]
            assert len(frames)==1,lines
            frame=json.loads(frames[0]);assert ('error' in frame)==bool(args)
            assert lines[0]=='[warn] evaluator diagnostic'
    print('PASS actual static dispatch framing with logging evaluator and physical exception; no native library')
if __name__=='__main__':main()
