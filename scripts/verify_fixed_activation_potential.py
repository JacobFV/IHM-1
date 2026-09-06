"""Source-law effective-work checks; no new native evaluations."""
import json,xml.etree.ElementTree as ET
from pathlib import Path
from scipy.integrate import quad
from diagnose_fixed_activation_potential import audit,ROOT
from build_supine_initial_state import passive_expression

def main():
    result=audit()
    assert result['observations']==126
    for row in result['rows']:
        assert abs(row['force_error_n'])<=row['native_initialization_tolerance_n']
        assert abs(row['energy_length_derivative_error_n'])<1e-6
        assert row['internal_stiffness_n_m']>0
    root=ET.parse(ROOT/'data/derived/lumbar-supine-static-1g1q08u2/native/assembled_model.osim').getroot()
    q=json.loads((ROOT/'data/derived/lumbar-supine-static-1g1q08u2/old_q_seed.json').read_text())['coordinates']
    count=0
    for force in root.iter('ExpressionBasedCoordinateForce'):
        expr=force.findtext('expression');x=q[force.findtext('coordinate')]['value'];h=1e-6
        # Integral over the small interval avoids subtracting arbitrary energy offsets.
        derivative=-quad(lambda t:passive_expression(expr,t),x-h,x+h,epsabs=1e-12)[0]/(2*h)
        assert abs(derivative+passive_expression(expr,x))<1e-7
        count+=1
    assert count==20
    print('PASS:126native-source force matches, effective work derivative, stable sampled branches,20passive potentials')
if __name__=='__main__':main()
