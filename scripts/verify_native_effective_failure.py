"""Reproduce source normalization defect without changing native libraries."""
import math
from diagnose_native_effective_failure import fiber_primitive,audit

def main():
    lo=.1;lf=.14;fiso=900.;k=4.;e=.6;h=1e-7
    assert fiber_primitive(lf,lo,fiso,k,e)==0
    assert fiber_primitive(lf/lo,lo,fiso,k,e)>0
    derivative=(fiber_primitive((lf+h)/lo,lo,fiso,k,e)-fiber_primitive((lf-h)/lo,lo,fiso,k,e))/(2*h)
    passive=fiso*math.expm1(k*(lf/lo-1)/e)/math.expm1(k)
    assert math.isclose(derivative,passive,rel_tol=1e-8)
    r=audit();assert r['maximum_corrected_gradient_error']<r['maximum_original_gradient_error']
    assert r['maximum_corrected_gradient_error']>1e-4, 'Residual wrapped-path failure must remain visible'
    print('PASS: metres-versus-normalized energy defect reproduced; source correction preserves unresolved gate failure')
if __name__=='__main__':main()
