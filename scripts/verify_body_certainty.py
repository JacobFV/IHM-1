from ihm.assembly.certainty import Estimate,fuse,entity_certainty
import numpy as np

def e(id,group,mean=3,var=4,unit='m'):return Estimate(id,group,mean,var,unit,'analytic test')
a=fuse([e('original','shared'),e('derived','shared')]);assert a['variance']==4 and a['mean']==3
b=fuse([e('a','independent1'),e('b','independent2')]);assert b['variance']==2 and b['mean']==3
c=fuse([e('a','shared',2),e('b','shared',4)]);assert c['mean']==3 and c['variance']==4
for estimates in [[e('a','a'),e('a','b')],[e('a','a'),e('b','b',unit='mm')]]:
 try:fuse(estimates)
 except ValueError:pass
 else:raise AssertionError('unsupported evidence accepted')
for v in [None,0,-1,True,float('nan')]:
 try:e('x','x',var=v)
 except ValueError:pass
 else:raise AssertionError('unjustified variance accepted')
s=entity_certainty({'bounds_m':{'min':[-1,-1,-1],'max':[1,1,1]},'uncertainty':{'biological':{'status':'unquantified'}}})
assert s['biological']['status']=='unquantified'
rng=np.random.default_rng(11);x=rng.uniform(-1,1,10000)
assert np.max(np.abs(x.astype(np.float32).astype(float)-x))<=s['display_numerics']['gpu_coordinate_rounding_bound_m']
print('Verified dependence-aware fusion, exact units, unknown-variance rejection and separate display precision')
