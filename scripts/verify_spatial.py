"""Physical-frame and conservative-map behavioral checks."""
import numpy as np
from ihm.spatial import Frame, Registration, fit_landmarks, conservative_map

rng=np.random.default_rng(51);x=rng.normal(size=(12,3))
r=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
y=x@r.T+np.array([1,2,-1])
a=Frame('a','m','specimen1');b=Frame('b','mm','specimen1')
reg=fit_landmarks(x,y*1000,a,b,reference='analytic known transform')
assert np.allclose(reg.apply(x),y*1000,atol=1e-9)
assert reg.rms_residual_m<1e-12
assert np.allclose(reg.inverse().apply(y*1000),x)
for bad in [np.zeros((4,3)),np.arange(4)[:,None]*np.ones((1,3))]:
 try:fit_landmarks(bad,bad,a,a,reference='degenerate')
 except ValueError:pass
 else:raise AssertionError('rank deficient landmarks accepted')
try:Frame('unknown','source-units','case')
except ValueError:pass
else:raise AssertionError('unknown physical units accepted')
try:Registration(a,b,np.diag([2.,2,2,1]),0,'unscaled affine')
except ValueError:pass
else:raise AssertionError('unexplained scaling accepted')
w=conservative_map([[.8,.2,0],[0,.3,.7]])
assert np.allclose(w.intensive([7,7,7]),[7,7])
assert np.isclose(w.extensive([3,4]).sum(),7)
assert np.allclose(w.extensive([3,4]),[2.4,1.8,2.8])
print('verified: physical frame units, known registration/inverse, degenerate rejection, constant and extensive mapping')
from pathlib import Path
from ihm.spatial.vtk import read_arrays
import tempfile
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'points.vtp';p.write_text('<VTKFile><PolyData><Piece><Points><DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">0 0 0 1 2 3</DataArray></Points></Piece></PolyData></VTKFile>')
    d=read_arrays(p);assert np.array_equal(d['Points/Points'],[[0,0,0],[1,2,3]])
print('verified: VTK original array reader')
for construct in [lambda:conservative_map([[1.000009]]),lambda:Registration(a,b,np.diag([1,1,1,1.000009]),0,'invalid homogeneous row')]:
    try:construct()
    except ValueError:pass
    else:raise AssertionError('relative tolerance admitted a nonconservative/inhomogeneous mapping')
tri=np.array([[0,0,0],[1,0,0],[0,1,0.]])
weighted=fit_landmarks(tri,tri+100,a,a,reference='large finite weights',weights=[1e308]*3)
assert np.allclose(weighted.apply(tri),tri+100)
from ihm.spatial.vtk import surface
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'named.vtp'
    p.write_text('<VTKFile><PolyData><Piece><Points><DataArray type="Float64" Name="Array 0A144E50" NumberOfComponents="3">0 0 0 1 0 0 0 1 0</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity">0 1 2</DataArray><DataArray type="Int32" Name="offsets">3</DataArray></Polys></Piece></PolyData></VTKFile>')
    v,f=surface(p);assert v.shape==(3,3) and f.tolist()==[[0,1,2]]
    p.write_text(p.read_text().replace('>3</DataArray>','>2</DataArray>'))
    try:surface(p)
    except ValueError:pass
    else:raise AssertionError('malformed polygon offset accepted')
