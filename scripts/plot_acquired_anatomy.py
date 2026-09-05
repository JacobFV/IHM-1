"""render imported source anatomy and archived CFD; no invented geometry or flow."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def obj_surface(path):
    vertices=[];faces=[]
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('v '): vertices.append([float(x) for x in line.split()[1:4]])
            elif line.startswith('f '):
                f=[int(x.split('/')[0])-1 for x in line.split()[1:]]
                for j in range(1,len(f)-1): faces.append([f[0],f[j],f[j+1]])
    vertices=np.asarray(vertices);faces=np.asarray(faces)
    # Geometry decimation is visualization-only; original mesh is retained.
    mesh=trimesh.Trimesh(vertices=vertices,faces=faces,process=False)
    if len(faces)>1000: mesh=mesh.simplify_quadric_decimation(face_count=1000)
    return vertices,mesh.triangles


def main():
    atlas=json.loads(Path('data/derived/anatomy/bodyparts3d_index.json').read_text())
    out=Path('artifacts');out.mkdir(exist_ok=True)
    fig=plt.figure(figsize=(15,9),facecolor='#fafbfc')
    for number,concept,color in ((1,'bone organ','#91806a'),(2,'muscle organ','#aa4f5f')):
        ax=fig.add_subplot(1,3,number,projection='3d');bounds=[];n=0
        for mesh in atlas['meshes']:
            if not any(c['name']==concept for c in mesh['concepts']):continue
            vertices,faces=obj_surface(mesh['source_path']);bounds.append(vertices);n+=1
            ax.add_collection3d(Poly3DCollection(faces,facecolor=color,edgecolor='none',alpha=1.))
        points=np.concatenate(bounds);low=points.min(axis=0);high=points.max(axis=0)
        ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1]);ax.set_zlim(low[2],high[2]);ax.set_box_aspect(high-low)
        ax.view_init(elev=4,azim=-90);ax.set_axis_off()
        ax.set_title(f'BodyParts3D: {n} {concept} meshes\nNative atlas coordinates',fontsize=12,pad=15)
    ax=fig.add_subplot(1,3,3,projection='3d')
    file=Path('data/derived/vascular/0050_H_CERE_H_first_cfd_frame.npz')
    with np.load(file,allow_pickle=False) as d:
        points=d['Points/Points'];speed=np.linalg.norm(d['PointData/velocity_01010'],axis=1)
    sample=np.arange(0,len(points),3)
    sc=ax.scatter(*points[sample].T,c=speed[sample],s=1,cmap='turbo',alpha=.65,rasterized=True)
    extent=np.ptp(points,axis=0);ax.set_box_aspect(extent);ax.view_init(elev=15,azim=-70);ax.set_axis_off()
    ax.set_title('VMR cerebral case 0050\nArchived 3D simulated velocity, state 01010',fontsize=12,pad=15)
    cbar=fig.colorbar(sc,ax=ax,shrink=.4,pad=.02);cbar.set_label('Speed magnitude (source units; confirmation pending)',fontsize=9)
    fig.suptitle('Actual acquired anatomy and flow assets',fontsize=19,y=.96)
    fig.text(.5,.025,'Separate source frames and subjects • not registered into one human • mesh display decimated; original assets retained',ha='center',fontsize=11)
    fig.subplots_adjust(left=.02,right=.96,top=.86,bottom=.08,wspace=.05)
    fig.savefig(out/'acquired-anatomy-and-flow.png',dpi=170)
    print(out/'acquired-anatomy-and-flow.png')


if __name__=='__main__':main()
