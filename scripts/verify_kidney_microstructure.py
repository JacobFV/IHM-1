#!/usr/bin/env python3
"""Check acquired source integrity, graph topology and measured statistics."""
import json
from pathlib import Path
import sys
import hashlib
import zlib
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/"scripts"))
import numpy as np
from collect_kidney_microstructure import collect, slab_metadata
from ihm.anatomy.kidney_graph import graph_statistics, validate_graph

def main():
    stats = collect(ROOT,offline=True)
    assert (stats["node_count"],stats["edge_count"],stats["point_count"]) == (148,147,19924)
    graph = dict(np.load(ROOT/"data/derived/microstructure/kidney/example_graph.npz"))
    # Connectivity, not just n-1 edges, establishes a single tree.
    visited, stack = set(), [0]
    adjacent = [[] for _ in graph["nodes_mm"]]
    for a,b in graph["edges"]:
        adjacent[a].append(b); adjacent[b].append(a)
    while stack:
        n = stack.pop()
        if n not in visited:
            visited.add(n); stack.extend(adjacent[n])
    assert len(visited)==148
    assert stats["global"]["path_over_chord_tortuosity"]["min"] >= 1-1e-5
    assert stats["radius_statistics"] is None and stats["murray_residual_statistics"] is None
    assert not stats["usable_for_population_priors"]
    for condition in stats["conditioned"].values():
        assert sum(x["segment_length_mm"]["count"] for x in condition.values())==147
    # Analytic curved polyline catches unit/path/chord confusion.
    fixture={"nodes_mm":np.array([[0.,0,0],[3,4,0]]),"edges":np.array([[0,1]]),
             "points_mm":np.array([[0.,0,0],[3,0,0],[3,4,0]]),"point_offsets":np.array([0,3]),
             "thickness_native":np.ones(3)}
    result=graph_statistics(fixture)
    assert result["global"]["segment_length_mm"]["mean"]==7
    assert result["global"]["path_over_chord_tortuosity"]["mean"]==1.4
    fixture["edges"][0,1]=2
    try: validate_graph(fixture)
    except ValueError: pass
    else: raise AssertionError("Out-of-range node accepted")
    slab_path=ROOT/"data/raw/microstructure/kidney/segmentation_19685382/slab_manifest.json"
    if slab_path.exists():
        from PIL import Image
        slab=json.loads(slab_path.read_text())
        pairs={}
        foreground=[]
        masks={}
        for member in slab["files"]:
            path=ROOT/member["relative_path"]
            content=path.read_bytes()
            assert len(content)==member["bytes"]
            assert hashlib.sha256(content).hexdigest()==member["sha256"]
            assert zlib.crc32(content)==member["crc32"]
            with Image.open(path) as im: array=np.array(im)
            assert array.shape==(1303,912)
            pairs.setdefault(path.stem,{})[path.parent.name]=array.shape
            if path.parent.name=="labels":
                assert set(np.unique(array)).issubset({0,255})
                foreground.append(int(np.count_nonzero(array)))
                masks[int(path.stem)]=array!=0
            else: assert array.dtype==np.uint16
        assert len(pairs)==slab["slice_count"]
        assert all(set(p)=={"images","labels"} for p in pairs.values())
        assert np.all(np.diff(sorted(int(n) for n in pairs))==1)
        assert sum(foreground)>0
        validation={"slice_count":len(pairs),"array_shape_zyx":[len(pairs),1303,912],
                    "foreground_voxels":sum(foreground),"mask_values":[0,255],
                    "image_dtype":"uint16","sha256_and_crc32_verified":True,
                    "voxel_spacing_um":None,"physical_morphometry_ready":False}
        evidence=ROOT/"data/raw/microstructure/kidney/segmentation_19685382/metadata/evidence_receipt.json"
        if evidence.exists():
            metadata=slab_metadata(ROOT)
            assert slab["voxel_spacing_um"]==metadata["voxel_spacing_um"]
            assert slab["donor"]==metadata["donor"]
            assert slab["metadata_evidence"]==metadata["metadata_evidence"]
            voxel_volume_mm3=float(np.prod(np.asarray(metadata["voxel_spacing_um"])/1000))
            volume=np.stack([masks[i] for i in sorted(masks)])
            faces={"z_min":volume[0],"z_max":volume[-1],"y_min":volume[:,0],
                   "y_max":volume[:,-1],"x_min":volume[:,:,0],"x_max":volume[:,:,-1]}
            border_counts={face:int(np.count_nonzero(mask)) for face,mask in faces.items()}
            validation.update(metadata)
            validation.update({"physical_morphometry_ready":False,
                "labeled_volume_measurement_ready":True,"voxel_volume_mm3":voxel_volume_mm3,
                "labeled_vessel_volume_mm3":float(sum(foreground)*voxel_volume_mm3),
                "slab_box_volume_mm3":float(volume.size*voxel_volume_mm3),
                "labeled_fraction_of_slab_box":float(sum(foreground)/volume.size),
                "boundary_face_foreground_voxels":border_counts,
                "boundary_censored":any(border_counts.values()),
                "volume_interpretation":"Binary labeled vessel volume inside the retained slab at the published voxel spacing; not in-vivo lumen volume or whole-organ vessel volume.",
                "density_denominator":"entire rectangular slab including background, not segmented kidney tissue",
                "unsupported_inferences":["complete vessel lengths","capillary completeness","in-vivo radii","organ-wide density","population priors"]})
        (ROOT/"data/derived/microstructure/kidney/slab_validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    print(json.dumps({"verified":True,"scope":stats["scope"],"edges":147,"full_graph_available":False,"paired_slab_verified":slab_path.exists()}))

if __name__=="__main__": main()
