#!/usr/bin/env python3
"""Acquire pinned public kidney example bytes; never substitute it for final data."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import urllib.error
import urllib.request
import struct
import tempfile
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from ihm.anatomy.kidney_graph import read_amira, graph_statistics

COMMIT = "e35024b33018e010a057a3e69147459e7cadcc4e"
GRAPH_DIGESTS = {
    "Test.am": "349544f1070225529bfa4631231258baa896641f18afd4a62d5ddfcf3fbb0277",
    "Test2.am": "6073c2cce231a848f9da7add0151832002dc3b31d717f0f6d558c7ac7d0248ae",
}
BASE = f"https://raw.githubusercontent.com/HiPCTProject/Skeleton_analysis/{COMMIT}/"
SOURCES = {
    "Test.am": BASE+"Initial_Ordering/Test.am",
    "Test2.am": BASE+"Initial_Ordering/Test2.am",
    "upstream_reader.m": BASE+"Initial_Ordering/ultimate_amira_read.m",
    "upstream_ordering.m": BASE+"Initial_Ordering/run_ordering.m",
    "supplement.pdf": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs44303-025-00090-2/MediaObjects/44303_2025_90_MOESM1_ESM.pdf",
    "final_record_response.json": "https://zenodo.org/api/records/15275990",
    "earlier_record.json": "https://zenodo.org/api/records/7708966",
    "next_segmentation_record.json": "https://zenodo.org/api/records/19685382",
}

ARCHIVE_URL = "https://zenodo.org/records/19685382/files/final-data-repository-for-paper.zip"
ARCHIVE_SIZE = 51177854763


def fetch_range(start, end):
    request = urllib.request.Request(ARCHIVE_URL, headers={"Range":f"bytes={start}-{end}"})
    with urllib.request.urlopen(request,timeout=180) as response:
        if response.code != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end}/{ARCHIVE_SIZE}":
            raise ValueError("Server did not honor exact archive range")
        data = response.read()
    if len(data) != end-start+1:
        raise ValueError("Truncated archive range")
    return data


def collect_slab(root: Path, start: int, count: int):
    """Retain consecutive original 2D TIFFs and matching dense human labels."""
    raw = root/"data/raw/microstructure/kidney/segmentation_19685382"
    raw.mkdir(parents=True,exist_ok=True)
    tail_path=raw/"archive_central_directory.bin"
    if not tail_path.exists():
        tail_path.write_bytes(fetch_range(ARCHIVE_SIZE-6500000,ARCHIVE_SIZE-1))
    tail=tail_path.read_bytes()
    with tempfile.TemporaryFile() as sparse:
        sparse.seek(ARCHIVE_SIZE-len(tail)); sparse.write(tail)
        with zipfile.ZipFile(sparse) as archive:
            entries=archive.infolist()
    inventory=[dict(name=x.filename,size=x.file_size,compressed_size=x.compress_size,offset=x.header_offset,crc32=x.CRC,method=x.compress_type) for x in entries]
    (raw/"archive_inventory.json").write_text(json.dumps(inventory,indent=2)+"\n")
    selected=[]
    for category in ("images","labels"):
        prefix=f"final-data-repository-for-paper/competition-data/train/kidney_1_dense/{category}/"
        candidates=sorted([x for x in entries if x.filename.startswith(prefix) and x.file_size>0],key=lambda x:x.filename)
        if start < 0 or count < 1 or start+count > len(candidates):
            raise ValueError("Slab outside available slices")
        selected.extend(candidates[start:start+count])
    def acquire(entry):
        relative=Path(*Path(entry.filename).parts[-2:])
        target=raw/"kidney_1_dense"/relative
        if target.exists():
            content=target.read_bytes()
        else:
            # ZIP local headers are 30 bytes plus name and extra fields.
            header=fetch_range(entry.header_offset,entry.header_offset+1023)
            if header[:4]!=b"PK\x03\x04": raise ValueError("Invalid ZIP local header")
            name_length,extra_length=struct.unpack_from("<HH",header,26)
            offset=entry.header_offset+30+name_length+extra_length
            compressed=fetch_range(offset,offset+entry.compress_size-1)
            if entry.compress_type != zipfile.ZIP_DEFLATED: raise ValueError("Unsupported compression")
            content=zlib.decompress(compressed,-15)
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(content)
        if len(content)!=entry.file_size or zlib.crc32(content)!=entry.CRC:
            raise ValueError(f"ZIP size/CRC mismatch: {entry.filename}")
        return {"archive_member":entry.filename,"relative_path":str(target.relative_to(root)),"bytes":len(content),"sha256":hashlib.sha256(content).hexdigest(),"crc32":entry.CRC}
    with ThreadPoolExecutor(max_workers=4) as pool:
        files=list(pool.map(acquire,selected))
    receipt={"source_doi":"10.5281/zenodo.19685382","archive_url":ARCHIVE_URL,
             "archive_bytes":ARCHIVE_SIZE,"archive_md5_publisher":"7fee4bd2677f3f2447f910fb61a34862",
             "full_archive_md5_verified":False,"member_crc32_verified":True,"license":"CC-BY-4.0",
             "source_tier":"measured_other_donor","dataset":"kidney_1_dense","slice_start_index":start,
             "slice_count":count,"coverage":"contiguous slab, not whole kidney",
             "voxel_spacing_um":None,"spacing_status":"Requires paper/metadata verification; do not infer from filename",
             "files":files}
    (raw/"slab_manifest.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return receipt


def collect(root: Path, offline: bool = False) -> dict:
    raw = root / "data/raw/microstructure/kidney"
    derived = root / "data/derived/microstructure/kidney"
    raw.mkdir(parents=True, exist_ok=True)
    derived.mkdir(parents=True, exist_ok=True)
    receipt_path = raw / "acquisition.json"
    previous = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    receipt = {}
    for name, url in SOURCES.items():
        path = raw / name
        if offline:
            if name not in previous:
                raise ValueError(f"No recorded acquisition for {name}")
            status = previous[name]["http_status"]
        else:
            try:
                response = urllib.request.urlopen(url, timeout=90)
            except urllib.error.HTTPError as error:
                response = error
            with response:
                status = response.code
                path.write_bytes(response.read())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if name in GRAPH_DIGESTS and digest != GRAPH_DIGESTS[name]:
            raise ValueError(f"Pinned source digest mismatch: {name}")
        if offline and digest != previous[name]["sha256"]:
            raise ValueError(f"Checksum mismatch: {name}")
        receipt[name] = {"url":url,"http_status":status,"bytes":path.stat().st_size,"sha256":digest}
    receipt_path.write_text(json.dumps(receipt,indent=2)+"\n")
    graph = read_amira(raw / "Test.am")
    # Test2 has a corrupted UTF-8 unit label. Borrow units only after checking
    # that coordinates/connectivity/thickness equal the explicitly µm Test.am.
    ordered = read_amira(raw / "Test2.am", coordinate_unit="um")
    for key in ("nodes_mm","edges","points_mm","point_offsets","thickness_native"):
        if not np.array_equal(graph[key],ordered[key]):
            raise ValueError(f"Ordered example changed {key}")
    for key in ("strahler","topo"):
        graph[key] = ordered[key]
    np.savez_compressed(derived / "example_graph.npz",**graph)
    stats = graph_statistics(graph)
    stats.update({"schema_version":1,"source_id":"hipct-kidney-arterial",
                  "source_tier":"measured_other_donor","scope":"reduced_artery1_example",
                  "usable_for_population_priors":False,"complete_published_graph_acquired":False,
                  "source_commit":COMMIT,"input_sha256":{n:receipt[n]["sha256"] for n in ("Test.am","Test2.am")}})
    (derived / "statistics.json").write_text(json.dumps(stats,indent=2,allow_nan=False)+"\n")
    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline",action="store_true")
    parser.add_argument("--slab-slices",type=int,default=0)
    parser.add_argument("--slab-start",type=int,default=1100)
    args = parser.parse_args()
    stats = collect(ROOT,args.offline)
    print(json.dumps({k:stats[k] for k in ("scope","node_count","edge_count","point_count","complete_published_graph_acquired")},indent=2))
    if args.slab_slices:
        receipt=collect_slab(ROOT,args.slab_start,args.slab_slices)
        print(json.dumps({"slab_slices":receipt["slice_count"],"retained_bytes":sum(x["bytes"] for x in receipt["files"])}))
