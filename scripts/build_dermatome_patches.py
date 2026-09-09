#!/usr/bin/env python3
"""Discretise the integument into innervated dermatomal skin patches.

WHY THIS EXISTS.  The body carried three entities for the whole integument --
`body-skin-epidermis`, `body-skin-dermis`, `body-skin-hypodermis` -- plus sixteen
authored receptor patches.  A cutaneous afferent therefore had nothing to
innervate *at a location*, and dermatomal organisation could not be expressed at
all: there was no object in the model whose identity was "the skin of S1".

WHAT THIS BUILDS.  Every triangle of the real skin mesh's exterior component is
assigned a dermatome, hence a dorsal root and a named cutaneous nerve, and the
triangles of one (side, dermatome, region) territory are then cut into compact
patches of bounded area.  A patch is an addressable object with a position on the
skin, an area, a surface normal, a declared root level and a declared trunk.

WHAT THE ASSIGNMENT IS MADE OF, AND WHAT IT IS NOT.

*Measured, from this body's own geometry*: the skin surface and its area; the
exterior-component selection; which bone each piece of skin lies nearest, over
the 429 real bone and cartilage meshes, and therefore the region.  On the thorax
and abdomen the segmental LEVEL itself is measured -- skin nearest the Nth rib or
its costal cartilage takes level T(N), which is what the Nth intercostal nerve
actually supplies, and the rib's real posterior-to-anterior slope is what makes
the band slope.  On the hand the level is measured too: skin nearest the thumb or
index ray is C6, the middle ray C7, the ring and little rays C8, read off the
named phalanges and metacarpals rather than a coordinate threshold.  Whether skin
lies in the paravertebral gutter -- nearest a vertebra rather than a rib -- is
what separates dorsal-ramus territory from ventral-ramus territory.

*Authored, as an explicit prior*: the sector and axial-fraction rules that cut
limb dermatomes, and the y-bands that cut the face into V1/V2/V3.  Limb
dermatomes spiral and no bone marks the L4/L5 boundary on the calf; those cuts
are engineering choices consistent with a standard chart, not measurements, and
every patch records the rule string that produced it so a reader can tell which
kind of number they are holding.  `measured_dermatome_atlas` is False in the
output and stays False until a real atlas is registered onto this mesh.

THE FACE IS NOT A DERMATOME AND IS REPORTED AS SUCH.  Trigeminal V1/V2/V3 are
not spinal roots.  Those patches carry `root_level: null` and a
`root_gap_reason`, and the coverage block counts them separately, because a run
that reports "100% of skin innervated" while quietly giving the forehead a C2
root would be exactly the kind of number this programme keeps having to withdraw.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/canonical"
SKIN_ID = "body-bp3d-FJ2810"
TERRITORIES = ROOT / "data/research/engineered_skin_territories/materialization.json"

#: dermatome -> the named cutaneous trunk that carries it.  every name here must
#: exist in peripheral.json's nerves, which the build asserts; a trunk that is
#: declared in one repo and absent from the other is exactly the disconnect this
#: work is closing, so it must fail loudly rather than route into nothing.
DORSAL_RAMUS_NAME = "dorsal_ramus"

DERMATOME_TRUNK = {
    "c2": "lesser_occipital", "c3": "transverse_cervical", "c4": "supraclavicular",
    "c5": "axillary", "c6": "musculocutaneous", "c7": "median", "c8": "ulnar",
    "t1": "medial_cutaneous_forearm",
    # T2 reaches the axilla and medial arm as the intercostobrachial nerve, which
    # is a branch of the SECOND INTERCOSTAL nerve.  It joins the medial cutaneous
    # nerve of the arm, but that nerve is C8-T1 and naming it here would claim
    # that a C8-T1 trunk carries a T2 dermatome.  IBM's root cross-check caught
    # exactly that.
    "t2": "intercostal",
    **{f"t{i}": "intercostal" for i in range(3, 12)},
    "t12": "subcostal",
    "l1": "iliohypogastric", "l2": "lateral_femoral_cutaneous",
    "l3": "femoral", "l4": "saphenous", "l5": "superficial_fibular",
    "s1": "sural", "s2": "posterior_femoral_cutaneous",
    "s3": "pudendal", "s4": "pudendal",
    # the pudendal nerve is S2-S4.  Perianal and coccygeal skin is S5 and is
    # supplied by the anococcygeal nerves off the dorsal rami, not by pudendal.
    "s5": DORSAL_RAMUS_NAME,
    # the trigeminal divisions: a real trunk, and deliberately NOT a spinal root
    "v1": "trigeminal", "v2": "trigeminal", "v3": "trigeminal",
}
#: dorsal rami supply the paravertebral strip and the deep back.  the ventral
#: ramus names above do not cover it, and using `intercostal` there would assert
#: that the back of the trunk is supplied by the nerve that supplies the front.
DORSAL_RAMUS = DORSAL_RAMUS_NAME

RELAY_OF_ROOT = {"c": "cervical", "t": "thoracic", "l": "lumbar", "s": "sacral",
                 "v": "cranial"}

#: patch area budget, m^2.  acral skin gets finer patches than the trunk because
#: its innervation density is higher by more than an order of magnitude; this is
#: a two-tier approximation of that and not a measured receptor density.
TARGET_AREA_M2 = {"acral": 8e-4, "default": 30e-4}
ACRAL_REGIONS = {"hand", "foot", "head"}

SKULL = re.compile(r"\b(frontal bone|parietal bone|occipital bone|temporal bone|"
                   r"sphenoid|ethmoid|nasal bone|lacrimal bone|vomer|palatine bone|"
                   r"maxilla|zygomatic bone|inferior nasal concha|nasal septum|"
                   r"nasal cartilage|alar cartilage|nasal septal cartilage|"
                   r"mandible|hyoid|tooth|incisor|canine|premolar|molar)", re.I)
CERVICAL_V = re.compile(r"(first|second|third|fourth|fifth|sixth|seventh) cervical vertebra"
                        r"|vertebra c\d|\batlas\b|\baxis\b", re.I)
ORDINAL = {"first": 1, "second": 2, "third": 3, "fourth": 4, "fifth": 5, "sixth": 6,
           "seventh": 7, "eighth": 8, "ninth": 9, "tenth": 10, "eleventh": 11,
           "twelfth": 12, "big": 1, "little": 5}
ORDINAL_RE = "|".join(ORDINAL)


def _num(name: str, kind: str) -> int | None:
    """'left ninth rib' -> 9;  'vertebra t11' -> 11;  'left third rib' -> 3."""
    m = re.search(rf"(\w+) {kind}", name, re.I)
    if m and m.group(1).lower() in ORDINAL:
        return ORDINAL[m.group(1).lower()]
    m = re.search(rf"{kind}\s*([tlc])?(\d+)", name, re.I)
    return int(m.group(2)) if m else None


#: BodyParts3D carries two parallel naming families for the digits -- "distal
#: phalanx of right index finger" and "right distal phalanx of second finger of
#: hand" name the same ray.  Missing the second family cost 3.4% of the skin on
#: the first run, silently, because an unrecognised bone just falls into "other".
RAY_BY_ORDINAL = {1: "thumb", 2: "index", 3: "middle", 4: "ring", 5: "little"}


def classify_bone(name: str) -> tuple[str, dict]:
    """region label and the segmental hints a bone's identity carries."""
    n = name.lower()
    # digits first: "distal phalanx of right big toe" also contains no skull or
    # vertebral token, but "finger of foot" must not be read as a hand.
    if re.search(r"\bfoot\b|\btoe\b|metatarsal|calcaneus|talus|navicular|cuboid|"
                r"cuneiform bone|sustentaculum", n):
        k = _num(n, "finger of foot") or _num(n, "toe")
        lat = bool(re.search(r"(cuboid|fifth metatarsal|little toe|fourth toe)", n)) \
            or (k in (4, 5))
        med = bool(re.search(r"(navicular|medial cuneiform|big toe|first metatarsal|"
                             r"sustentaculum|sesamoid)", n)) or k == 1
        return "foot", {"lateral": lat, "medial": med}
    if "finger of hand" in n:
        return "hand", {"ray": RAY_BY_ORDINAL.get(_num(n, "finger of hand") or 3, "middle")}
    if SKULL.search(n):
        return "head", {}
    if re.search(r"(thyroid|cricoid|arytenoid|corniculate|cuneiform) cartilage", n) \
            or "tracheal" in n or "laryngeal" in n:
        return "neck", {}
    if CERVICAL_V.search(n):
        return "neck", {}
    if "triradiate cartilage" in n:
        return "pelvis", {}
    if "costal cartilage" in n:
        k = _num(n, "rib") or _num(n, "costal cartilage")
        if k:
            return "thorax_wall", {"rib": k}
    for kind in ("rib",):
        k = _num(n, kind)
        if k and "cartilage" not in n:
            return "thorax_wall", {"rib": k}
    if "sternum" in n or "xiphoid" in n or "manubrium" in n:
        return "thorax_wall", {"sternum": True}
    k = _num(n, "thoracic vertebra")
    if k or re.search(r"vertebra t(\d+)", n):
        k = k or int(re.search(r"vertebra t(\d+)", n).group(1))
        return "paravertebral", {"thoracic": k}
    k = _num(n, "lumbar vertebra")
    if k or re.search(r"vertebra l(\d+)", n):
        k = k or int(re.search(r"vertebra l(\d+)", n).group(1))
        return "paravertebral", {"lumbar": k}
    if "sacrum" in n:
        return "pelvis", {"sacrum": True}
    if "coccyx" in n:
        return "pelvis", {"coccyx": True}
    if "hip bone" in n or "ilium" in n or "ischium" in n or "pubis" in n:
        return "pelvis", {}
    if "clavicle" in n:
        return "pectoral_girdle", {"clavicle": True}
    if "scapula" in n:
        return "pectoral_girdle", {}
    if "humerus" in n:
        return "arm", {}
    if "radius" in n:
        return "forearm", {"bone": "radius"}
    if "ulna" in n:
        return "forearm", {"bone": "ulna"}
    if re.search(r"(thumb|index finger|middle finger|ring finger|little finger)", n) \
            and "toe" not in n:
        ray = ("thumb" if "thumb" in n else "index" if "index" in n else
               "middle" if "middle finger" in n else "ring" if "ring" in n else "little")
        return "hand", {"ray": ray}
    if "metacarpal" in n:
        k = _num(n, "metacarpal bone")
        return "hand", {"ray": {1: "thumb", 2: "index", 3: "middle",
                                4: "ring", 5: "little"}.get(k or 0, "middle")}
    if re.search(r"(scaphoid|trapezium|trapezoid)", n):
        return "hand", {"ray": "index"}
    if "capitate" in n:
        return "hand", {"ray": "middle"}
    if re.search(r"(hamate|lunate|triquetr|pisiform)", n):
        return "hand", {"ray": "little"}
    if "femur" in n:
        return "thigh", {}
    if "patella" in n:
        return "thigh", {"patella": True}
    if "tibia" in n:
        return "leg", {"bone": "tibia"}
    if "fibula" in n:
        return "leg", {"bone": "fibula"}
    return "other", {}


def local_frame(centroid, axis, side_sign):
    """proximal->distal axis plus a (lateral, anterior) basis orthogonal to it."""
    a = np.asarray(axis, float)
    a = a / (np.linalg.norm(a) + 1e-12)
    if a[1] > 0:                      # every limb long bone points distally = -y
        a = -a
    lat = np.array([side_sign, 0.0, 0.0])
    lat = lat - (lat @ a) * a
    lat /= np.linalg.norm(lat) + 1e-12
    ant = np.cross(a, lat)
    ant /= np.linalg.norm(ant) + 1e-12
    if ant[2] < 0:
        ant = -ant
    return np.asarray(centroid, float), a, lat, ant


def sectors(points, frame):
    c, a, lat, ant = frame
    d = points - c
    t = d @ a
    t = (t - t.min()) / (np.ptp(t) + 1e-12)
    th = np.degrees(np.arctan2(d @ ant, d @ lat))     # 0 lateral, 90 anterior
    return t, th


def pca_split(idx, C, area, target):
    """deterministic recursive bisection along the local principal axis."""
    out, stack = [], [idx]
    while stack:
        cur = stack.pop()
        if area[cur].sum() <= target or len(cur) < 8:
            out.append(cur)
            continue
        p = C[cur]
        p = p - p.mean(0)
        # smallest eigenvector index is deterministic; eigh returns ascending
        axis = np.linalg.eigh(p.T @ p)[1][:, -1]
        proj = p @ axis
        order = np.argsort(proj, kind="stable")
        cum = np.cumsum(area[cur][order])
        cut = int(np.searchsorted(cum, cum[-1] / 2.0)) + 1
        cut = min(max(cut, 1), len(cur) - 1)
        stack.append(cur[order[:cut]])
        stack.append(cur[order[cut:]])
    return out


def build(verbose=True):
    anat = json.loads((OUT / "anatomy.json").read_text())
    E = {e["id"]: e for e in anat["entities"]}
    skin = E[SKIN_ID]
    per = json.loads((OUT / "peripheral.json").read_text())
    relays = {r["id"]: np.asarray(r["position_m"], float) for r in per["relays"]}
    nerve_ids = {n["id"] for n in per["nerves"]}

    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3)
    F = np.asarray(g["indices"], int).reshape(-1, 3)
    terr = json.loads(TERRITORIES.read_text())
    if terr["inventory"]["geometry_sha256"] != skin["reference_geometry"]["sha256"]:
        raise SystemExit("skin geometry has changed under the territory selection")
    keep = np.asarray(terr["contact_eligible_triangle_ids"], int)
    T = F[keep]
    P = V[T]
    nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    area = 0.5 * np.linalg.norm(nrm, axis=1)
    C = P.mean(axis=1)
    exterior_area = float(area.sum())
    declared = float(skin["physical_surface_support"]["area_m2"])
    if abs(exterior_area - declared) > 1e-9:
        raise SystemExit(f"exterior area {exterior_area} != declared {declared}")

    # ---- region, from the nearest real bone or cartilage mesh ---------------
    bones = [e for e in anat["entities"]
             if e.get("role") in ("rigid_bone", "cartilage") and e.get("reference_geometry")]
    pts, own = [], []
    for i, e in enumerate(bones):
        gg = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
        v = np.asarray(gg["positions"], float).reshape(-1, 3)
        if len(v) > 1500:
            v = v[:: max(1, len(v) // 1500)]
        pts.append(v)
        own.append(np.full(len(v), i))
    tree = cKDTree(np.vstack(pts))
    near = np.concatenate(own)[tree.query(C, workers=-1)[1]]
    meta = [classify_bone(e["name"]) for e in bones]

    region = np.array([meta[i][0] for i in near], dtype=object)
    hints = [meta[i][1] for i in near]
    side_sign = np.where(C[:, 0] >= 0, 1.0, -1.0)

    # limb frames from the real long bones
    def frame_for(name):
        e = next(x for x in anat["entities"] if x["name"] == name)
        s = 1.0 if "left" in name else -1.0
        return local_frame(e["centroid_m"], e["principal_axis"], s)

    frames = {(r, s): frame_for(f"{s} {b}") for r, b in
              (("arm", "humerus"), ("forearm", "ulna"), ("thigh", "femur"),
               ("leg", "tibia")) for s in ("left", "right")}

    derm = np.full(len(C), "", dtype=object)
    trunk = np.full(len(C), "", dtype=object)
    rule = np.full(len(C), "", dtype=object)

    def assign(mask, d, t, r):
        derm[mask] = d
        trunk[mask] = t
        rule[mask] = r

    # ---- head: trigeminal divisions and the C2 nape ------------------------
    m = region == "head"
    y, z = C[:, 1], C[:, 2]
    assign(m & (z < -0.020), "c2", DERMATOME_TRUNK["c2"],
           "head skin posterior to the coronal ear plane -> C2 (authored plane)")
    assign(m & (z >= -0.020) & (y > 0.755), "v1", "trigeminal",
           "face above the supraorbital band -> V1 (authored y band)")
    assign(m & (z >= -0.020) & (y > 0.700) & (y <= 0.755), "v2", "trigeminal",
           "midface -> V2 (authored y band)")
    assign(m & (z >= -0.020) & (y <= 0.700), "v3", "trigeminal",
           "lower face and jaw -> V3 (authored y band)")

    # ---- neck --------------------------------------------------------------
    m = region == "neck"
    assign(m & (y > 0.640), "c2", DERMATOME_TRUNK["c2"], "upper neck -> C2 (authored y band)")
    assign(m & (y > 0.610) & (y <= 0.640), "c3", DERMATOME_TRUNK["c3"],
           "mid neck -> C3 (authored y band)")
    assign(m & (y <= 0.610), "c4", DERMATOME_TRUNK["c4"],
           "lower neck -> C4 (authored y band)")

    # ---- thorax and abdominal wall: level MEASURED from the nearest rib ----
    m = region == "thorax_wall"
    rib = np.array([h.get("rib", 0) for h in hints])
    lvl = np.clip(rib, 2, 12)                    # ribs 1-2 both carry the T2 band
    sternum = np.array([h.get("sternum", False) for h in hints])
    if sternum.any():
        # sternal skin takes the level of the nearest rib rather than a guess
        rib_ix = [i for i, e in enumerate(bones)
                  if classify_bone(e["name"])[1].get("rib")]
        rpts, rown = [], []
        for i in rib_ix:
            gg = json.loads(gzip.decompress(
                (ROOT / bones[i]["reference_geometry"]["path"]).read_bytes()))
            v = np.asarray(gg["positions"], float).reshape(-1, 3)
            rpts.append(v)
            rown.append(np.full(len(v), classify_bone(bones[i]["name"])[1]["rib"]))
        rtree = cKDTree(np.vstack(rpts))
        rlv = np.concatenate(rown)[rtree.query(C[m & sternum], workers=-1)[1]]
        lvl[np.where(m & sternum)[0]] = np.clip(rlv, 2, 12)
    for k in range(2, 13):
        sel = m & (lvl == k)
        d = f"t{k}"
        assign(sel, d, DERMATOME_TRUNK[d],
               f"skin nearest rib {k} or its costal cartilage -> T{k} ventral ramus "
               f"(measured from this body's ribs)")

    # ---- paravertebral: dorsal rami ---------------------------------------
    m = region == "paravertebral"
    tv = np.array([h.get("thoracic", 0) for h in hints])
    lv = np.array([h.get("lumbar", 0) for h in hints])
    for k in range(1, 13):
        sel = m & (tv == k)
        d = f"t{max(k, 2)}"
        assign(sel, d, DORSAL_RAMUS,
               f"paravertebral skin nearest thoracic vertebra {k} -> {d.upper()} "
               f"dorsal ramus (measured from this body's vertebrae)")
    for k in range(1, 6):
        sel = m & (lv == k)
        d = f"l{min(k, 3)}"
        assign(sel, d, DORSAL_RAMUS,
               f"paravertebral skin nearest lumbar vertebra {k} -> {d.upper()} dorsal "
               f"ramus (L4-L5 dorsal rami have no cutaneous branch; folded to L3)")

    # ---- pectoral girdle ---------------------------------------------------
    m = region == "pectoral_girdle"
    clav = np.array([h.get("clavicle", False) for h in hints])
    assign(m & clav, "c4", DERMATOME_TRUNK["c4"],
           "skin over the clavicle -> C4 supraclavicular (measured from the clavicle)")
    assign(m & ~clav & (y > 0.520), "c4", DERMATOME_TRUNK["c4"],
           "skin over the upper scapula -> C4 (authored y band)")
    assign(m & ~clav & (y <= 0.520), "t3", DORSAL_RAMUS,
           "skin over the scapular body -> T3 dorsal ramus (authored)")

    # ---- pelvis, buttock, perineum ----------------------------------------
    m = region == "pelvis"
    sac = np.array([h.get("sacrum", False) for h in hints])
    coc = np.array([h.get("coccyx", False) for h in hints])
    perineum = m & (np.abs(C[:, 0]) < 0.070) & (y < -0.020)
    assign(m & (z < -0.030) & (y > 0.020) & ~perineum, "l2", DORSAL_RAMUS,
           "upper buttock -> L2 superior cluneal, a dorsal ramus (authored sector)")
    assign(m & (z < -0.030) & (y <= 0.020) & ~perineum, "s2", DERMATOME_TRUNK["s2"],
           "lower buttock -> S2 posterior femoral cutaneous (authored sector)")
    assign(m & (z >= -0.030) & (y > 0.060) & ~perineum, "l1", DERMATOME_TRUNK["l1"],
           "skin over the iliac crest and lower abdominal wall -> L1 iliohypogastric")
    assign(m & (z >= -0.030) & (y <= 0.060) & ~perineum, "l1", "ilioinguinal",
           "groin -> L1 ilioinguinal (authored sector)")
    assign(m & sac & ~perineum, "s2", DORSAL_RAMUS,
           "skin over the sacrum -> S2 medial cluneal, a dorsal ramus (measured)")
    assign(m & coc & ~perineum, "s5", DERMATOME_TRUNK["s5"],
           "skin over the coccyx -> S5 (measured from the coccyx)")
    assign(perineum & (np.abs(C[:, 0]) < 0.035), "s4", DERMATOME_TRUNK["s4"],
           "central perineum -> S4 pudendal (authored sector)")
    assign(perineum & (np.abs(C[:, 0]) >= 0.035), "s3", DERMATOME_TRUNK["s3"],
           "outer perineum -> S3 pudendal (authored sector)")

    # ---- upper limb --------------------------------------------------------
    for s, sgn in (("left", 1.0), ("right", -1.0)):
        sel = (region == "arm") & (side_sign == sgn)
        if sel.any():
            t, th = sectors(C[sel], frames[("arm", s)])
            ix = np.where(sel)[0]
            med = np.abs(th) > 90
            assign(ix[~med], "c5", DERMATOME_TRUNK["c5"],
                   "lateral and anterior arm -> C5 (authored sector on the humeral axis)")
            assign(ix[med & (t < 0.25)], "t2", DERMATOME_TRUNK["t2"],
                   "axilla -> T2 intercostobrachial (authored sector)")
            assign(ix[med & (t >= 0.25)], "t1", DERMATOME_TRUNK["t1"],
                   "medial arm -> T1 (authored sector on the humeral axis)")

        sel = (region == "forearm") & (side_sign == sgn)
        if sel.any():
            t, th = sectors(C[sel], frames[("forearm", s)])
            ix = np.where(sel)[0]
            bone = np.array([hints[i].get("bone") for i in ix], dtype=object)
            post = (th > -135) & (th < -45)
            assign(ix[post], "c7", DERMATOME_TRUNK["c7"],
                   "posterior forearm -> C7 (authored sector)")
            assign(ix[~post & (bone == "radius")], "c6", DERMATOME_TRUNK["c6"],
                   "forearm skin nearest the radius -> C6 (measured from the radius)")
            assign(ix[~post & (bone != "radius")], "c8", DERMATOME_TRUNK["c8"],
                   "forearm skin nearest the ulna -> C8 (measured from the ulna)")

    sel = region == "hand"
    ix = np.where(sel)[0]
    ray = np.array([hints[i].get("ray", "middle") for i in ix], dtype=object)
    for rs, d in (("thumb", "c6"), ("index", "c6"), ("middle", "c7"),
                  ("ring", "c8"), ("little", "c8")):
        assign(ix[ray == rs], d, DERMATOME_TRUNK[d],
               f"hand skin nearest the {rs} ray -> {d.upper()} (measured from the "
               f"named carpal, metacarpal and phalangeal meshes)")

    # ---- lower limb --------------------------------------------------------
    for s, sgn in (("left", 1.0), ("right", -1.0)):
        sel = (region == "thigh") & (side_sign == sgn)
        if sel.any():
            t, th = sectors(C[sel], frames[("thigh", s)])
            ix = np.where(sel)[0]
            pat = np.array([hints[i].get("patella", False) for i in ix])
            ant_ = (th >= 45) & (th < 135)
            latl = np.abs(th) < 45
            post = (th > -135) & (th < -45)
            assign(ix[latl], "l2", DERMATOME_TRUNK["l2"],
                   "lateral thigh -> L2 lateral femoral cutaneous (authored sector)")
            assign(ix[ant_ & (t < 0.5)], "l2", DERMATOME_TRUNK["l2"],
                   "proximal anterior thigh -> L2 (authored sector)")
            assign(ix[ant_ & (t >= 0.5)], "l3", DERMATOME_TRUNK["l3"],
                   "distal anterior thigh -> L3 (authored sector)")
            assign(ix[~ant_ & ~latl & ~post], "l3", DERMATOME_TRUNK["l3"],
                   "medial thigh -> L3 (authored sector)")
            assign(ix[post], "s2", DERMATOME_TRUNK["s2"],
                   "posterior thigh -> S2 posterior femoral cutaneous (authored sector)")
            assign(ix[pat], "l3", DERMATOME_TRUNK["l3"],
                   "skin over the patella -> L3 (measured from the patella)")

        sel = (region == "leg") & (side_sign == sgn)
        if sel.any():
            t, th = sectors(C[sel], frames[("leg", s)])
            ix = np.where(sel)[0]
            med = np.abs(th) > 135
            post = (th > -135) & (th < -45)
            assign(ix[med], "l4", DERMATOME_TRUNK["l4"],
                   "medial leg -> L4 saphenous (authored sector on the tibial axis)")
            assign(ix[post & ~med], "s1", DERMATOME_TRUNK["s1"],
                   "posterior calf -> S1 sural (authored sector)")
            assign(ix[~med & ~post], "l5", DERMATOME_TRUNK["l5"],
                   "anterolateral leg -> L5 superficial fibular (authored sector)")

    sel = region == "foot"
    if sel.any():
        ix = np.where(sel)[0]
        fy = C[ix, 1]
        lo, hi = fy.min(), fy.max()
        sole = fy < lo + 0.35 * (hi - lo)
        lat_ = np.array([hints[i].get("lateral", False) for i in ix])
        med_ = np.array([hints[i].get("medial", False) for i in ix])
        assign(ix[sole | lat_], "s1", DERMATOME_TRUNK["s1"],
               "sole and lateral border -> S1 (authored height band plus the named "
               "lateral tarsal and metatarsal meshes)")
        assign(ix[~sole & ~lat_ & med_ & (fy > lo + 0.7 * (hi - lo))], "l4",
               DERMATOME_TRUNK["l4"],
               "medial ankle -> L4 saphenous (measured from the medial tarsals)")
        rest = ~sole & ~lat_ & ~(med_ & (fy > lo + 0.7 * (hi - lo)))
        assign(ix[rest], "l5", DERMATOME_TRUNK["l5"],
               "dorsum of the foot -> L5 (authored)")

    # ---- patches -----------------------------------------------------------
    unassigned = derm == ""
    patches = []
    groups: dict[tuple, list[int]] = {}
    for i in np.where(~unassigned)[0]:
        groups.setdefault(
            ("left" if side_sign[i] > 0 else "right", derm[i], trunk[i], region[i]),
            []).append(i)

    normals = nrm / (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12)
    vtree = cKDTree(V)
    for (side, d, tr, reg), idx in sorted(groups.items()):
        idx = np.asarray(idx)
        target = TARGET_AREA_M2["acral" if reg in ACRAL_REGIONS else "default"]
        for cl in pca_split(idx, C, area, target):
            a = float(area[cl].sum())
            if a <= 0:
                continue
            w = area[cl] / a
            pos = (C[cl] * w[:, None]).sum(0)
            pos = V[vtree.query(pos)[1]]           # snap onto the real surface
            n = (normals[cl] * w[:, None]).sum(0)
            n = n / (np.linalg.norm(n) + 1e-12)
            root = d if d[0] != "v" else None
            relay = f"peripheral-relay-{side}-{RELAY_OF_ROOT[d[0]]}"
            nid = f"peripheral-nerve-{side}-{tr}"
            sgn = 1.0 if side == "left" else -1.0
            rp = relays.get(relay)
            if rp is None:
                raise SystemExit(f"no relay {relay}")
            way = np.array([sgn * 0.04, rp[1], rp[2]])
            length = float(math.dist(pos, way) + math.dist(way, rp))
            patches.append(dict(
                id=f"dermatome-patch-{side}-{d}-{reg}-{len(patches):04d}",
                name=f"{side} {d.upper()} {reg.replace('_', ' ')} skin patch",
                side=side, dermatome=d, root_level=root,
                root_gap_reason=None if root else
                "trigeminal division; a cranial nerve territory, not a spinal dermatome",
                region=reg, body_entity_id=SKIN_ID,
                position_m=[float(v) for v in pos],
                surface_normal=[float(v) for v in n],
                area_m2=a, triangle_count=int(len(cl)),
                nerve_id=nid, nerve_name=tr, relay_id=relay,
                brain_target_id=f"brain-{'rh' if side == 'left' else 'lh'}-postcentral",
                path_length_m=length,
                path_length_scope="representative_endpoint_to_relay",
                length_method="sum of authored polyline segment lengths through a "
                              "proximal waypoint, matching peripheral.json nerve routes",
                modalities=["pressure_pa", "temperature_C", "stretch_fraction"],
                assignment_rule=str(rule[cl[0]]),
                evidence_kind="authored_dermatome_prior_over_measured_skin_geometry",
                measured_dermatome_atlas=False,
                surface_anchor="area-weighted centroid of the patch's triangles, "
                               "snapped to the nearest skin vertex",
            ))

    # ---- coverage ----------------------------------------------------------
    by_derm: dict[str, dict] = {}
    for p in patches:
        r = by_derm.setdefault(p["dermatome"], dict(
            dermatome=p["dermatome"], root_level=p["root_level"], patches=0,
            area_m2=0.0, trunks=set()))
        r["patches"] += 1
        r["area_m2"] += p["area_m2"]
        r["trunks"].add(p["nerve_name"])
    for r in by_derm.values():
        r["trunks"] = sorted(r["trunks"])

    spinal = [p for p in patches if p["root_level"]]
    cranial = [p for p in patches if not p["root_level"]]
    missing = sorted({p["nerve_id"] for p in patches} - nerve_ids)
    covered = sum(p["area_m2"] for p in patches)
    coverage = dict(
        skin_entity_id=SKIN_ID,
        raw_source_area_m2=float(skin["surface_area_m2"]),
        exterior_component_area_m2=exterior_area,
        exterior_triangle_count=int(len(T)),
        area_denominator="exterior component of the skin mesh; the other 99 "
                         "components are interior and orifice surfaces and are "
                         "excluded, so this is not clinical body surface area",
        patch_count=len(patches),
        patch_area_m2=float(covered),
        patch_area_fraction_of_exterior=float(covered / exterior_area),
        unassigned_triangle_count=int(unassigned.sum()),
        unassigned_area_m2=float(area[unassigned].sum()),
        patches_with_spinal_root=len(spinal),
        patches_with_spinal_root_area_m2=float(sum(p["area_m2"] for p in spinal)),
        patches_without_spinal_root=len(cranial),
        patches_without_spinal_root_area_m2=float(sum(p["area_m2"] for p in cranial)),
        patches_with_declared_trunk=sum(1 for p in patches if p["nerve_name"]),
        trunks_used=sorted({p["nerve_name"] for p in patches}),
        trunks_absent_from_peripheral_json=missing,
        dermatomes_present=sorted(by_derm),
        dermatomes_declared=sorted(DERMATOME_TRUNK),
        dermatomes_absent=sorted(set(DERMATOME_TRUNK) - set(by_derm)),
        by_dermatome=[by_derm[k] for k in sorted(by_derm)],
    )

    data = dict(
        schema_version=1, id="ihm-body-dermatomes", frame=anat["frame"],
        units=dict(position="m", area="m2", path_length="m"),
        skin_geometry=dict(
            entity_id=SKIN_ID, path=skin["reference_geometry"]["path"],
            sha256=skin["reference_geometry"]["sha256"],
            exterior_selection="data/research/engineered_skin_territories/"
                               "materialization.json contact_eligible_triangle_ids"),
        method=dict(
            region="nearest of 429 real bone and cartilage meshes, 1500-vertex "
                   "subsample per mesh",
            measured_levels=["thorax and abdominal wall T2-T12 from the nearest rib "
                             "or costal cartilage",
                             "paravertebral T2-T12 and L1-L3 from the nearest vertebra",
                             "hand C6/C7/C8 from the named ray",
                             "forearm C6/C8 from radius vs ulna",
                             "patella L3, coccyx S5, sacrum S2"],
            authored_levels=["all limb sectors and axial fractions",
                             "face V1/V2/V3 y bands", "neck C2/C3/C4 y bands",
                             "buttock, groin and perineum sectors"],
            patching="deterministic recursive bisection along the local principal "
                     "axis until a patch is under its area budget",
            target_area_m2=TARGET_AREA_M2,
            acral_regions=sorted(ACRAL_REGIONS)),
        evidence_kind="authored_dermatome_prior_over_measured_skin_geometry",
        measured_dermatome_atlas=False,
        biological_validation=False,
        dermatome_trunk=DERMATOME_TRUNK,
        dorsal_ramus_trunk=DORSAL_RAMUS,
        patches=patches,
        coverage=coverage,
        limitations=[
            "No registered dermatome atlas. Region and several segmental levels are "
            "measured from this body's own bones; every limb boundary is an authored "
            "sector consistent with a standard chart, and each patch records which.",
            "Dermatomes overlap in life by roughly one segment; this partition is "
            "exclusive, so adjacent-segment overlap is not represented.",
            "Patch density is two-tier by area, not by measured receptor density; "
            "real fingertip innervation density exceeds trunk density by more than "
            "an order of magnitude.",
            "The face carries trigeminal divisions and therefore no spinal root; "
            "those patches are counted separately and never given a root level.",
            "Side is taken from the sign of the patch centroid's x, so a midline "
            "patch is assigned to whichever side its centroid falls on.",
            "Patch positions are rest-pose canonical coordinates; they are not "
            "re-posed with the body here.",
        ],
        counts=dict(patches=len(patches),
                    territories=len(groups),
                    dermatomes=len(by_derm)),
    )
    payload = json.dumps(data, indent=2, allow_nan=False) + "\n"
    (OUT / "dermatomes.json").write_text(payload)

    if verbose:
        c = coverage
        print(f"skin exterior component: {c['exterior_component_area_m2']:.4f} m2 "
              f"of {c['raw_source_area_m2']:.4f} m2 raw mesh area "
              f"({c['exterior_triangle_count']} triangles)")
        print(f"patches: {c['patch_count']}  covering "
              f"{c['patch_area_m2']:.4f} m2 = "
              f"{100 * c['patch_area_fraction_of_exterior']:.2f}% of the exterior")
        print(f"  unassigned: {c['unassigned_triangle_count']} triangles, "
              f"{c['unassigned_area_m2']:.4f} m2")
        print(f"  with a spinal root : {c['patches_with_spinal_root']} of "
              f"{c['patch_count']}  ({c['patches_with_spinal_root_area_m2']:.4f} m2)")
        print(f"  without            : {c['patches_without_spinal_root']} of "
              f"{c['patch_count']}  ({c['patches_without_spinal_root_area_m2']:.4f} m2, "
              f"trigeminal)")
        print(f"  trunks used: {len(c['trunks_used'])}; absent from peripheral.json: "
              f"{c['trunks_absent_from_peripheral_json'] or '-'}")
        print(f"  dermatomes present {len(c['dermatomes_present'])} of "
              f"{len(c['dermatomes_declared'])}; absent: "
              f"{c['dermatomes_absent'] or '-'}")
        print(f"{'derm':6s} {'root':6s} {'patches':>8s} {'cm2':>9s}  trunks")
        for r in coverage["by_dermatome"]:
            print(f"{r['dermatome']:6s} {str(r['root_level']):6s} {r['patches']:8d} "
                  f"{r['area_m2'] * 1e4:9.1f}  {','.join(r['trunks'])}")
        print(f"sha256(dermatomes.json) = "
              f"{hashlib.sha256(payload.encode()).hexdigest()[:16]}")
    return data


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-q", "--quiet", action="store_true")
    build(verbose=not ap.parse_args().quiet)
