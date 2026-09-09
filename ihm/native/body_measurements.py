"""Measure the mechanical body's SHAPE, not just its size.

``model_scaling.head_marker_height_m`` answers "how big"; every quantity here
answers "what proportions", which is what an anisotropic parameter changes and
what a single isotropic factor provably cannot.  Each is read by forward
kinematics from the model's own joint frames at its declared default pose --
the same walk that measures stature -- so a parameter that claims to widen the
pelvis has to show up here or it did not widen anything.

Naming, because one of these is easy to overclaim.  The frontal-plane angle of
the femoral mechanical axis is reported as ``femoral_obliquity_deg`` and NOT as
the Q-angle.  The clinical Q-angle is measured from ASIS through the patella
centre to the tibial tuberosity, and this model carries no ASIS marker and no
tibial tuberosity landmark.  Femoral obliquity is the anatomical quantity that
pelvic breadth actually moves, and it is a component of the Q-angle rather than
a synonym for it.
"""
from pathlib import Path
import math
import xml.etree.ElementTree as ET

import numpy as np

from .model_scaling import _vec3, default_pose_frames, head_marker_height_m


def joint_frames(model_path, pose=None):
    """World position of each joint's parent-side and child-side frame origins.

    Returns ``{joint: {'parent': xyz, 'child': xyz, 'parent_body': .., 'child_body': ..}}``.
    The two coincide at the default pose for a joint whose coordinates are all
    zero there; they are kept separate because the walker knee's are not.
    """
    root = ET.parse(model_path).getroot()
    pose = pose or default_pose_frames(model_path)
    out = {}
    for joint in root.find('.//JointSet/objects'):
        frames = {frame.get('name'): frame for frame in (joint.find('frames') or [])}
        record = {}
        for side in ('parent', 'child'):
            socket = joint.findtext('socket_%s_frame' % side) or ''
            frame = frames.get(socket.rsplit('/', 1)[-1])
            if frame is not None:
                body = (frame.findtext('socket_parent') or '').rsplit('/', 1)[-1]
                offset = _vec3(frame.findtext('translation'))
            else:
                body = socket.rsplit('/', 1)[-1]
                offset = np.zeros(3)
            body = 'ground' if body in ('', 'ground') else body
            origin, rotation = pose[body]
            record[side] = origin + rotation @ offset
            record[side + '_body'] = body
        out[joint.get('name')] = record
    return out


def segment_masses(model_path):
    root = ET.parse(model_path).getroot()
    return {body.get('name'): float(body.findtext('mass')) for body in root.iter('Body')}


def measure(model_path):
    """Every shape quantity this module knows how to read, in SI units."""
    model_path = Path(model_path)
    pose = default_pose_frames(model_path)
    centres = joint_frames(model_path, pose)
    masses = segment_masses(model_path)
    total_mass = sum(masses.values())

    def centre(joint, side='parent'):
        return centres[joint][side]

    hip_r, hip_l = centre('hip_r'), centre('hip_l')
    knee_r = centre('walker_knee_r')
    ankle_r = centre('ankle_r')
    acromial_r, acromial_l = centre('acromial_r'), centre('acromial_l')

    femur = knee_r - hip_r
    # Model axes: +x anterior, +y superior, +z to the subject's right.  The
    # frontal plane is y-z, so the femur's obliquity in it is the angle its
    # mechanical axis makes with vertical.
    obliquity = math.degrees(math.atan2(abs(femur[2]), abs(femur[1])))

    out = {
        'stature_m': head_marker_height_m(model_path),
        'hip_joint_separation_m': float(np.linalg.norm(hip_r - hip_l)),
        'acromial_separation_m': float(np.linalg.norm(acromial_r - acromial_l)),
        'femur_length_m': float(np.linalg.norm(femur)),
        'tibia_length_m': float(np.linalg.norm(ankle_r - knee_r)),
        'femoral_obliquity_deg': obliquity,
        'total_mass_kg': total_mass,
        # Stance geometry. A wider pelvis at the model's default pose does not
        # tilt the femur -- it translates the whole leg laterally -- so the
        # obliquity below does not move. What moves is where the foot lands.
        # The hip adduction needed to put the foot back where it was is the
        # honest Q-angle-relevant consequence, and it is derived in `compare`.
        'ankle_lateral_offset_m': float(ankle_r[2]),
        'hip_to_ankle_length_m': float(np.linalg.norm(ankle_r - hip_r)),
        'hip_lateral_offset_m': float(hip_r[2]),
    }
    out['shoulder_over_hip_breadth'] = (out['acromial_separation_m']
                                        / out['hip_joint_separation_m'])
    out['hip_separation_over_stature'] = out['hip_joint_separation_m'] / out['stature_m']
    out['femur_over_stature'] = out['femur_length_m'] / out['stature_m']
    out['tibia_over_stature'] = out['tibia_length_m'] / out['stature_m']
    out['leg_over_stature'] = ((out['femur_length_m'] + out['tibia_length_m'])
                              / out['stature_m'])
    out['segment_mass_fraction'] = {name: mass / total_mass
                                    for name, mass in sorted(masses.items())}
    return out


def compare(before, after):
    """Relative change of every scalar measurement, plus mass-fraction shift.

    ``hip_adduction_for_matched_stance_deg`` is the extra hip adduction that
    would put the foot back on its original line of progression after the pelvis
    moved.  It is the femoral obliquity change that pelvic breadth actually
    buys, and it is reported instead of a "Q-angle" the model cannot measure.
    """
    out = {}
    for key, value in before.items():
        if key == 'segment_mass_fraction':
            out[key] = {name: after[key][name] / value[name] - 1 for name in value}
        else:
            out[key] = after[key] / value - 1 if value else float('nan')
    shift = after['hip_lateral_offset_m'] - before['hip_lateral_offset_m']
    out['hip_adduction_for_matched_stance_deg'] = math.degrees(
        math.asin(max(-1.0, min(1.0, shift / after['hip_to_ankle_length_m']))))
    return out
