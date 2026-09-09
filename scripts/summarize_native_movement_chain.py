#!/usr/bin/env python3
"""Summarise a native movement chain: contact, load transfer and forward travel.

Reads only what the evaluator already retained. Reports the measurements that
decide whether a step happened, and never asserts that one did.
"""
import argparse, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def summarise(directory):
    chain = json.loads((Path(directory) / 'chain_report.json').read_text())
    rows, first_com, out = [], None, []
    for stage in chain['stages']:
        frames = json.loads((Path(stage['directory']) / 'trajectory.json').read_text())
        if not frames:
            continue
        if first_com is None:
            first_com = np.array(frames[0]['balance']['com_position_ground_m'])
        for f in frames:
            rows.append(f)
        segments = {}
        for f in frames:
            segments.setdefault(f['segment'], []).append(f)
        for name, block in segments.items():
            last = block[-1]
            com = np.array(last['balance']['com_position_ground_m'])
            out.append({
                'stage': stage['stage'], 'segment': name,
                'time_span_s': [block[0]['time_s'], last['time_s']],
                'end_right_clearance_m': last['per_foot_clearance_m']['r'],
                'end_left_clearance_m': last['per_foot_clearance_m']['l'],
                'min_right_clearance_m': min(b['per_foot_clearance_m']['r'] for b in block),
                'end_com_support_margin_m': last['balance']['com_support_margin_m'],
                'min_com_support_margin_m': min(b['balance']['com_support_margin_m'] for b in block),
                'end_active_foot_contacts': len(last['balance']['active_foot_contact_names']),
                'end_right_foot_contacts': sum(1 for n in last['balance']['active_foot_contact_names'] if n.endswith('_r')),
                'max_clipped_count': max(b['clipped_count'] for b in block),
                'max_state_reference_error_norm': max(b['state_reference_error_norm'] for b in block),
                'end_state_reference_error_norm': last['state_reference_error_norm'],
                'end_com_forward_travel_m': float(com[0] - first_com[0]),
                'end_com_speed_m_s': float(np.linalg.norm(last['balance']['com_velocity_ground_m_s'])),
                'end_pelvis_ty_m': last['coordinates']['pelvis_ty']['value'],
                'end_pelvis_tx_m': last['coordinates']['pelvis_tx']['value'],
                'end_knee_angle_l_rad': last['coordinates']['knee_angle_l']['value'],
                'end_foot_separation_forward_m': float(last['per_foot_centroid_ground_m']['r'][0] - last['per_foot_centroid_ground_m']['l'][0]),
            })
    final = rows[-1] if rows else None
    return {
        'chain': str(directory),
        'aborted': chain['aborted'],
        'all_stages_completed_horizon': chain['all_stages_completed_horizon'],
        'single_native_process': chain['single_native_process'],
        'chain_final_time_s': chain['chain_final_time_s'],
        'wall_seconds': chain['wall_seconds'],
        'segments': out,
        'final_com_forward_travel_m': out[-1]['end_com_forward_travel_m'] if out else None,
        'final_right_clearance_m': final['per_foot_clearance_m']['r'] if final else None,
        'walking_demonstrated': False,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', nargs='+')
    p.add_argument('--output')
    a = p.parse_args()
    result = [summarise(d) for d in a.directory]
    text = json.dumps(result if len(result) > 1 else result[0], indent=2)
    if a.output:
        Path(a.output).write_text(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
