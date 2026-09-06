"""Reproduce declared manual pixel picks; no fitted or invented material data.

Figure1b from Hong et al., Biology2022,11,1030, CC BY4.0.
Pixels are line-center estimates on the retained698x661 JPEG. Origin is a
zero-load engineering boundary anchor. See manifest for material/test limits.
"""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parent
IMAGE_SHA = '71c3630b2cd251667a4bf29832827ea318a8e7e5f6a7090a95c062e2ddef999e'
X = [395, 415, 435, 474, 513, 553, 592, 632, 651, 657]
Y = {
    'SM': [221, 213.5, 212, 211.5, 211.5, 211, 208, 203, 198.5, 198],
    'MM': [221, 206, 200, 199, 197, 194, 188.5, 176, 165, 157.5],
    'HM': [221, 173, 161.5, 158, 151, 140.5, 125.5, 86, 43, 30.5],
}


def main():
    image = ROOT / 'hong2022-figure1.jpg'
    assert hashlib.sha256(image.read_bytes()).hexdigest() == IMAGE_SHA
    rows = []
    for material, ys in Y.items():
        previous = -1
        for x, y in zip(X, ys):
            strain = .7 * (x - 395) / (671 - 395)
            pressure = 25000 * (221 - y) / (221 - 27)
            assert 0 <= strain < 1 and pressure >= previous
            previous = pressure
            rows.append(dict(material=material,compressive_strain=strain,stress_pa=pressure,
                             image_x_px=x,image_y_px=y,
                             basis='zero_load_boundary_anchor' if x == 395 else 'manual_figure_line_center'))
    target = ROOT / 'hong2022-compression.csv'
    with target.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print(json.dumps({'passed': True, 'rows':len(rows), 'maximum_strain':max(r['compressive_strain'] for r in rows),
                      'csv_sha256':hashlib.sha256(target.read_bytes()).hexdigest()}))


if __name__ == '__main__':main()
