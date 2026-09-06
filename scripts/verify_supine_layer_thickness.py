"""Protect explicit shell thickness from exterior/raw-area migration."""
import math
from build_supine_surface_contact import layer_thickness

def main():
    raw=3.502598974493317;outer=1.7804602548390722
    explicit=[];legacy=[]
    for h in (.0001,.0015,.005):
        entity=dict(volume_m3=outer*h,shell={'thickness_m':h,'prior_source':'fixture engineering prior'},physical_surface_support={'area_m2':outer})
        value,basis=layer_thickness(entity,raw);explicit.append(value);assert basis['method']=='explicit_shell_thickness'
        del entity['shell']
        try:layer_thickness(entity,raw,'legacy_raw_full_skin_area_volume_ratio')
        except ValueError:pass
        else:raise AssertionError('New surface support silently used legacy area ratio')
        old={'volume_m3':raw*h}
        try:layer_thickness(old,raw)
        except ValueError:pass
        else:raise AssertionError('Unmarked legacy ratio accepted')
        legacy.append(layer_thickness(old,raw,'legacy_raw_full_skin_area_volume_ratio')[0])
    assert math.isclose(sum(explicit),.0066) and math.isclose(sum(legacy),.0066)
    assert math.isclose(.0066*outer/raw,.0033549480735623674)
    print('PASS: explicit6.6mm preserved; new missing thickness rejected; legacy requires named basis')
if __name__=='__main__':main()
