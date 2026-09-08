"""Short real native/IBM/physiology run with an attached bedroom owner."""
from pathlib import Path
import json,tempfile,time
from ihm.assembly.embodied import EmbodiedRuntime
root=Path(__file__).resolve().parents[1]
output=Path(tempfile.mkdtemp(prefix='environment-embodied-',dir=root/'data/derived'))
body=None;start=time.monotonic()
try:
    body=EmbodiedRuntime.from_workspace(root,output/'body',environment='supine',environment_selection={'scene':'bedroom','objects':[]})
    initial=body.snapshot()
    assert initial['environment_state']['time_s']==0
    for _ in range(3):frame=body.step({})
    assert abs(frame['time_s']-.06)<1e-10
    assert abs(frame['environment_state']['time_s']-frame['mechanics']['time_s'])<1e-10
    assert frame['environment_state']['contact_count']>0
    report={'passed':True,'duration_s':frame['time_s'],'contact_count':frame['environment_state']['contact_count'],
            'wall_s':time.monotonic()-start,'scope':'Three actual articulated/IBM/native physiology exchanges; no long-horizon validation'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'output':str(output),**report},indent=2))
except BaseException as error:
    (output/'report.json').write_text(json.dumps({'passed':False,'error':str(error),'wall_s':time.monotonic()-start},indent=2)+'\n')
    raise
finally:
    if body is not None:body.close()
