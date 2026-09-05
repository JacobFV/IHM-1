"""Plot saved native BioGears output, without fabricating or resampling waveforms."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source=Path('data/derived/physiology/biogears_native_run/native_multisystem.csv')
data=pd.read_csv(source)
channels=[('ArterialPressure(mmHg)','Arterial pressure (mmHg)'),('TotalLungVolume(mL)','Lung volume (mL)'),
          ('ArterialOxygenPressure(mmHg)','Arterial oxygen pressure (mmHg)'),('ArterialBloodPH','Arterial pH'),
          ('GlomerularFiltrationRate(mL/min)','Glomerular filtration (mL/min)'),('CoreTemperature(degC)','Core temperature (°C)')]
t=data['#Time(s)'].to_numpy(dtype=float)
assert len(t)>2500 and np.all(np.diff(t)>0)
fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True)
for axis,(column,label) in zip(axes.flat,channels,strict=True):
    y=data[column].to_numpy(dtype=float);assert np.isfinite(y).all()
    axis.plot(t,y,lw=1,color='#227e98');axis.set_ylabel(label)
    axis.grid(alpha=.18);axis.spines[['top','right']].set_visible(False)
    axis.ticklabel_format(axis='y',style='plain',useOffset=False)
for axis in axes[-1]:axis.set_xlabel('Seconds after engine stabilization')
fig.suptitle('Native BioGears • resting StandardMale\nPublished engine parameterization; execution check, not independent calibration',fontsize=14)
fig.tight_layout();Path('artifacts').mkdir(exist_ok=True)
fig.savefig('artifacts/native-multisystem-physiology.png',dpi=170)
metrics={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':len(t),'sample_interval_s':float(np.median(np.diff(t))),'independently_validated':False}
for column,label,prominence in [('ArterialPressure(mmHg)','heart',10),('TotalLungVolume(mL)','breathing',100)]:
    peaks,_=find_peaks(data[column].to_numpy(dtype=float),prominence=prominence)
    assert len(peaks)>2
    metrics[label+'_cycles_detected']=len(peaks)
    metrics[label+'_rate_per_min']=float(60/np.mean(np.diff(t[peaks])))
Path('artifacts/native-multisystem-metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print(json.dumps(metrics))
