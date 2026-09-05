from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.skin_transport import build_skin_transport

if __name__=='__main__':
    data=build_skin_transport(Path(__file__).resolve().parents[1])
    print('Built isolated forearm skin transport:',data['model']['attachment']['unit_id'])
