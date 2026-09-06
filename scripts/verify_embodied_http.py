"""Real HTTP routing checks with a fake body owner; never launches native code."""
import json,threading,unittest
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from ihm.app import create_server

class Sessions:
    ident='a'*32
    def __init__(self):self.calls=[];self.closed=False
    def list(self):return {'sessions':[{'id':self.ident,'status':'initializing','closed':False}]}
    def create(self,data):self.calls.append(('create',data));return {'id':self.ident,'status':'initializing','closed':False}
    def command(self,ident,op,data=None):
        self.calls.append((ident,op,data))
        if op=='step':raise RuntimeError('uncertain owner command')
        return {'id':ident,'status':'closed' if op=='close' else 'initializing','closed':op=='close'}
    def close(self):self.closed=True

class Tests(unittest.TestCase):
    def test_pending_lifecycle_routing_and_error(self):
        server=create_server(Path(__file__).resolve().parents[1],port=0)
        owner=server.embodied=Sessions()
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        base=f'http://127.0.0.1:{server.server_port}/api/embodied/sessions'
        def request(url,data=None):
            req=Request(url,data=None if data is None else json.dumps(data).encode(),headers={'Content-Type':'application/json'})
            with urlopen(req,timeout=3) as response:return response.status,json.load(response)
        try:
            self.assertEqual(request(base,{})[0],201)
            self.assertEqual(request(base)[1]['sessions'][0]['status'],'initializing')
            self.assertEqual(request(base+'/'+owner.ident)[1]['status'],'initializing')
            with self.assertRaises(HTTPError) as error:request(base+'/'+owner.ident+'/step',{'sequence':0})
            self.assertEqual(error.exception.code,503)
            self.assertTrue(request(base+'/'+owner.ident+'/close',{})[1]['closed'])
        finally:server.shutdown();server.server_close();thread.join(timeout=2)
        self.assertTrue(owner.closed)

if __name__=='__main__':unittest.main()
