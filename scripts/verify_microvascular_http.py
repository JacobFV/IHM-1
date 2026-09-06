"""Actual loopback route with bounded fake patch owner; no geometry/native run."""
import json,threading,unittest
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from ihm.app import create_server

class Tests(unittest.TestCase):
    def test_patch_route_keeps_generation_separate_from_body_commands(self):
        server=create_server(Path(__file__).resolve().parents[1],port=0)
        class Patches:
            def materialize(self,data):
                if data.get('busy'):
                    from ihm.app.microvascular_patch import PatchRequestError
                    raise PatchRequestError(503,'busy','One query at a time')
                return {'schema':'fixture','request':data}
        server.microvascular_patches=Patches()
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        address=f'http://127.0.0.1:{server.server_port}/api/body/microvascular-patch'
        try:
            request=Request(address,data=json.dumps({'entity_id':'body-bp3d-FJ1442'}).encode(),headers={'Content-Type':'application/json'})
            with urlopen(request,timeout=3) as response:
                self.assertEqual(response.status,200)
                self.assertEqual(json.load(response)['request']['entity_id'],'body-bp3d-FJ1442')
            busy=Request(address,data=b'{"busy":true}',headers={'Content-Type':'application/json'})
            with self.assertRaises(HTTPError) as error:urlopen(busy,timeout=3)
            self.assertEqual(error.exception.code,503)
            self.assertEqual(json.load(error.exception)['code'],'busy')
            bad=Request(address,data=b'{}',headers={'Content-Type':'text/plain'})
            with self.assertRaises(HTTPError) as error:urlopen(bad,timeout=3)
            self.assertEqual(error.exception.code,415)
        finally:server.shutdown();server.server_close();thread.join(timeout=2)

if __name__=='__main__':unittest.main()
