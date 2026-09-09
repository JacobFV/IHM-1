"""Immutable owner-registered surface assets, including the real HTTP boundary."""
import gzip,hashlib,json,threading,unittest
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from ihm.app import create_server
from ihm.app.surface_assets import SurfaceAssetRegistry,register_body_surface_assets

class Assets(unittest.TestCase):
    def test_integrity_limits_atomicity_and_owner_capture(self):
        raw=b'{"weights":[1,0]}';digest=hashlib.sha256(raw).hexdigest();r=SurfaceAssetRegistry(maximum_asset_bytes=100,maximum_total_bytes=200)
        class Binding:
            def asset_records(self):return {digest:raw}
        class Mechanics:surface_binding=Binding()
        class Body:plant=Mechanics()
        register_body_surface_assets(Body(),r)
        self.assertIs(r.get(digest),raw);self.assertEqual(gzip.decompress(r.get(digest,compressed=True)),raw)
        r.register({digest:raw})
        with self.assertRaisesRegex(ValueError,'integrity'):r.register({'0'*64:raw})
        with self.assertRaisesRegex(ValueError,'envelope'):r.register({hashlib.sha256(b'x'*101).hexdigest():b'x'*101})
        with self.assertRaises(ValueError):r.get('../../etc/passwd')
        with self.assertRaises(KeyError):r.get('0'*64)
        small=SurfaceAssetRegistry(maximum_total_bytes=1)
        with self.assertRaisesRegex(ValueError,'capacity'):small.register({digest:raw})
        with self.assertRaises(KeyError):small.get(digest)
    def test_http_serves_only_retained_exact_bytes(self):
        server=create_server(Path(__file__).resolve().parents[1],port=0);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        raw=json.dumps({'weights':[.5,.5]*1000},separators=(',',':')).encode();digest=hashlib.sha256(raw).hexdigest();server.embodied.surface_assets.register({digest:raw})
        base=f'http://127.0.0.1:{server.server_port}/api/surface-binding/'
        try:
            for compressed in (False,True):
                with urlopen(Request(base+digest,headers={'Accept-Encoding':'gzip' if compressed else 'identity'})) as response:
                    data=response.read();self.assertEqual(gzip.decompress(data) if compressed else data,raw)
                    self.assertEqual(response.headers['ETag'],'"'+digest+'"');self.assertIn('immutable',response.headers['Cache-Control'])
            for target,status in [('0'*64,404),('not-a-digest',400),('%2e%2e/secret',400)]:
                with self.assertRaises(HTTPError) as result:urlopen(base+target)
                self.assertEqual(result.exception.code,status)
        finally:server.shutdown();server.server_close()

if __name__=='__main__':unittest.main()
