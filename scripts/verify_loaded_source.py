"""Source receipt verification without executing native owners or module bodies."""
import importlib.util
from pathlib import Path
import sys,tempfile,types,unittest
from unittest.mock import patch
from ihm.assembly.interactive_scene import _loaded_source


class Tests(unittest.TestCase):
    def test_native_and_sensorimotor_generated_methods_are_receipted(self):
        from ihm.native import session
        from ihm.assembly import sensorimotor
        for module in (session,sensorimotor):
            receipt=_loaded_source(module)
            self.assertTrue(receipt['generated_code_sha256'])
            self.assertIn('dataclasses',receipt['code_generators']['modules'])
            self.assertEqual(receipt['code_generators']['python_version'],sys.version)

    def fixture(self,directory):
        path=Path(directory)/'receipt_fixture.py'
        path.write_text('from dataclasses import dataclass\n'
                        'IMPORT_COUNT=globals().get("IMPORT_COUNT",0)+1\n'
                        '@dataclass(frozen=True)\nclass Sample:\n'
                        '    value: int=3\n'
                        '    def double(self):return self.value*2\n')
        spec=importlib.util.spec_from_file_location('receipt_fixture',path)
        module=importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules,{'receipt_fixture':module}):spec.loader.exec_module(module)
        return module

    def test_no_module_reexecution_and_stale_user_method_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            module=self.fixture(directory);_loaded_source(module)
            self.assertEqual(module.IMPORT_COUNT,1)
            path=Path(module.__file__);path.write_text(path.read_text().replace('self.value*2','self.value*3'))
            with self.assertRaisesRegex(ValueError,'Loaded code differs'):_loaded_source(module)

    def test_generated_monkeypatch_and_defaults_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            module=self.fixture(directory)
            original=module.Sample.__init__
            replacement=types.FunctionType(original.__code__,original.__globals__,argdefs=(99,),closure=original.__closure__)
            replacement.__qualname__=original.__qualname__
            with patch.object(module.Sample,'__init__',replacement):
                with self.assertRaisesRegex(ValueError,'Loaded code differs'):_loaded_source(module)
            replacement=lambda self:'forged'
            replacement.__qualname__=module.Sample.__repr__.__qualname__
            with patch.object(module.Sample,'__repr__',replacement):
                with self.assertRaisesRegex(ValueError,'Loaded code differs'):_loaded_source(module)
            original=module.Sample.__repr__
            def cell(value):return (lambda:value).__closure__[0]
            closure=tuple(cell(lambda self:'forged') if name=='user_function' else value
                          for name,value in zip(original.__code__.co_freevars,original.__closure__))
            replacement=types.FunctionType(original.__code__,original.__globals__,closure=closure)
            replacement.__qualname__=original.__qualname__
            with patch.object(module.Sample,'__repr__',replacement):
                with self.assertRaisesRegex(ValueError,'Loaded code differs'):_loaded_source(module)


if __name__=='__main__':unittest.main()
