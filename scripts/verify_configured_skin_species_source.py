#!/usr/bin/env python3
"""Bounded preparation/seed identity checks; does not execute native transport."""
import json
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
import verify_configured_skin_species as verifier
from verify_configured_skin_species import prepare,STAGE,SNAPSHOT,ROOT,sha


class Checks(unittest.TestCase):
    def test_pinned_complete_native_species_seed(self):
        receipt=prepare()
        self.assertEqual(receipt['species_count'],28)
        values=json.loads(SNAPSHOT.read_text())['values']
        prefix='tissue.regional_skin.aggregate.species.'
        names=sorted(k[len(prefix):-8] for k in values if k.startswith(prefix) and k.endswith('.mass_ug'))
        self.assertEqual(receipt['species_names'],names)
        self.assertIn('Albumin',names)
        self.assertFalse(receipt['native_executed'])
        for name,digest in receipt['prepared_sha256'].items():self.assertEqual(sha(STAGE/name),digest)

    def test_changed_frozen_header_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            stage=Path(directory)/'frozen'
            with patch.object(verifier,'STAGE',stage):
                verifier.prepare()
                (stage/'native_configured_regional_skin.h').write_text('changed')
                with self.assertRaisesRegex(ValueError,'Frozen preparation differs'):
                    verifier.prepare()


if __name__=='__main__':unittest.main()
