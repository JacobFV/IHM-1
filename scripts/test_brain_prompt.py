"""exercise all six paths of POST /api/brain/prompt, and assert the honesty fields.

the point of these tests is not that the endpoint returns 200.  it is that a
response CANNOT reach a UI without carrying, in the payload:

  * the top-5 with cosine similarities and the corpus entry each one came from,
  * the chance level of the bank actually used,
  * an explicit statement that the output is a RETRIEVAL and not a generation,
  * the measured fact that similarity is not confidence -- the median top-1
    cosine falls only 0.706 -> 0.686 under a scramble that costs 63 points of
    accuracy,
  * and the measured accuracy of the path that produced it, so an audio result at
    6.4% cannot be presented as if it were a visual one at 63.5%.

so most of what follows asserts on the SHAPE of the payload.  a test that only
checked the retrieval worked would let the interface drop the caveats and still
pass, and the caveats are the part of this endpoint that is load-bearing.

these run against the real checkpoints and corpora in IBM-1, on CPU.  the first
test to touch each head pays its load (tens of seconds); the service is shared
across the class so that is paid once.

run:  .venv/bin/python -m unittest scripts.test_brain_prompt -v
"""
from __future__ import annotations

import base64
import io
import json
import unittest
import wave
from pathlib import Path

import numpy as np

from ihm.app.brain_prompt import (PATHS, UNSUPPORTED, BrainPrompt, describe_handle,
                                  handle, ibm_root)

ROOT = Path(__file__).resolve().parents[1]
IBM = ibm_root(ROOT)

REQUIRED = [IBM / 'ckpt/visual_contrastive_v2.pt',
            IBM / 'ckpt/audio_contrastive_v1.pt',
            IBM / 'data/derived/things-paired/images_test.npy',
            IBM / 'data/derived/libribrain-paired/row_index.npz',
            IBM / 'data/derived/libribrain-paired/cochleagram_scale.json',
            IBM / 'out/introspect/report.json',
            IBM / 'out/audio_prompt_bank.json']


@unittest.skipUnless(all(p.exists() for p in REQUIRED),
                     'IBM-1 materializations are not present: '
                     + ', '.join(str(p) for p in REQUIRED if not p.exists()))
class BrainPromptTest(unittest.TestCase):
    service = None

    @classmethod
    def setUpClass(cls):
        cls.service = BrainPrompt(ROOT)

    def call(self, body, expect=200):
        payload, status = handle(ROOT, body, service=self.service)
        self.assertEqual(status, expect, payload.get('error', payload))
        return payload

    # ------------------------------------------------------------ the honesty

    def assert_honest(self, p, out_kind):
        """what every response must carry, whichever path produced it."""
        self.assertEqual(p['schema'], 'ihm.brain-prompt.v1')
        self.assertEqual(p['output_kind'], out_kind)

        # 2. the output is a retrieval, said explicitly and not by implication
        self.assertIs(p['output_is_retrieval'], True)
        self.assertIs(p['generated'], False)
        self.assertIn('RETRIEVAL', p['output_is'])

        # 1. top-5, with cosine similarities, the chance level, and the entry
        self.assertEqual(len(p['results']), 5)
        self.assertEqual(p['retrieval']['top_k'], 5)
        self.assertGreater(p['retrieval']['n_bank_entries'], 1)
        self.assertAlmostEqual(p['retrieval']['chance'],
                               1.0 / p['retrieval']['n_bank_entries'], places=9)
        self.assertAlmostEqual(p['retrieval']['chance_pct'],
                               100 * p['retrieval']['chance'], places=9)
        sims = []
        for rank, r in enumerate(p['results'], 1):
            self.assertEqual(r['rank'], rank)
            self.assertIsInstance(r['cosine_similarity'], float)
            self.assertLessEqual(abs(r['cosine_similarity']), 1.0 + 1e-6)
            self.assertIn('corpus_entry', r)
            self.assertIn('corpus', r['corpus_entry'])
            sims.append(r['cosine_similarity'])
        self.assertEqual(sims, sorted(sims, reverse=True), 'results are not ranked')
        # both are rounded to 6 dp in the payload, so compare at 5
        self.assertAlmostEqual(p['top1_cosine'], sims[0], places=5)
        self.assertAlmostEqual(p['margin_over_rank2'], sims[0] - sims[1], places=5)

        # 3. similarity is not confidence, measured, in the payload
        c = p['similarity_is_not_confidence']
        self.assertIn('similarity is NOT confidence', c['headline'])
        self.assertAlmostEqual(c['intact_median_top1_cosine'], 0.7059, places=3)
        self.assertAlmostEqual(c['scrambled_median_top1_cosine'], 0.6860, places=3)
        self.assertAlmostEqual(c['accuracy_cost_of_that_scramble'], 0.63, places=2)
        self.assertLess(c['intact_median_top1_cosine']
                        - c['scrambled_median_top1_cosine'], 0.05,
                        'the whole point is that the cosine barely moves')
        self.assertIn('margin', c['use_instead'])

        # 4. the measured accuracy OF THE PATH USED
        m = p['measured']
        self.assertEqual(m['path'], p['path'])
        self.assertEqual(m['n_bank_entries'], p['retrieval']['n_bank_entries'])
        self.assertAlmostEqual(m['chance'], p['retrieval']['chance'], places=9)
        self.assertIsInstance(m['top1'], float)
        self.assertGreater(m['top1'], m['chance'])
        self.assertTrue(p['caveats'], 'a response with no caveats')

        # whether a correct answer is even present in the bank
        self.assertIn('query_target_in_bank', p['retrieval'])
        self.assertIsInstance(p['retrieval']['query_target_in_bank'], bool)
        if not p['retrieval']['query_target_in_bank']:
            self.assertIn('DOES NOT apply',
                          p['retrieval']['query_target_in_bank_means'])
        return p

    def assert_visual_accuracy(self, p):
        self.assertAlmostEqual(p['measured']['top1'], 0.635, places=3)
        self.assertEqual(p['retrieval']['n_bank_entries'], 200)
        self.assertAlmostEqual(p['retrieval']['chance'], 0.005, places=9)

    def assert_audio_accuracy(self, p):
        """the audio number must be present, weak, and named as weak."""
        m = p['measured']
        self.assertLess(m['top1'], 0.10)
        self.assertGreater(m['top1'], m['chance'])
        self.assertIn('top1_sd', m)
        self.assertEqual(m['top1_pools_averaged'], 8)
        # it does not beat its own dynamics-free control, and says so
        self.assertLess(m['top1'], m['dynamics_free_control']['top1'])
        self.assertIn('does NOT beat', m['verdict'])
        self.assertIn('not comparable', m['compare_to_vision'])
        self.assertLess(m['top1'], 0.2 * 0.635,
                        'the audio path is being reported as near-visual')
        self.assertTrue(any('WEAK' in c for c in p['caveats']))
        # a single pool of 200 has sd ~2.8 points; the served pool is labelled
        self.assertIn('single pool', m['served_bank_note'])

    # ------------------------------------------------------------- the six paths

    def test_1_image_to_image(self):
        p = self.call({'modality': 'image', 'want': 'image',
                       'payload': {'corpus_index': 0}})
        self.assert_honest(p, 'image')
        self.assert_visual_accuracy(p)
        self.assertEqual(p['path'], 'image->image')
        for r in p['results']:
            raw = base64.b64decode(r['image_png_base64'])
            self.assertEqual(raw[:8], b'\x89PNG\r\n\x1a\n')
            self.assertIn('concept', r['corpus_entry'])
            self.assertIn('image_path', r['corpus_entry'])
            self.assertIn('MEASURED evoked EEG', r['corpus_entry']['bank_is'])
        self.assertEqual(p['results'][0]['corpus_entry']['concept'],
                         'aircraft carrier')

    def test_2_image_to_text(self):
        p = self.call({'modality': 'image', 'want': 'text',
                       'payload': {'corpus_index': 0}})
        self.assert_honest(p, 'text')
        self.assert_visual_accuracy(p)
        self.assertEqual(p['results'][0]['text'], 'aircraft carrier')
        for r in p['results']:
            self.assertNotIn('image_png_base64', r)
        # the label is not an independent decode and the payload has to say so
        self.assertTrue(any('SAME NUMBER by construction' in c
                            for c in p['caveats']))

    def test_3_text_to_image(self):
        p = self.call({'modality': 'text', 'want': 'image',
                       'payload': {'text': 'aircraft carrier'}})
        self.assert_honest(p, 'image')
        self.assert_visual_accuracy(p)
        lookup = p['vocabulary_lookup']
        self.assertEqual(lookup['matched_concept'], 'aircraft carrier')
        self.assertEqual(lookup['match'], 'exact')
        self.assertEqual(base64.b64decode(lookup['image_png_base64'])[:8],
                         b'\x89PNG\r\n\x1a\n')
        # the exact half and the fallible half are reported separately
        self.assertIn('exact', p['measured']['vocabulary_lookup']['accuracy'])
        self.assertTrue(any('can miss' in c for c in p['caveats']))

    def test_4_audio_to_audio(self):
        q = int(self.service.audio().bank_rows[0])
        p = self.call({'modality': 'audio', 'want': 'audio',
                       'payload': {'corpus_index': q},
                       'resynthesise_top_k': 2})
        self.assert_honest(p, 'audio')
        self.assert_audio_accuracy(p)
        self.assertIn('MEASURED MEG', p['measured']['bank'])
        for r in p['results'][:2]:
            raw = base64.b64decode(r['audio_wav_base64'])
            self.assertEqual(raw[:4], b'RIFF')
            with wave.open(io.BytesIO(raw), 'rb') as w:
                self.assertEqual(w.getframerate(), 16000)
                self.assertEqual(w.getnchannels(), 1)
                self.assertGreater(w.getnframes(), 8000)
            # the resynthesis names itself, and names what it destroyed
            how = r['resynthesis']
            self.assertIn('Griffin-Lim', how['method'])
            self.assertIn('vocoding', how['this_is'])
            self.assertEqual(len(how['irreversible']), 2)
        for r in p['results'][2:]:
            self.assertNotIn('audio_wav_base64', r)
        # the retrieved entry is addressed in the corpus, not just scored
        e = p['results'][0]['corpus_entry']
        for key in ('run', 'chapter_time_s', 'stimulus_rows', 'meg_rows'):
            self.assertIn(key, e)

    def test_5_audio_to_text(self):
        q = int(self.service.audio().bank_rows[3])
        p = self.call({'modality': 'audio', 'want': 'text',
                       'payload': {'corpus_index': q}})
        self.assert_honest(p, 'text')
        self.assert_audio_accuracy(p)
        for r in p['results']:
            self.assertIsInstance(r['text'], str)
            self.assertNotIn('audio_wav_base64', r)
            self.assertIn('forced aligner', r['corpus_entry']['words_are'])
        self.assertTrue(any('FORCED ALIGNER' in c.upper() for c in p['caveats']))

    def test_6_text_to_audio(self):
        p = self.call({'modality': 'text', 'want': 'audio',
                       'payload': {'text': 'murder'}})
        self.assert_honest(p, 'audio')
        self.assert_audio_accuracy(p)
        lookup = p['vocabulary_lookup']
        self.assertEqual(lookup['matched_word'], 'murder')
        self.assertGreaterEqual(lookup['occurrences_in_corpus'], 1)
        self.assertIn('murder', [w.lower()
                                 for w in lookup['used']['words_in_window']])
        self.assertEqual(base64.b64decode(lookup['audio_wav_base64'])[:4], b'RIFF')
        self.assertIn('exact', p['measured']['vocabulary_lookup']['accuracy'])

    # ------------------------------------------------------ inputs and refusals

    def test_uploaded_image_is_flagged_off_distribution(self):
        from PIL import Image
        buf = io.BytesIO()
        Image.fromarray(np.full((96, 96, 3), 128, np.uint8)).save(buf, format='PNG')
        p = self.call({'modality': 'image', 'want': 'text',
                       'payload': {'image_base64':
                                   base64.b64encode(buf.getvalue()).decode()}})
        self.assert_honest(p, 'text')
        self.assertEqual(p['input']['kind'], 'uploaded_image')
        self.assertIn('does not transfer', p['input']['caveat'])

    def test_corpus_audio_query_is_scored_against_its_own_target(self):
        """a pool with no right answer in it cannot carry an accuracy figure."""
        ab = self.service.audio()
        q = int(ab.bank_rows[7]) + 4001          # deliberately NOT a bank row
        p = self.call({'modality': 'audio', 'want': 'text',
                       'payload': {'corpus_index': q}})
        self.assertIs(p['retrieval']['query_target_in_bank'], True)
        self.assertEqual(p['retrieval']['n_bank_entries'], 200)
        self.assertIn(q, [r['corpus_entry']['bank_row'] for r in p['results']]
                      + [q])                      # the target is in the pool
        rows, bank, in_bank = ab.pool_for(q)
        self.assertEqual(len(rows), 200)
        self.assertIn(q, rows)
        self.assertIs(in_bank, True)

    def test_uploaded_audio_has_no_target_in_the_bank(self):
        ab = self.service.audio()
        q = int(ab.bank_rows[0])
        x, _ = ab.resynthesise(ab.cochleagram_of(q), ab.entry(q)['chapter_file'])
        buf = io.BytesIO()
        with wave.open(buf, 'wb') as w:
            w.setnchannels(1); w.setsampwidth(2); w.setframerate(ab.FS)
            w.writeframes((np.clip(x, -1, 1) * 32767).astype('<i2').tobytes())
        p = self.call({'modality': 'audio', 'want': 'text',
                       'payload': {'audio_base64':
                                   base64.b64encode(buf.getvalue()).decode()}})
        self.assertIs(p['retrieval']['query_target_in_bank'], False)
        self.assertIn('DOES NOT apply',
                      p['retrieval']['query_target_in_bank_means'])

    def test_uploaded_audio_round_trips(self):
        """a resynthesised window, fed back in, is accepted and flagged."""
        ab = self.service.audio()
        q = int(ab.bank_rows[0])
        x, _ = ab.resynthesise(ab.cochleagram_of(q), ab.entry(q)['chapter_file'])
        buf = io.BytesIO()
        with wave.open(buf, 'wb') as w:
            w.setnchannels(1); w.setsampwidth(2); w.setframerate(ab.FS)
            w.writeframes((np.clip(x, -1, 1) * 32767).astype('<i2').tobytes())
        p = self.call({'modality': 'audio', 'want': 'text',
                       'payload': {'audio_base64':
                                   base64.b64encode(buf.getvalue()).decode()}})
        self.assert_honest(p, 'text')
        self.assertEqual(p['input']['kind'], 'uploaded_audio')
        self.assertIn('does not transfer', p['input']['caveat'])

    def test_resynthesis_round_trip_is_what_it_claims(self):
        """the number the payload quotes for the vocoder is the number measured."""
        ab = self.service.audio()
        rs = []
        for q in ab.bank_rows[:20]:
            q = int(q)
            coch = ab.cochleagram_of(q)
            x, how = ab.resynthesise(coch, ab.entry(q)['chapter_file'])
            buf = io.BytesIO()
            with wave.open(buf, 'wb') as w:
                w.setnchannels(1); w.setsampwidth(2); w.setframerate(ab.FS)
                w.writeframes((np.clip(x, -1, 1) * 32767).astype('<i2').tobytes())
            back, _ = ab.cochleagram_from_wav(buf.getvalue())
            rs.append(float(np.corrcoef(coch[-249:].ravel(),
                                        back[-249:].ravel())[0, 1]))
        quoted = how['measured_round_trip']
        self.assertEqual(quoted['n_windows'], len(rs))
        self.assertAlmostEqual(float(np.mean(rs)), quoted['pearson_r'], places=2)
        self.assertGreaterEqual(min(rs), quoted['min'] - 0.02)

    def test_training_split_window_is_labelled(self):
        ab = self.service.audio()
        p = self.call({'modality': 'audio', 'want': 'text',
                       'payload': {'corpus_index': ab.ctx + 10}})
        self.assertIn('TRAINING SPLIT', p['input']['split'])

    def test_unsupported_pairs_are_refused_not_faked(self):
        for pair in UNSUPPORTED:
            body = self.call({'modality': pair[0], 'want': pair[1], 'payload': {}},
                             expect=501)
            self.assertIn('why', body)
            self.assertIn('declared_paths', body)
            self.assertEqual(len(body['declared_paths']), 6)

    def test_word_outside_the_corpus_is_refused(self):
        body = self.call({'modality': 'text', 'want': 'audio',
                          'payload': {'text': 'zzzznotaword'}}, expect=404)
        self.assertIn('vocabulary_is', body)
        self.assertIn('cannot answer', body['vocabulary_is'])

    def test_concept_outside_the_corpus_is_refused(self):
        body = self.call({'modality': 'text', 'want': 'image',
                          'payload': {'text': 'zzzznotaconcept'}}, expect=404)
        self.assertIn('vocabulary_is', body)
        self.assertEqual(body['vocabulary_size'], 200)

    def test_bad_requests(self):
        for body, status in [({'modality': 'image'}, 400),
                             ({'modality': 'x', 'want': 'text'}, 400),
                             ({'modality': 'image', 'want': 'text',
                               'payload': {'corpus_index': 10 ** 6}}, 400),
                             ({'modality': 'image', 'want': 'text',
                               'payload': {'image_base64': 'not base64!!'}}, 400),
                             ({'modality': 'text', 'want': 'image',
                               'payload': {'text': 'x' * 300}}, 400)]:
            out, got = handle(ROOT, body, service=self.service)
            self.assertEqual(got, status, out)
            self.assertIn('error', out)

    # --------------------------------------------------------------- the catalog

    def test_catalog_lists_the_six_and_their_numbers(self):
        body, status = describe_handle(ROOT, service=self.service)
        self.assertEqual(status, 200)
        self.assertEqual({p['path'] for p in body['paths']},
                         {'%s->%s' % p for p in PATHS})
        for row in body['paths']:
            self.assertIn('top1', row)
            self.assertIn('chance', row)
        self.assertIs(body['output_is_retrieval'], True)
        self.assertIs(body['generated'], False)
        self.assertIn('does NOT beat', body['audio_is_much_weaker_than_vision'])
        self.assertEqual(len(body['refused']), 3)

    def test_every_payload_is_json_serialisable_by_the_server(self):
        """the workbench serialises with allow_nan=False; a NaN would 500."""
        for body in ({'modality': 'image', 'want': 'text',
                      'payload': {'corpus_index': 7}},
                     {'modality': 'audio', 'want': 'audio',
                      'payload': {'corpus_index':
                                  int(self.service.audio().bank_rows[1])}}):
            p = self.call(body)
            json.dumps(p, allow_nan=False, separators=(',', ':'))


if __name__ == '__main__':
    unittest.main()
