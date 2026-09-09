import test from 'node:test';
import assert from 'node:assert/strict';
import {readout, barFractions, wantsFor, STANDING_CAVEATS} from '../src/prompt-bar.js';

// The envelope the workbench actually serves.
const entry = (rank, cosine, concept) => ({
 rank, cosine_similarity: cosine, text: concept,
 corpus_entry: {bank_index: rank - 1, concept,
  bank_is: 'the group-mean MEASURED evoked EEG of this image, embedded by eeg_head -- never the image itself'},
});
const response = {
 schema: 'ihm.brain-prompt.v1', path: 'image->image', output_kind: 'image',
 output_is_retrieval: true, generated: false,
 output_is: 'a RETRIEVAL from the corpus.',
 similarity_is_not_confidence: {headline: 'similarity is NOT confidence',
  measured: 'the median top-1 cosine falls only 0.7059 -> 0.6860 under shuffle_sites.',
  consequence: 'a wrong answer looks almost exactly as certain as a right one.',
  use_instead: 'the margin over rank 2 separates them.'},
 retrieval: {n_bank_entries: 200, chance: 0.005, chance_pct: 0.5,
  ranked_by: 'cosine similarity between the cortical embedding and the embedded MEASURED neural response'},
 top1_cosine: 0.8007, margin_over_rank2: 0.0648,
 results: [entry(1, 0.8007, 'aircraft carrier'), entry(2, 0.7359, 'duffel bag'),
  entry(3, 0.6720, 'sailboat'), entry(4, 0.6340, 'cd player'),
  entry(5, 0.6290, 'chopsticks'), entry(6, 0.6000, 'sixth')],
};

test('a retrieval stays labelled a retrieval whatever the endpoint says', () => {
 const view = readout(response);
 assert.equal(view.retrieval, true);
 assert.equal(view.results.length, 5, 'top five, as introspect.py prints them');
 assert.equal(view.chance, 0.005);
 assert.equal(view.bank_entries, 200);
 assert.equal(view.modality, 'image');
 assert.equal(view.output_modality, 'image');
 assert.equal(view.results[0].label, 'aircraft carrier');
 assert.ok(Math.abs(view.margin - 0.0648) < 1e-9, 'the margin over rank 2 is the readable quantity');
 // The bank holds measured neural responses, never the images.
 assert.match(view.corpus, /MEASURED/);
 // A server that forgets the flag does not get to withdraw it: this path is a
 // nearest-entry lookup in the corpus and the label is not the server's.
 const quiet = readout({...response, output_is_retrieval: false});
 assert.equal(quiet.retrieval, true);
 assert.equal(quiet.declared_retrieval, false);
 // The caveat that a similarity is not a confidence is always carried.
 for (const line of STANDING_CAVEATS) assert.ok(view.caveats.includes(line));
 assert.ok(view.caveats.some(c => /not a confidence/i.test(c)));
});

test('a missing chance level is reported missing, never defaulted to a number', () => {
 const view = readout({results: response.results.slice(0, 2)});
 assert.equal(view.chance, null);
 assert.equal(view.results.length, 2);
 const broken = readout({results: [{label: 'x', cosine_similarity: 'not a number'}]});
 assert.equal(broken.results[0].similarity, null);
 assert.equal(broken.margin, null);
 assert.equal(readout(null).error, 'The prompt endpoint returned nothing readable.');
 assert.equal(readout({results: []}).empty, true);
});

test('a refusal keeps what the endpoint said about the corpus', () => {
 const refused = readout({error: 'image -> audio is not one of the six declared paths',
  why: 'the two corpora share no entry', declared_paths: ['image->image', 'image->text']});
 assert.match(refused.error, /six declared paths/);
 assert.match(refused.why, /share no entry/);
 assert.deepEqual(refused.declared_paths, ['image->image', 'image->text']);
 const missing = readout({error: 'no such concept', example_concepts: ['anchor', 'apple']});
 assert.deepEqual(missing.candidates, ['anchor', 'apple']);
});

test('only the six declared paths are offered', () => {
 // image->audio and audio->image have no materialization joining the corpora,
 // and text->text is a lookup returning its own word. None is offered.
 assert.deepEqual(wantsFor('image'), ['image', 'text']);
 assert.deepEqual(wantsFor('audio'), ['audio', 'text']);
 assert.deepEqual(wantsFor('text'), ['image', 'audio']);
 assert.deepEqual(wantsFor('nonsense'), []);
});

test('similarity bars are scaled to the top result, not to a probability', () => {
 const view = readout(response);
 const fractions = barFractions(view.results);
 assert.equal(fractions[0], 1);
 assert.ok(Math.abs(fractions[1] - 0.7359 / 0.8007) < 1e-9);
 // Every bar is a fraction of the widest cosine drawn, so none of them can be
 // read as "82% sure".
 assert.ok(fractions.every(f => f >= 0 && f <= 1));
 assert.deepEqual(barFractions([{similarity: null}]), [0]);
});
