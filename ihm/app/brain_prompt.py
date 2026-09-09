"""POST /api/brain/prompt -- image, audio or text in; image, audio or text out.

WHAT THIS IS.  six paths, each one a materialization that already exists in IBM-1,
wired to an HTTP route.  nothing here is new science and nothing here calls out of
the process: no integrations, no tools, no retrieval from outside the corpus.

    image -> image   drive the optic nerve, run the cortical dynamics, embed the
    image -> text    state, match it against the bank of MEASURED THINGS-EEG2
                     evoked responses, and return what the matched entries carry
                     -- their image, and the concept name in their path.
                     `scripts/introspect.py` in IBM-1 is the reference; this is
                     that decode, served.
    text  -> image   look the concept up in the corpus vocabulary (exact, no
    text  -> audio   model), then drive the cortex with the entry it names and
                     report what comes back, so the retrieval is visible and not
                     hidden behind a lookup that cannot fail.
    audio -> audio   cochleagram in, drive the temporal port, run, embed, match
    audio -> text    against the bank of MEASURED LibriBrain MEG windows, and
                     return the cochleagram paired to the match -- resynthesised
                     to a waveform -- or the words the aligner puts in its span.

EVERY OUTPUT IS A RETRIEVAL.  the brain does not draw an image and does not speak.
it finds the corpus entry its cortical state is nearest to, and that entry's
image, words or sound is what comes back.  the response says so in a field
(`output_is_retrieval`, `generated`) rather than leaving it to the caller's
caption, because showing a returned image without that label would be the single
most misleading thing this app could do.

SIMILARITY IS NOT CONFIDENCE, and this is measured, not asserted.  on the
designated visual test set the median top-1 cosine is 0.7059 for the intact
cortical state and 0.6860 when the read sites are permuted -- a 0.020 fall for an
ablation that takes accuracy from 63.50% to 0.50%.  a wrong answer looks almost
exactly as certain as a right one.  the margin over rank 2 is the quantity that
separates them, so both travel in every response, in
`similarity_is_not_confidence`, and the caller cannot render a result without
having been handed the disclaimer.

THE AUDIO PATHS ARE MUCH WEAKER THAN THE VISUAL ONES, and the response says how
much.  re-measured at the checkpoint over 8 held-out pools of 200
(`scripts/eval_audio_prompt_bank.py`), the auditory head reaches 6.37% +/- 1.45
against a measured dynamics-free control ceiling of 7.12% +/- 1.71 -- it does not
beat its own control.  the visual head reaches 63.50%.  `measured` carries the
accuracy of the path ACTUALLY USED, so an audio result cannot be dressed as a
visual one.
"""
from __future__ import annotations

import base64
import io
import json
import os
import threading
import wave
from pathlib import Path

SCHEMA = 'ihm.brain-prompt.v1'

#: the six declared paths.  anything else is refused rather than approximated.
PATHS = {('image', 'image'), ('image', 'text'), ('text', 'image'),
         ('text', 'audio'), ('audio', 'audio'), ('audio', 'text')}
UNSUPPORTED = {
    ('image', 'audio'): 'no materialization joins the visual corpus to the auditory one: '
                        'THINGS-EEG2 carries no sound and LibriBrain carries no images, so '
                        'an image cannot be turned into audio without inventing the link.',
    ('audio', 'image'): 'the converse of the same gap; the two corpora share no entry.',
    ('text', 'text'): 'a vocabulary lookup returning the word it was given is not a decode, '
                      'and there is no measurement that would make it one.',
}

TOP_K = 5


def ibm_root(root):
    """where IBM-1 is.  the checkpoints and corpora live there, not here."""
    override = os.environ.get('IBM_ROOT')
    return Path(override).resolve() if override else Path(root).resolve().parent / 'IBM-1'


class PromptError(Exception):
    def __init__(self, message, status=400, detail=None):
        super().__init__(message)
        self.status = status
        self.detail = detail or {}


# --------------------------------------------------------------------- signals

def _griffin_lim(mag, n_fft, hop, iters, seed=0):
    """phase from magnitude, by iterated STFT/ISTFT.  mag: (frames, bins)."""
    import numpy as np
    window = np.hanning(n_fft).astype(np.float32)
    n = (len(mag) - 1) * hop + n_fft
    rng = np.random.default_rng(seed)
    angle = np.exp(2j * np.pi * rng.random(mag.shape)).astype(np.complex64)

    def istft(spec):
        frames = np.fft.irfft(spec, n=n_fft, axis=1).astype(np.float32) * window
        out = np.zeros(n, np.float32)
        norm = np.zeros(n, np.float32)
        for i in range(len(frames)):
            out[i * hop:i * hop + n_fft] += frames[i]
            norm[i * hop:i * hop + n_fft] += window ** 2
        return out / np.maximum(norm, 1e-8)

    def stft(x):
        idx = np.arange(n_fft)[None, :] + hop * np.arange(len(mag))[:, None]
        pad = np.zeros(n, np.float32)
        pad[:min(len(x), n)] = x[:n]
        return np.fft.rfft(pad[idx] * window[None, :], axis=1)

    x = istft(mag * angle)
    for _ in range(iters):
        s = stft(x)
        x = istft(mag * np.exp(1j * np.angle(s)))
    peak = float(np.abs(x).max())
    return x / peak if peak > 0 else x


def _wav_base64(x, rate):
    import numpy as np
    pcm = np.clip(x, -1.0, 1.0)
    pcm = (pcm * 32767.0).astype('<i2')
    buf = io.BytesIO()
    with wave.open(buf, 'wb') as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(rate)
        w.writeframes(pcm.tobytes())
    return base64.b64encode(buf.getvalue()).decode()


def _png_base64(arr):
    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format='PNG')
    return base64.b64encode(buf.getvalue()).decode()


def _decode_b64(value, field, limit=8 << 20):
    if not isinstance(value, str):
        raise PromptError(f'{field} must be a base64 string')
    try:
        raw = base64.b64decode(value, validate=True)
    except Exception:
        raise PromptError(f'{field} is not valid base64')
    if not raw or len(raw) > limit:
        raise PromptError(f'{field} must be 1-{limit} bytes once decoded')
    return raw


# ------------------------------------------------------------ the visual brain

class VisualBrain:
    """`ckpt/visual_contrastive_v2.pt`, loaded once, as introspect.py loads it.

    the normalisation of the evoked bank is fitted on the TRAINING split only.
    fitting it on the test set leaks the test set into the query, and the number
    this endpoint reports is the designated-test-set number, so it has to be the
    same normalisation the measurement used.
    """

    ONSET, KEEP = 20, 50
    TEMP = 0.07                     # the training INITIALISATION; see caveats

    def __init__(self, ibm):
        import importlib.util
        import numpy as np
        import torch

        self.np, self.torch = np, torch
        self.ibm = ibm
        d = ibm / 'data/derived/things-paired'
        ckpt = ibm / 'ckpt/visual_contrastive_v2.pt'
        for p in (ckpt, d / 'images_test.npy', d / 'evoked_test_groupmean.npy',
                  d / 'image_paths_test.npy', d / 'evoked_training_groupmean.npy'):
            if not p.exists():
                raise PromptError(f'visual materialization unavailable: {p} is missing',
                                  503)

        spec = importlib.util.spec_from_file_location(
            'ibm_pretrain_video_loop', str(ibm / 'scripts/pretrain_video_loop.py'))
        P = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(P)

        blob = torch.load(str(ckpt), map_location='cpu', weights_only=False)
        cfg, sd = blob['config'], blob['model']
        n_sites, embed = sd['dyn.embed'].shape
        k = sd['dyn.idx'].shape[1]

        tr = np.load(str(d / 'evoked_training_groupmean.npy'), mmap_mode='r')
        ntr = int(min(len(tr), len(np.load(str(d / 'images_training.npy'), mmap_mode='r')))
                  * 0.8)
        samp = np.asarray(tr[:ntr:7]).astype(np.float32)
        med = np.median(samp, 0)
        iqr = ((np.percentile(samp, 75, 0) - np.percentile(samp, 25, 0)) / 1.349).clip(1e-9)

        self.images = np.load(str(d / 'images_test.npy'))
        self.paths = np.load(str(d / 'image_paths_test.npy'))
        ev = np.load(str(d / 'evoked_test_groupmean.npy'))
        dec = cfg['decimate']
        y = np.clip((ev.astype(np.float32) - med) / iqr, -6, 6)[..., self.ONSET:self.ONSET + self.KEEP:dec]

        dyn = P.CorticalDynamics(n_sites, embed, k, 'cpu', long_range=cfg['long_range'])
        self.model = P.VisualContrastiveLoop(dyn, n_sensors=ev.shape[1], n_times=y.shape[-1])
        self.model.load_state_dict(sd)
        self.model.eval()
        self.dyn, self.cfg, self.step = dyn, cfg, blob['step']
        self.readout = cfg.get('n_steps', 4) * cfg['substeps']
        self.h = cfg['dt'] / cfg['substeps']

        with torch.no_grad():
            self.bank = self.model.embed_eeg(torch.from_numpy(y))   # MEASURED EEG
        self.concepts = [self._concept(p) for p in self.paths]
        self.chance = 1.0 / len(self.images)

    @staticmethod
    def _concept(p):
        """'.../00001_aircraft_carrier/aircraft_carrier_06s.jpg' -> 'aircraft carrier'.

        read from the directory the dataset itself names the concept in, never
        reconstructed from an ordering that merely has the right length.
        """
        base = os.path.basename(os.path.dirname(str(p)))
        return base.split('_', 1)[1].replace('_', ' ') if '_' in base else base

    def embed_image(self, img_uint8):
        """(H, W, 3) uint8 -> the 128-d cortical embedding at the trained readout."""
        torch, np = self.torch, self.np
        x = torch.from_numpy(np.array(img_uint8[None], copy=True))
        x = (x.permute(0, 3, 1, 2).float() / 127.5) - 1
        with torch.no_grad():
            drive = torch.zeros(1, self.dyn.n)
            drive[:, :self.model.port] = self.model.to_cortex(self.model.enc(x))
            s = self.dyn.init_state(1, x.device)
            w = self.dyn.edge_weights()
            for _ in range(self.readout):
                s = self.dyn.step(s, drive, self.h, w)
            r = s[1][:, self.model.read_idx]
            return torch.nn.functional.normalize(self.model.cortex_head(r), dim=-1)

    def retrieve(self, z, k=TOP_K):
        torch = self.torch
        sim = (z @ self.bank.T)[0]
        order = torch.argsort(sim, descending=True)
        top = sim[order]
        return ([int(v) for v in order[:k]], [float(v) for v in top[:k]],
                float(top[0]), float(top[0] - top[1]))

    def image_shape(self):
        return tuple(int(v) for v in self.images.shape[1:])


# ------------------------------------------------------------- the audio brain

class AudioBrain:
    """`ckpt/audio_contrastive_v1.pt` over `data/derived/libribrain-paired` (v3).

    the retrieval bank is the MEASURED MEG of held-out windows, exactly as the
    visual bank is measured EEG.  a bank of cochleagram embeddings would make the
    query a member of its own corpus and the "decode" an identity lookup that
    returned the input, which is the failure introspect.py names first.

    the bank is ONE pool of 200 drawn from the held-out tail with a fixed seed, so
    the chance level the response quotes (0.5%) is the chance level of the bank in
    use.  its single-pool accuracy is reported too, and labelled: a single pool of
    200 has sd ~2.8 points, so the 8-pool mean is the number to read.
    """

    SEED = 20260909
    GAP = 2500                       # the builder's 10 s guard band
    FS, RATE, NFFT = 16000, 250, 1024

    def __init__(self, ibm):
        import importlib.util
        import numpy as np
        import torch

        self.np, self.torch = np, torch
        self.ibm = ibm
        d = ibm / 'data/derived/libribrain-paired'
        ckpt = ibm / 'ckpt/audio_contrastive_v1.pt'
        need = [ckpt, d / 'cochleagram_250hz.npy', d / 'meg_250hz.npy',
                d / 'row_index.npz', d / 'row_index.json', d / 'row_words.json',
                d / 'cochleagram_scale.json']
        for p in need:
            if not p.exists():
                raise PromptError(
                    f'auditory materialization unavailable: {p} is missing.  '
                    'row_index/row_words/cochleagram_scale come from '
                    'scripts/build_libribrain_index.py and '
                    'scripts/build_cochleagram_scale.py in IBM-1.', 503)

        spec = importlib.util.spec_from_file_location(
            'ibm_pretrain_video_loop_a', str(ibm / 'scripts/pretrain_video_loop.py'))
        P = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(P)

        blob = torch.load(str(ckpt), map_location='cpu', weights_only=False)
        cfg, sd = blob['config'], blob['model']
        self.cfg, self.step = cfg, blob['step']
        self.win, self.ctx, self.pool = cfg['win'], cfg['ctx'], cfg['pool']

        self.stim = np.load(str(d / 'cochleagram_250hz.npy'), mmap_mode='r')
        self.meg = np.load(str(d / 'meg_250hz.npy'), mmap_mode='r')
        self.n = min(len(self.stim), len(self.meg))
        self.ntr = int(self.n * 0.8)

        idx = np.load(str(d / 'row_index.npz'))
        self.row_run, self.row_t = idx['run'], idx['chapter_time']
        self.index_meta = json.loads((d / 'row_index.json').read_text())
        self.words = json.loads((d / 'row_words.json').read_text())
        self.scale = json.loads((d / 'cochleagram_scale.json').read_text())
        self.run_stems = [r['stem'] for r in self.index_meta['runs']]
        self.run_chapter = {r['stem']: r['chapter'] for r in self.index_meta['runs']}
        # word -> the rows it is spoken over, for the text -> audio lookup
        self.vocab = {}
        for ri, stem in enumerate(self.run_stems):
            for w in self.words[stem]:
                self.vocab.setdefault(w['word'].lower(), []).append((ri, w['t'], w['dur']))

        n_sites, embed = sd['dyn.embed'].shape
        dyn = P.CorticalDynamics(n_sites, embed, sd['dyn.idx'].shape[1], 'cpu',
                                 long_range=cfg['long_range'])
        self.model = P.AudioContrastiveLoop(dyn, n_sensors=self.meg.shape[1],
                                            n_bands=self.stim.shape[1],
                                            stim_len=self.ctx, n_times=self.win)
        self.model.load_state_dict(sd)
        self.model.eval()
        self.dyn = dyn

        lo, hi = self.ntr + self.GAP, self.n - self.win - 1
        self.bank_rows = [int(v) for v in
                          np.random.default_rng(self.SEED).integers(lo, hi, self.pool)]
        with torch.no_grad():
            y = np.stack([np.clip(np.asarray(self.meg[q:q + self.win]).astype(np.float32),
                                  -6, 6).T for q in self.bank_rows])
            self.bank = self.model.embed_meg(torch.from_numpy(y))   # MEASURED MEG
        self.chance = 1.0 / self.pool

        # the ERB filterbank the corpus was built with, rebuilt from its declaration
        f = np.fft.rfftfreq(self.NFFT, 1 / self.FS)
        cf = np.array(self.scale['erb_centres_hz'], np.float64)
        W = np.exp(-0.5 * ((f[None, :] - cf[:, None]) / (0.12 * cf[:, None])) ** 2)
        self.W = (W / W.sum(1, keepdims=True)).astype(np.float32)

    # ---- corpus addressing

    def row_span(self, q):
        """the cochleagram window paired to bank row q: stim[q-ctx : q].

        this is the pairing `train_audio_contrastive.py` trains: the stimulus
        window ENDS where the MEG window begins, because the cortical response
        lags the sound.  reproduced exactly rather than re-centred.
        """
        return q - self.ctx, q

    def entry(self, q):
        a, b = self.row_span(q)
        run = int(self.row_run[a])
        stem = self.run_stems[run]
        t0, t1 = float(self.row_t[a]), float(self.row_t[b - 1])
        spoken = [w['word'] for w in self.words[stem]
                  if t0 <= w['t'] <= t1 or t0 <= w['t'] + w['dur'] <= t1]
        return {
            'bank_row': int(q),
            'corpus': 'data/derived/libribrain-paired (v3)',
            'source': 'LibriBrain: one participant, MEG at 250 Hz while listening '
                      'to "A Study in Scarlet"',
            'run': stem,
            'chapter_file': self.run_chapter[stem],
            'chapter_time_s': [round(t0, 3), round(t1, 3)],
            'stimulus_rows': [int(a), int(b)],
            'meg_rows': [int(q), int(q + self.win)],
            'words_in_window': spoken,
            'words_are': 'the forced aligner\'s output for this run, read from the '
                         'dataset\'s own events file.  a word onset here is where an '
                         'aligner put it, not when a phoneme occurred; the source '
                         'card grades the alignments stream as context for that reason.',
        }

    def cochleagram_of(self, q):
        a, b = self.row_span(q)
        return self.np.asarray(self.stim[a:b]).astype('float32')

    # ---- the model

    def embed(self, coch):
        """(ctx, bands) -> the 128-d cortical embedding, temporal port driven."""
        torch = self.torch
        x = torch.from_numpy(coch[None]).permute(0, 2, 1)
        with torch.no_grad():
            z, _ = self.model.embed_audio(x, substeps=self.cfg['substeps'],
                                          dt=self.cfg['dt'],
                                          n_steps=self.cfg['n_steps'])
        return z

    def embed_meg_row(self, q):
        np, torch = self.np, self.torch
        y = np.clip(np.asarray(self.meg[q:q + self.win]).astype(np.float32), -6, 6).T
        with torch.no_grad():
            return self.model.embed_meg(torch.from_numpy(y[None]))

    def pool_for(self, q):
        """the 200-way pool this query is scored against.

        the measured task is a 200-way match in which the query's own MEG is one
        of the candidates.  a fixed bank that does not contain the query's target
        is a DIFFERENT and easier-looking task -- nothing in it can be right, so
        no ranking from it can be wrong either, and quoting 6.37% beside it would
        be quoting an accuracy for a question that was not asked.

        so a corpus query is scored against 199 fixed held-out distractors plus
        its own paired MEG.  the pool size, and therefore the chance level, is the
        200 the measurement used.  an UPLOADED clip has no paired MEG in this
        corpus and gets the fixed 200 with `query_target_in_bank` false.
        """
        if q is None or int(q) in set(self.bank_rows):
            return self.bank_rows, self.bank, q is not None
        rows = self.bank_rows[:self.pool - 1] + [int(q)]
        bank = self.torch.cat([self.bank[:self.pool - 1], self.embed_meg_row(int(q))])
        return rows, bank, True

    def retrieve(self, z, rows, bank, k=TOP_K):
        torch = self.torch
        sim = (z @ bank.T)[0]
        order = torch.argsort(sim, descending=True)
        top = sim[order]
        return ([rows[int(v)] for v in order[:k]],
                [float(v) for v in top[:k]],
                float(top[0]), float(top[0] - top[1]))

    # ---- resynthesis

    def resynthesise(self, coch, chapter, iters=48):
        """a stored cochleagram window back to a waveform.

        the chain, and what each link costs:

          1. undo the per-chapter standardisation.  EXACT -- the two scalars the
             corpus builder divided out and did not save were recomputed by
             `scripts/build_cochleagram_scale.py` and are read back here.
          2. expm1(a)/20 -> the 64 ERB band magnitudes.  exact.
          3. 64 bands -> 513 FFT bins through the filterbank's normalised
             transpose.  UNDERDETERMINED: the forward map threw away 449
             dimensions per frame and no inverse can put them back.  this
             spreads each band over its own filter shape, which is a choice.
          4. Griffin-Lim, 48 iterations, for the phase.  the STFT phase was
             never stored, so it is INVENTED here -- consistent with the
             magnitudes, but not the phase the microphone recorded.

        so the waveform is a vocoding of the retrieved corpus entry, and the
        response says so.  it is not the stimulus audio and it is not generated
        speech.
        """
        np = self.np
        sc = self.scale['chapters'][chapter]
        a = coch.astype(np.float64) * sc['std'] + sc['mean']
        bands = np.expm1(np.clip(a, 0.0, 30.0)) / 20.0
        spread = self.W.sum(0)
        mag = (bands @ self.W) / np.maximum(spread, 1e-6)
        x = _griffin_lim(mag.astype(np.float32), self.NFFT, self.FS // self.RATE, iters)
        return x, {
            'method': 'inverse ERB filterbank (normalised transpose) + Griffin-Lim',
            'griffin_lim_iters': iters,
            'sample_rate_hz': self.FS,
            'standardisation': 'exact; per-chapter mean and sd recomputed by '
                               'scripts/build_cochleagram_scale.py',
            'irreversible': ['the 513 -> 64 band projection: 449 spectral '
                             'dimensions per frame were discarded at corpus build '
                             'and no inversion recovers them',
                             'the STFT phase was never stored; Griffin-Lim invents '
                             'a phase consistent with the magnitudes, not the '
                             'recorded one'],
            'peak_normalised': True,
            'measured_round_trip': {
                'what': 'cochleagram the resynthesis and correlate it against the '
                        'cochleagram it came from, over the 20 bank windows',
                'pearson_r': 0.9163, 'sd': 0.0370, 'min': 0.8223, 'n_windows': 20,
                'asserted_by': 'scripts/test_brain_prompt.py',
                'read_it_as': 'the SPECTRAL ENVELOPE survives the round trip.  it '
                              'says nothing about intelligibility, and it is a '
                              'self-consistency check on the vocoder, not a '
                              'measure of the retrieval.'},
            'this_is': 'a vocoding of a RETRIEVED corpus cochleagram, not '
                       'generated speech and not the stimulus recording',
        }

    def cochleagram_from_wav(self, raw):
        """an uploaded mono/stereo PCM WAV -> a (ctx, 64) window in corpus units."""
        np = self.np
        with wave.open(io.BytesIO(raw), 'rb') as w:
            if w.getsampwidth() != 2:
                raise PromptError('audio upload must be 16-bit PCM WAV')
            frames = w.readframes(w.getnframes())
            ch, rate = w.getnchannels(), w.getframerate()
        x = np.frombuffer(frames, '<i2').astype(np.float32) / 32768.0
        if ch > 1:
            x = x.reshape(-1, ch).mean(1)
        if rate != self.FS:                      # linear resample onto 16 kHz
            m = int(round(len(x) * self.FS / rate))
            x = np.interp(np.linspace(0, len(x) - 1, m), np.arange(len(x)), x
                          ).astype(np.float32)
        hop = self.FS // self.RATE
        n = (len(x) - self.NFFT) // hop
        if n < 50:
            raise PromptError(
                f'audio upload is too short: {len(x) / self.FS:.2f} s gives '
                f'{max(n, 0)} cochleagram frames and the head reads a '
                f'{self.ctx / self.RATE:.1f} s window.  at least 0.27 s is needed '
                f'before the window is mostly the padding.')
        idx = np.arange(self.NFFT)[None, :] + hop * np.arange(n)[:, None]
        S = np.abs(np.fft.rfft(x[idx] * np.hanning(self.NFFT)[None, :], axis=1))
        a = np.log1p((S @ self.W.T) * 20.0).astype(np.float32)
        a = (a - a.mean()) / (a.std() + 1e-6)
        pad = max(self.ctx - n, 0)
        if pad:
            # the corpus window is a fixed 250 frames.  a shorter clip is held at
            # the RIGHT of the window and the head of it is silence, because the
            # window ends where the MEG the head is matched against begins.
            a = np.concatenate([np.full((pad, a.shape[1]), a.min(), np.float32), a])
        return a[-self.ctx:], {
            'frames': int(n), 'used': 'the last %d frames (%.1f s)' % (
                self.ctx, self.ctx / self.RATE),
            'silence_padded_frames': int(pad),
            'padding_note': ('%d of the %d frames the head reads are silence '
                             'padding, not sound you supplied' % (pad, self.ctx))
                            if pad else 'none',
            'standardisation': 'per-clip mean and sd, as the corpus builder '
                               'standardises per chapter.  a short clip\'s two '
                               'moments are not a chapter\'s, so an uploaded clip '
                               'sits at a slightly different operating point from '
                               'the corpus the head was fitted on.',
        }


# ------------------------------------------------------------------ the service

class BrainPrompt:
    """loads the two heads lazily and holds them; each load is seconds, not ms."""

    def __init__(self, root):
        self.root = Path(root)
        self.ibm = ibm_root(root)
        self.lock = threading.Lock()
        self._visual = None
        self._audio = None
        self._measured = None

    def measured(self):
        if self._measured is None:
            vis = self.ibm / 'out/introspect/report.json'
            aud = self.ibm / 'out/audio_prompt_bank.json'
            for p, script in ((vis, 'scripts/introspect.py'),
                              (aud, 'scripts/eval_audio_prompt_bank.py')):
                if not p.exists():
                    raise PromptError(
                        f'no measurement to report: {p} is missing.  run {script} '
                        'in IBM-1.  this endpoint does not serve a result without '
                        'the accuracy of the path that produced it.', 503)
            self._measured = (json.loads(vis.read_text()), json.loads(aud.read_text()))
        return self._measured

    def visual(self):
        with self.lock:
            if self._visual is None:
                self._visual = VisualBrain(self.ibm)
            return self._visual

    def audio(self):
        with self.lock:
            if self._audio is None:
                self._audio = AudioBrain(self.ibm)
            return self._audio


# ---------------------------------------------------------------- the payloads

def _visual_measured(report, path):
    full = report['ablations']['full']
    scr = report['ablations']['shuffle_sites']
    return {
        'path': path,
        'materialization': 'ckpt/visual_contrastive_v2.pt (IBM-1), train step %d'
                           % report['train_step'],
        'set': report['set'],
        'bank': report['retrieval_bank'],
        'n_bank_entries': report['n'],
        'chance': report['chance'],
        'top1': full['top1'],
        'top5': full['top5'],
        'x_chance': report['x_chance'],
        'median_top1_cosine': full['median_top1_sim'],
        'median_margin_over_rank2': report['median_margin_over_rank2'],
        'median_margin_when_correct': report['median_margin_when_correct'],
        'median_margin_when_wrong': report['median_margin_when_wrong'],
        'dynamics_free_control_top1_same_set': report['dynamics_free_control_top1_same_set'],
        'control_note': 'a separately-trained control with no dynamics reaches '
                        '66.50% on this same set against this model\'s 63.50%, a '
                        '0.89 sd difference.  the two are indistinguishable; this '
                        'is a model that HAS a cortex, not evidence for having one.',
        'ablations': {k: {'top1': v['top1'], 'median_top1_cosine': v['median_top1_sim'],
                          'median_margin': v['median_margin']}
                      for k, v in report['ablations'].items()},
        'scramble_control': {
            'arm': 'shuffle_sites (the same state values, the wrong topography)',
            'top1': scr['top1'],
            'median_top1_cosine': scr['median_top1_sim'],
        },
    }


def _audio_measured(report, path, visual_top1=None):
    return {
        'path': path,
        'materialization': 'ckpt/audio_contrastive_v1.pt (IBM-1), train step %d'
                           % report['train_step'],
        'corpus': report['corpus'],
        'split': report['split'],
        'bank': 'the MEASURED MEG of %d held-out 1 s windows, embedded by meg_head. '
                'the query is the cortical state; the bank is never the cochleagrams.'
                % report['pool'],
        'n_bank_entries': report['pool'],
        'chance': report['chance'],
        'top1': report['top1_mean_8_pools'],
        'top1_sd': report['top1_sd_8_pools'],
        'top1_pools_averaged': 8,
        'x_chance': report['top1_x_chance'],
        'served_bank_top1_single_pool': report['served_bank_top1_single_pool'],
        'served_bank_note': report['served_bank_note'],
        'pool_composition': 'a corpus query is scored against 199 fixed held-out '
                            'distractors plus its OWN paired MEG, so the pool size '
                            'and chance level are the measurement\'s.  the '
                            'measurement itself scored 200 windows against each '
                            'other simultaneously; this is one query against 200 '
                            'candidates.  same chance, same difficulty class, not '
                            'the identical procedure.',
        'dynamics_free_control': report['dynamics_free_control'],
        'verdict': report['verdict'],
        'compare_to_vision': 'the visual paths reach %.2f%% top-1 on a 200-way set '
                             'with the same 0.5%% chance level.  this path reaches '
                             '%.2f%%.  they are not comparable results and must not '
                             'be shown as if they were.'
                             % (100 * (visual_top1 if visual_top1 is not None else 0.635),
                                100 * report['top1_mean_8_pools']),
    }


def _similarity_caveat(report):
    full = report['ablations']['full']
    scr = report['ablations']['shuffle_sites']
    return {
        'headline': 'similarity is NOT confidence',
        'measured': 'the median top-1 cosine falls only %.4f -> %.4f under the '
                    'shuffle_sites ablation, which costs %.0f points of accuracy '
                    '(%.2f%% -> %.2f%%).'
                    % (full['median_top1_sim'], scr['median_top1_sim'],
                       100 * (full['top1'] - scr['top1']),
                       100 * full['top1'], 100 * scr['top1']),
        'consequence': 'a wrong answer looks almost exactly as certain as a right '
                       'one.  do not render the cosine as a confidence bar.',
        'use_instead': 'the margin over rank 2 separates them: %.4f when the '
                       'retrieval is correct against %.4f when it is wrong.'
                       % (report['median_margin_when_correct'],
                          report['median_margin_when_wrong']),
        'intact_median_top1_cosine': full['median_top1_sim'],
        'scrambled_median_top1_cosine': scr['median_top1_sim'],
        'accuracy_cost_of_that_scramble': full['top1'] - scr['top1'],
    }


IN_BANK = {
    True: 'the neural response paired to this query IS one of the bank entries, so '
          'a correct answer exists in the bank and the measured accuracy applies '
          'as stated.',
    False: 'the neural response paired to this query is NOT in the bank.  no bank '
           'entry is the right answer, so the ranking below is the nearest of the '
           'bank to a query the bank does not contain, and the measured accuracy '
           'DOES NOT apply to it -- there is nothing here for it to be right '
           'about.  read the top-5 as neighbours, not as candidates.',
}


def _envelope(path, out_kind, results, chance, n_bank, measured, sim_caveat, extra,
              target_in_bank=True):
    payload = {
        'schema': SCHEMA,
        'path': path,
        'output_kind': out_kind,
        # -- the two fields the UI is not permitted to omit
        'output_is_retrieval': True,
        'generated': False,
        'output_is': 'a RETRIEVAL from the corpus.  the brain did not draw, write '
                     'or speak this.  it found the corpus entry its cortical state '
                     'is nearest to, and what you are looking at is that entry.',
        'similarity_is_not_confidence': sim_caveat,
        'measured': measured,
        'retrieval': {
            'n_bank_entries': n_bank,
            'chance': chance,
            'query_target_in_bank': bool(target_in_bank),
            'query_target_in_bank_means': IN_BANK[bool(target_in_bank)],
            'chance_pct': 100 * chance,
            'ranked_by': 'cosine similarity between the cortical embedding and the '
                         'embedded MEASURED neural response of each bank entry',
            'top_k': len(results),
        },
        'results': results,
    }
    payload.update(extra)
    return payload


# ------------------------------------------------------------------- the paths

def _visual_results(vb, order, sims, want_image, k=TOP_K):
    import numpy as np
    torch = vb.torch
    logits = torch.tensor(sims) / VisualBrain.TEMP
    probs = torch.softmax(logits, 0)
    out = []
    for rank, (i, s) in enumerate(zip(order, sims), 1):
        entry = {
            'rank': rank,
            'cosine_similarity': round(float(s), 6),
            'softmax_over_top5_at_tau_0.07': round(float(probs[rank - 1]), 6),
            'corpus_entry': {
                'bank_index': int(i),
                'corpus': 'THINGS-EEG2 designated test set (200 images, 80 '
                          'repetitions each), data/derived/things-paired',
                'bank_is': 'the group-mean MEASURED evoked EEG of this image, '
                           'embedded by eeg_head -- never the image itself',
                'concept': vb.concepts[i],
                'image_path': str(vb.paths[i]),
            },
            'text': vb.concepts[i],
        }
        if want_image:
            entry['image_png_base64'] = _png_base64(vb.images[i])
            entry['image_shape'] = list(vb.image_shape())
        out.append(entry)
    return out


def _audio_results(ab, rows, sims, resynth_k):
    out = []
    for rank, (q, s) in enumerate(zip(rows, sims), 1):
        entry = ab.entry(q)
        item = {
            'rank': rank,
            'cosine_similarity': round(float(s), 6),
            'corpus_entry': entry,
            'text': ' '.join(entry['words_in_window']),
        }
        if rank <= resynth_k:
            coch = ab.cochleagram_of(q)
            x, how = ab.resynthesise(coch, entry['chapter_file'])
            item['audio_wav_base64'] = _wav_base64(x, ab.FS)
            item['audio_duration_s'] = round(len(x) / ab.FS, 4)
            item['resynthesis'] = how
        out.append(item)
    return out


def _lookup_concept(vb, text):
    key = ' '.join(str(text).lower().replace('_', ' ').split())
    if not key:
        raise PromptError('text payload is empty')
    exact = [i for i, c in enumerate(vb.concepts) if c.lower() == key]
    if exact:
        return exact[0], 'exact'
    sub = [i for i, c in enumerate(vb.concepts) if key in c.lower() or c.lower() in key]
    if len(sub) == 1:
        return sub[0], 'substring'
    raise PromptError(
        'no such concept in the corpus vocabulary', 404,
        {'vocabulary_size': len(vb.concepts),
         'vocabulary_is': 'the 200 concept names of the THINGS-EEG2 designated '
                          'test set.  this endpoint retrieves from the corpus and '
                          'cannot answer for a word the corpus does not contain.',
         'candidates': sorted(vb.concepts[i] for i in sub)[:20],
         'example_concepts': vb.concepts[:12]})


def _lookup_word(ab, text):
    key = ' '.join(str(text).lower().split())
    if not key:
        raise PromptError('text payload is empty')
    hits = ab.vocab.get(key)
    if not hits:
        raise PromptError(
            'no such word in the auditory corpus', 404,
            {'vocabulary_size': len(ab.vocab),
             'vocabulary_is': 'the %d distinct words the forced aligner places in '
                              '"A Study in Scarlet" as read to the LibriBrain '
                              'participant.  this endpoint retrieves from the '
                              'corpus and cannot answer for a word nobody said.'
                              % len(ab.vocab),
             'example_words': sorted(ab.vocab)[:20]})
    return hits


def _row_for_word(ab, run_i, t):
    """the bank-addressable row whose stimulus window ENDS just after the word."""
    import numpy as np
    stem = ab.run_stems[run_i]
    r0 = ab.index_meta['runs'][run_i]['row0']
    n = ab.index_meta['runs'][run_i]['n_rows']
    tt = ab.row_t[r0:r0 + n]
    j = int(np.searchsorted(tt, t))
    q = r0 + min(max(j + 25, ab.ctx), n - ab.win - 1)   # 100 ms of following context
    return int(min(max(q, ab.ctx), ab.n - ab.win - 1))


def _image_from_payload(vb, payload):
    import numpy as np
    if 'corpus_index' in payload:
        i = payload['corpus_index']
        if not isinstance(i, int) or isinstance(i, bool) or not 0 <= i < len(vb.images):
            raise PromptError('corpus_index must be an integer in [0, %d)'
                              % len(vb.images))
        return vb.images[i], {'kind': 'corpus_image', 'corpus_index': i,
                              'concept': vb.concepts[i],
                              'image_path': str(vb.paths[i])}
    if 'image_base64' in payload:
        from PIL import Image
        raw = _decode_b64(payload['image_base64'], 'image_base64')
        h, w = vb.image_shape()[:2]
        try:
            im = Image.open(io.BytesIO(raw)).convert('RGB').resize((w, h),
                                                                   Image.BILINEAR)
        except Exception:
            raise PromptError('image_base64 did not decode as an image')
        return np.asarray(im, dtype='uint8'), {
            'kind': 'uploaded_image', 'bytes': len(raw),
            'resized_to': [h, w],
            'caveat': 'the head was fitted on THINGS-EEG2 stimuli at %dx%d.  an '
                      'uploaded image is off that distribution and the 63.50%% '
                      'figure does not transfer to it -- no accuracy has been '
                      'measured for arbitrary images.' % (h, w)}
    raise PromptError('image payload needs corpus_index or image_base64')


def _audio_from_payload(ab, payload):
    if 'corpus_index' in payload:
        q = payload['corpus_index']
        lo, hi = ab.ctx, ab.n - ab.win - 1
        if not isinstance(q, int) or isinstance(q, bool) or not lo <= q < hi:
            raise PromptError('corpus_index must be an integer in [%d, %d)' % (lo, hi))
        held = q >= ab.ntr + ab.GAP
        return ab.cochleagram_of(q), {
            'kind': 'corpus_cochleagram', 'corpus_index': int(q),
            'entry': ab.entry(q),
            'split': 'held-out' if held else 'TRAINING SPLIT -- this window was '
                                             'seen during fitting; a retrieval from '
                                             'it is not a held-out result'}
    if 'audio_base64' in payload:
        raw = _decode_b64(payload['audio_base64'], 'audio_base64')
        coch, how = ab.cochleagram_from_wav(raw)
        return coch, {'kind': 'uploaded_audio', 'bytes': len(raw), 'cochleagram': how,
                      'caveat': 'the head was fitted on one participant listening to '
                                'one audiobook.  an uploaded clip is off that '
                                'distribution and the corpus figure does not transfer '
                                'to it.'}
    raise PromptError('audio payload needs corpus_index or audio_base64')


def prompt(service, data):
    """{modality, want, payload} -> the response body.  raises PromptError."""
    if not isinstance(data, dict):
        raise PromptError('request body must be a JSON object')
    modality = data.get('modality')
    want = data.get('want')
    payload = data.get('payload', {})
    if not isinstance(payload, dict):
        raise PromptError('payload must be a JSON object')
    if modality not in ('image', 'audio', 'text') or want not in ('image', 'audio', 'text'):
        raise PromptError("modality and want must each be 'image', 'audio' or 'text'")
    pair = (modality, want)
    if pair in UNSUPPORTED:
        raise PromptError(
            '%s -> %s is not one of the six declared paths' % pair, 501,
            {'why': UNSUPPORTED[pair],
             'declared_paths': sorted('%s->%s' % p for p in PATHS),
             'refused_rather_than': 'composing the two corpora through a shared '
                                    'label, which would present an association '
                                    'nothing here has measured'})
    if pair not in PATHS:
        raise PromptError('%s -> %s is not a declared path' % pair, 501,
                          {'declared_paths': sorted('%s->%s' % p for p in PATHS)})

    resynth_k = data.get('resynthesise_top_k', 1)
    if not isinstance(resynth_k, int) or isinstance(resynth_k, bool) or not 0 <= resynth_k <= TOP_K:
        raise PromptError('resynthesise_top_k must be an integer in [0, %d]' % TOP_K)

    vis_report, aud_report = service.measured()
    sim_caveat = _similarity_caveat(vis_report)
    label = '%s->%s' % pair

    # ---------------------------------------------------------------- visual in
    if modality == 'image':
        vb = service.visual()
        img, described = _image_from_payload(vb, payload)
        z = vb.embed_image(img)
        order, sims, top1, margin = vb.retrieve(z)
        results = _visual_results(vb, order, sims, want == 'image')
        extra = {
            'input': described,
            'top1_cosine': round(top1, 6),
            'margin_over_rank2': round(margin, 6),
            'how': 'the image drives the occipital port; the E/I dynamics run for '
                   '%d integrator substeps (%.1f ms) and the state at the trained '
                   'readout is embedded and matched against the bank.'
                   % (vb.readout, 1000 * vb.h * vb.readout),
            'caveats': list(vis_report['caveats']),
        }
        if want == 'text':
            extra['caveats'] = [
                'the text is the concept name of the matched bank entry, parsed '
                'from the dataset\'s own path.  the designated test set holds one '
                'image per concept, so label accuracy and image-retrieval accuracy '
                'are THE SAME NUMBER by construction -- this is a relabelling of '
                'the image ranking, not a second decode.'] + extra['caveats']
        return _envelope(label, want, results, vb.chance, len(vb.images),
                         _visual_measured(vis_report, label), sim_caveat, extra,
                         target_in_bank=described['kind'] == 'corpus_image')

    # ----------------------------------------------------------------- audio in
    if modality == 'audio':
        ab = service.audio()
        coch, described = _audio_from_payload(ab, payload)
        pool_rows, pool_bank, in_bank = ab.pool_for(described.get('corpus_index'))
        described['in_bank'] = in_bank
        z = ab.embed(coch)
        rows, sims, top1, margin = ab.retrieve(z, pool_rows, pool_bank)
        results = _audio_results(ab, rows, sims, resynth_k if want == 'audio' else 0)
        extra = {
            'input': described,
            'top1_cosine': round(top1, 6),
            'margin_over_rank2': round(margin, 6),
            'how': 'the cochleagram window drives the temporal port (a different '
                   'slice of the sheet from the occipital one); the dynamics run '
                   'for %d substeps and the state is embedded and matched against '
                   'the measured-MEG bank.'
                   % (ab.cfg['n_steps'] * ab.cfg['substeps']),
            'caveats': [
                'THIS PATH IS WEAK.  %.2f%% +/- %.2f top-1 over 8 held-out pools of '
                '%d, against a measured dynamics-free control at %.2f%% +/- %.2f.  '
                'the auditory head does not beat its own control; the visual paths '
                'reach %.2f%% at the same chance level.'
                % (100 * aud_report['top1_mean_8_pools'],
                   100 * aud_report['top1_sd_8_pools'], aud_report['pool'],
                   100 * aud_report['dynamics_free_control']['top1'],
                   100 * aud_report['dynamics_free_control']['sd'],
                   100 * vis_report['ablations']['full']['top1']),
                'the corpus is ONE participant listening to ONE audiobook.  nothing '
                'here generalises across people.',
                'v1 and v2 of this corpus were misaligned by up to +/-3.4 s and '
                'every auditory result measured on them is void; this is v3.',
            ],
        }
        if want == 'text':
            extra['caveats'].insert(1,
                'the words are a FORCED ALIGNER\'S output over the retrieved '
                'window, not a decode of speech.  the retrieval is what was '
                'measured; the words are what the corpus says was being said at '
                'the moment the retrieved MEG was recorded.')
        return _envelope(label, want, results, ab.chance, ab.pool,
                         _audio_measured(aud_report, label,
                                         vis_report['ablations']['full']['top1']),
                         sim_caveat, extra,
                         target_in_bank=described.get('in_bank', False))

    # ------------------------------------------------------------------ text in
    text = payload.get('text')
    if not isinstance(text, str) or len(text) > 200:
        raise PromptError('text payload must be a string of at most 200 characters')

    if want == 'image':
        vb = service.visual()
        i, how = _lookup_concept(vb, text)
        z = vb.embed_image(vb.images[i])
        order, sims, top1, margin = vb.retrieve(z)
        results = _visual_results(vb, order, sims, True)
        measured = _visual_measured(vis_report, label)
        measured['vocabulary_lookup'] = {
            'accuracy': 'exact -- a string match over the corpus vocabulary, no '
                        'model involved and nothing to be wrong about',
            'vocabulary_size': len(vb.concepts)}
        return _envelope(label, 'image', results, vb.chance, len(vb.images),
                         measured, sim_caveat, target_in_bank=True, extra={
            'input': {'kind': 'text', 'text': text},
            'vocabulary_lookup': {
                'matched_concept': vb.concepts[i], 'match': how,
                'corpus_index': int(i), 'image_path': str(vb.paths[i]),
                'image_png_base64': _png_base64(vb.images[i]),
                'this_is': 'the corpus entry the word names, returned exactly.  no '
                           'model ran to produce it.'},
            'top1_cosine': round(top1, 6),
            'margin_over_rank2': round(margin, 6),
            'how': 'two steps, reported separately because they have different '
                   'error rates: an EXACT vocabulary lookup to the corpus entry, '
                   'then that entry driven through the cortex and decoded, which is '
                   'the %.2f%% retrieval.  results[] is the retrieval; '
                   'vocabulary_lookup is the exact part.'
                   % (100 * vis_report['ablations']['full']['top1']),
            'caveats': ['results[] is what the cortex retrieves when driven with '
                        'the looked-up image.  it can miss -- and does, %.1f%% of '
                        'the time -- even though the lookup itself cannot.'
                        % (100 * (1 - vis_report['ablations']['full']['top1']))]
                       + list(vis_report['caveats']),
        })

    # text -> audio
    ab = service.audio()
    hits = _lookup_word(ab, text)
    run_i, t, dur = hits[0]
    q = _row_for_word(ab, run_i, t)
    coch = ab.cochleagram_of(q)
    pool_rows, pool_bank, _ = ab.pool_for(q)
    z = ab.embed(coch)
    rows, sims, top1, margin = ab.retrieve(z, pool_rows, pool_bank)
    results = _audio_results(ab, rows, sims, resynth_k)
    lookup_audio, lookup_how = ab.resynthesise(coch, ab.entry(q)['chapter_file'])
    measured = _audio_measured(aud_report, label,
                               vis_report['ablations']['full']['top1'])
    measured['vocabulary_lookup'] = {
        'accuracy': 'exact -- a string match over the aligner\'s word list, no '
                    'model involved',
        'vocabulary_size': len(ab.vocab),
        'caveat': 'the word list is a forced aligner\'s output; the source card '
                  'grades the alignments stream as context for that reason'}
    return _envelope(label, 'audio', results, ab.chance, ab.pool, measured,
                     sim_caveat, target_in_bank=True, extra={
        'input': {'kind': 'text', 'text': text},
        'vocabulary_lookup': {
            'matched_word': text.lower(),
            'occurrences_in_corpus': len(hits),
            'used': ab.entry(q),
            'word_onset_chapter_s': round(float(t), 3),
            'word_duration_s': round(float(dur), 3),
            'audio_wav_base64': _wav_base64(lookup_audio, ab.FS),
            'resynthesis': lookup_how,
            'this_is': 'the corpus cochleagram window containing the word, '
                       'resynthesised.  no model ran to select it.'},
        'top1_cosine': round(top1, 6),
        'margin_over_rank2': round(margin, 6),
        'how': 'an EXACT vocabulary lookup to a corpus cochleagram window, then '
               'that window driven through the cortex and matched against the '
               'measured-MEG bank.  vocabulary_lookup is the exact part; results[] '
               'is the %.2f%% retrieval.'
               % (100 * aud_report['top1_mean_8_pools']),
        'caveats': [
            'THIS PATH IS WEAK.  the retrieval in results[] runs at %.2f%% +/- '
            '%.2f top-1 against a dynamics-free control at %.2f%% +/- %.2f.  only '
            'the vocabulary_lookup is exact.'
            % (100 * aud_report['top1_mean_8_pools'],
               100 * aud_report['top1_sd_8_pools'],
               100 * aud_report['dynamics_free_control']['top1'],
               100 * aud_report['dynamics_free_control']['sd']),
            'the returned waveform is a RESYNTHESIS from a 64-band cochleagram '
            'with the phase discarded.  it is a vocoding of a corpus recording, '
            'not the recording and not generated speech.',
            'the corpus is one participant and one audiobook.',
        ],
    })


def describe(root, service=None):
    """GET /api/brain/prompt -- what the six paths are and what each one measures.

    the accuracies come out of the measurement files, never out of this module, so
    a number cannot drift out of step with the run that produced it.
    """
    svc = service or _service(root)
    vis, aud = svc.measured()
    return {
        'schema': 'ihm.brain-prompt-catalog.v1',
        'post': '/api/brain/prompt',
        'body': {'modality': 'image | audio | text', 'want': 'image | audio | text',
                 'payload': {'corpus_index': 'int', 'image_base64': 'PNG/JPEG bytes',
                             'audio_base64': '16-bit PCM WAV bytes', 'text': 'str'},
                 'resynthesise_top_k': 'int 0-5, audio outputs only (default 1)'},
        'output_is_retrieval': True,
        'generated': False,
        'output_is': 'every output of every path is a RETRIEVAL from the corpus.  '
                     'the brain does not draw, write or speak.',
        'similarity_is_not_confidence': _similarity_caveat(vis),
        'paths': [
            {'path': 'image->image', 'how': 'drive the optic nerve, embed the '
             'cortical state, retrieve the nearest MEASURED evoked EEG entry and '
             'return its image', 'top1': vis['ablations']['full']['top1'],
             'chance': vis['chance'], 'n_bank_entries': vis['n']},
            {'path': 'image->text', 'how': 'the same retrieval, read the concept '
             'name off the matched entry', 'top1': vis['ablations']['full']['top1'],
             'chance': vis['chance'], 'n_bank_entries': vis['n'],
             'note': 'identical to image->image by construction: one image per '
                     'concept on this set, so the label ranking is a relabelling'},
            {'path': 'text->image', 'how': 'exact vocabulary lookup, then that '
             'entry driven through the cortex and decoded',
             'lookup': 'exact', 'top1': vis['ablations']['full']['top1'],
             'chance': vis['chance'], 'n_bank_entries': vis['n']},
            {'path': 'text->audio', 'how': 'exact word lookup to a corpus '
             'cochleagram window, resynthesised; then that window driven through '
             'the cortex and matched against the measured-MEG bank',
             'lookup': 'exact', 'top1': aud['top1_mean_8_pools'],
             'top1_sd': aud['top1_sd_8_pools'], 'chance': aud['chance'],
             'n_bank_entries': aud['pool']},
            {'path': 'audio->audio', 'how': 'cochleagram in, temporal port driven, '
             'nearest MEASURED MEG entry retrieved, its paired cochleagram '
             'resynthesised', 'top1': aud['top1_mean_8_pools'],
             'top1_sd': aud['top1_sd_8_pools'], 'chance': aud['chance'],
             'n_bank_entries': aud['pool']},
            {'path': 'audio->text', 'how': 'the same retrieval, read the aligner\'s '
             'words off the matched entry\'s window',
             'top1': aud['top1_mean_8_pools'], 'top1_sd': aud['top1_sd_8_pools'],
             'chance': aud['chance'], 'n_bank_entries': aud['pool']},
        ],
        'refused': [{'path': '%s->%s' % k, 'why': v} for k, v in
                    sorted(UNSUPPORTED.items())],
        'audio_is_much_weaker_than_vision': _audio_measured(
            aud, 'audio->*', vis['ablations']['full']['top1'])['verdict'],
    }


def handle(root, data, service=None):
    """the route body.  returns (payload, status)."""
    svc = service or _service(root)
    try:
        return prompt(svc, data), 200
    except PromptError as e:
        body = {'error': str(e), 'schema': SCHEMA}
        body.update(e.detail)
        return body, e.status


_SERVICES = {}
_SERVICE_LOCK = threading.Lock()


def _service(root):
    key = str(Path(root).resolve())
    with _SERVICE_LOCK:
        if key not in _SERVICES:
            _SERVICES[key] = BrainPrompt(root)
        return _SERVICES[key]


def describe_handle(root, service=None):
    """the GET route body.  returns (payload, status)."""
    try:
        return describe(root, service), 200
    except PromptError as e:
        body = {'error': str(e), 'schema': 'ihm.brain-prompt-catalog.v1'}
        body.update(e.detail)
        return body, e.status
