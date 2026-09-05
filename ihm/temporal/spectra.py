"""Frequency axes in Hz; Laplace s and continuous state A in inverse seconds."""
import numpy as np
from scipy import signal
from scipy.integrate import trapezoid


def _signal(x, fs):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError('signal must be a finite vector with at least two samples')
    if not np.isfinite(fs) or fs <= 0:
        raise ValueError('sample rate must be positive and finite')
    return x


def _settings(x, fs, nperseg):
    x = _signal(x, fs)
    n = min(len(x), 1024) if nperseg is None else nperseg
    if not isinstance(n, (int, np.integer)) or not 2 <= n <= len(x):
        raise ValueError('nperseg must be an integer between 2 and signal length')
    return x, dict(fs=fs, window='hann', nperseg=n, noverlap=0, detrend='constant', scaling='density')


def spectral_estimate(x, fs, *, nperseg=None):
    """Non-overlapping Hann Welch density in signal-unit squared per Hz.

    Remainder samples are excluded. Window counts are not claimed independent
    for autocorrelated physiology; no nominal confidence intervals are invented.
    """
    x, kw = _settings(x, fs, nperseg)
    f, p = signal.welch(x, **kw)
    return dict(frequency_hz=f.tolist(), psd=p.tolist(), nperseg=kw['nperseg'],
                segments=len(x)//kw['nperseg'], dropped_samples=len(x)%kw['nperseg'],
                frequency_resolution_hz=float(fs/kw['nperseg']), window='periodic Hann',
                detrend='segment mean', overlap_samples=0)


def cross_spectrum(x, y, fs, *, nperseg=None):
    """CSD convention conj(X)*Y; undefined zero-power coherence is JSON null."""
    x, kw = _settings(x, fs, nperseg)
    y = _signal(y, fs)
    if x.shape != y.shape:
        raise ValueError('signals must share sample count and clock')
    f, pxy = signal.csd(x, y, **kw)
    _, pxx = signal.welch(x, **kw)
    _, pyy = signal.welch(y, **kw)
    denominator = pxx*pyy
    c = np.divide(abs(pxy)**2, denominator, out=np.full_like(pxx,np.nan), where=denominator>0)
    # One segment gives tautological coherence, not an estimate of coupling.
    if len(x)//kw['nperseg'] < 2:
        c[:] = np.nan
    return dict(frequency_hz=f.tolist(), real=pxy.real.tolist(), imag=pxy.imag.tolist(),
                coherence=[float(np.clip(v,0,1)) if np.isfinite(v) else None for v in c],
                convention='conj(X)*Y', segments=len(x)//kw['nperseg'])


def finite_laplace(time_s, x, s):
    """Trapezoidal integral from first to last sample of x(t)*exp(-s*t).

    Time is elapsed from the first sample, with no implicit detrending or tail.
    Finite-window transforms are entire functions, not estimates of system poles.
    """
    t = np.asarray(time_s,dtype=float)
    x = np.asarray(x)
    s = np.atleast_1d(np.asarray(s,dtype=complex))
    if t.ndim!=1 or len(t)<2 or x.shape!=t.shape or not np.isfinite(t).all() or not np.isfinite(x).all():
        raise ValueError('time and signal must be matching finite vectors')
    if not np.all(np.diff(t)>0) or s.ndim!=1 or not np.isfinite(s).all():
        raise ValueError('time must increase; Laplace points must be a finite vector')
    t=t-t[0]
    if np.max(-s.real*t[-1])>650:
        raise ValueError('Laplace exponential exceeds safe numerical range')
    return np.array([trapezoid(x*np.exp(-z*t),t) for z in s])


def resolvent(A, B, C, s, D=None):
    """Continuous-time LTI transfer C(sI-A)^-1 B+D; poles raise LinAlgError."""
    A,B,C=(np.asarray(a,dtype=complex) for a in (A,B,C))
    if A.ndim!=2 or A.shape[0]!=A.shape[1] or B.ndim!=2 or C.ndim!=2 or B.shape[0]!=len(A) or C.shape[1]!=len(A):
        raise ValueError('incompatible state-space dimensions')
    D=np.zeros((len(C),B.shape[1])) if D is None else np.asarray(D,dtype=complex)
    s=np.atleast_1d(np.asarray(s,dtype=complex))
    if D.shape!=(len(C),B.shape[1]) or s.ndim!=1 or not all(np.isfinite(a).all() for a in (A,B,C,D,s)):
        raise ValueError('invalid feedthrough or nonfinite state-space model')
    return np.array([C@np.linalg.solve(z*np.eye(len(A))-A,B)+D for z in s])
