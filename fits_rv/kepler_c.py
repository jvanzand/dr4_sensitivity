"""Loader for kepler.c: compiles it on first use (cc -O3 -shared -fPIC) into
fits_rv/_build/ and exposes `keplerian_rv_add` and `kepler_E`. If no C compiler is
available, `AVAILABLE` is False and rv_model.py falls back to its numpy solver
(about 2-3x slower, same results)."""

import ctypes
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'kepler.c')
BUILD = os.path.join(HERE, '_build')
LIB = os.path.join(BUILD, 'libkepler' + ('.dylib' if sys.platform == 'darwin' else '.so'))
AVAILABLE = False
_lib = None


def _compile():
    os.makedirs(BUILD, exist_ok=True)
    for cc in (os.environ.get('CC'), 'cc', 'clang', 'gcc'):
        if cc and shutil.which(cc):
            r = subprocess.run([cc, '-O3', '-shared', '-fPIC', '-o', LIB, SRC, '-lm'],
                               capture_output=True, text=True)
            if r.returncode == 0:
                return True
    return False


def _load():
    global _lib, AVAILABLE
    stale = not os.path.exists(LIB) or os.path.getmtime(LIB) < os.path.getmtime(SRC)
    if stale and not _compile():
        return
    try:
        lib = ctypes.CDLL(LIB)
    except OSError:
        return
    dp = np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS')
    lib.keplerian_rv.argtypes = [dp, ctypes.c_int] + [ctypes.c_double]*6 + [dp]
    lib.keplerian_rv.restype = None
    lib.kepler_E.argtypes = [dp, ctypes.c_int, ctypes.c_double, dp]
    lib.kepler_E.restype = None
    _lib = lib
    AVAILABLE = True


_load()


def keplerian_rv_add(out, t, tref, P, K, e, w, M0):
    """out += K [cos(nu + w) + e cos w] for the orbit (P, K, e, w, M0 at tref). `t` and
    `out` must be C-contiguous float64 arrays."""
    _lib.keplerian_rv(t, len(t), tref, P, K, e, w, M0, out)
    return out


def kepler_E(M, e):
    M = np.ascontiguousarray(M, dtype=np.float64)
    E = np.empty_like(M)
    _lib.kepler_E(M, len(M), float(e), E)
    return E
