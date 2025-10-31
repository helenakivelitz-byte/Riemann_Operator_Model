# src/diagnostics_utils.py
# -----------------------------------------------------------
# Utility functions for post-processing diagnostics:
#   - spectral zeta analysis
#   - bridge ratio computation
#   - Weyl law comparison
#   - Berry–Keating rescaling
# -----------------------------------------------------------

import json, os, math
import numpy as np
import mpmath as mp

# global numerical precision for mpmath
mp.mp.dps = 50

# -----------------------------------------------------------
# File + data handling
# -----------------------------------------------------------

def load_results():
    """Load raw and rescaled pipeline result JSONs."""
    with open("data/results.json", "r", encoding="utf-8") as f:
        RAW = json.load(f)
    with open("data/results_rescaled.json", "r", encoding="utf-8") as f:
        RES = json.load(f)
    return RAW, RES


def eigenvalues_from_results(RAW):
    """
    Extract eigenvalue array from result dictionary or list.
    Supports keys: eigenvalues, evals, E, E_raw.
    """
    if isinstance(RAW, dict):
        for key in ("eigenvalues", "evals", "E", "E_raw"):
            if key in RAW:
                return np.array(RAW[key], dtype=float)
    if isinstance(RAW, list):
        return np.array(RAW, dtype=float)
    return np.array([], dtype=float)


def ensure_dirs():
    """Ensure output directories exist."""
    os.makedirs("data/diagnostics", exist_ok=True)
    os.makedirs("data/figures", exist_ok=True)


# -----------------------------------------------------------
# Spectral analysis helpers
# -----------------------------------------------------------

def dense_s_grid(RAW, start=1.05, stop=1.50, step=0.01):
    """
    Build a dense s-grid for bridge ratio scans.
    If pipeline already provided zeta_s_values, reuse them.
    """
    if isinstance(RAW, dict) and "zeta_s_values" in RAW and RAW["zeta_s_values"]:
        svals = [float(s) for s in RAW["zeta_s_values"]]
        if len(svals) >= 10:  # sufficiently dense
            return np.array(sorted(svals))
    return np.arange(start, stop + 1e-12, step)


def spectral_zeta_truncated(E, s, K=None):
    """Compute ζ_H(s) = Σ_{n≤K} E_n^{-s} with truncation (defaults to full length)."""
    if K is None or K > len(E):
        K = len(E)
    s_mp = mp.mpf(s)
    total = mp.mpf('0')
    for En in E[:K]:
        if En > 0:
            total += mp.power(En, -s_mp)
    return total


def bridge_ratio_absdev(E, s, K=None, k=None):
    """
    Compute |Ψ(s) - 1| = |ζ_H(s)/ζ(s) - 1|.
    Accepts both K and k for backward compatibility.
    """
    K_use = K if K is not None else k
    if K_use is None:
        K_use = len(E)

    zH = spectral_zeta_truncated(E, s, K=K_use)
    zR = mp.zeta(s)
    if abs(zR) < mp.mpf("1e-30"):
        return float("inf")
    return float(abs(zH / zR - 1.0))


def weyl_count_leading(L, Emax):
    """1D Dirichlet Weyl-law approximation: N(E) ~ (L/π) √E."""
    if Emax <= 0:
        return 0.0
    return float((L / math.pi) * math.sqrt(Emax))


# -----------------------------------------------------------
# Berry–Keating rescaling utilities
# -----------------------------------------------------------

def quantile_rescale_BK(E, lam):
    """
    Perform monotone Berry–Keating rescaling:
    E'_n = N_BK^{-1}(N_emp(E_n)),  N_BK(E') = λ E' (log E' - 1).
    Inverse solved via mp.findroot with robust fallback.
    """
    E_sorted = np.sort(np.asarray(E, dtype=float))
    n = np.arange(1, len(E_sorted) + 1, dtype=float)

    Eprime = []
    for Nk in n:
        guess = max(10.0, Nk / lam)  # robust initial guess
        f = lambda x: lam * x * (mp.log(x) - 1) - Nk
        try:
            sol = mp.findroot(f, guess)
            Eprime.append(float(sol))
        except Exception:
            # Fallback: asymptotic approximation exp(1 + N/λ)
            Eprime.append(float(mp.e ** (1.0 + Nk / lam)))
    return np.array(Eprime, dtype=float)


# -----------------------------------------------------------
# Simple persistence helpers
# -----------------------------------------------------------

def save_json(obj, path):
    """Save object as formatted JSON."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


# -----------------------------------------------------------
# End of diagnostics_utils.py
# -----------------------------------------------------------
