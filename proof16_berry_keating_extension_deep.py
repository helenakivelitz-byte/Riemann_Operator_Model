#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Proof 16 — Berry–Keating Extension (optimized)

Ziel:
  1) Universality-Fit von N(E) ≈ a E log E + b E + c log E + d (oberes Quartil)
  2) BK-Vergleich: N_BK(E) = (E log E + k E) / (2π), skaliert mit λ := a_target / a
  3) Schnellere, rigorosere Korrelation (FFT, Phasen-Surrogate)
  4) Optionale DC-Dämpfung in |Trace exp(-iHt)|

Outputs:
  - proof16_universal_scaling_L{L}.png
  - proof16_counting_L{L}.png
  - proof16_trace_overlay_L{L}.png
  - proof16_rigorous_correlation_L{L}.png
  - proof16_summary.csv
"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# SciPy ist nur für curve_fit nötig. Fallback: linearer Ersatz.
try:
    from scipy.optimize import curve_fit
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


# -------------------------- Utilities --------------------------

def load_eigenvalues(path: Path) -> np.ndarray:
    eigs = np.loadtxt(path, dtype=float).ravel()
    eigs = eigs[np.isfinite(eigs)]
    eigs = eigs[eigs > 0]
    eigs.sort()
    return eigs


def empirical_counting(eigs: np.ndarray,
                       E_grid: np.ndarray) -> np.ndarray:
    """N_emp(E) = # { E_n ≤ E } über sortiertes Eigenwert-Array."""
    idx = np.searchsorted(eigs, E_grid, side="right")
    return idx.astype(float)


def density_histogram(eigs: np.ndarray,
                      bins: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Grobes Dichte-Histogramm ρ(E)."""
    hist, edges = np.histogram(eigs, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    width = np.diff(edges)
    rho = hist / (width + 1e-12)
    return centers, rho


def trace_quantum(eigs: np.ndarray, tmax: float, dt: float,
                  remove_dc: bool) -> Tuple[np.ndarray, np.ndarray]:
    """Spur e^{-iHt} = ∑_n e^{-i E_n t} und |·|."""
    t = np.arange(0.0, tmax, dt, dtype=float)
    tr = np.zeros_like(t, dtype=complex)
    for E in eigs:
        tr += np.exp(-1j * E * t)
    q_abs = np.abs(tr).astype(float)
    if remove_dc:
        k = min(100, len(q_abs)-1)
        base = np.median(q_abs[1:k]) if k > 1 else q_abs[1] if len(q_abs) > 1 else q_abs[0]
        q_abs[:3] = base
    return t, q_abs


def trace_bk_proxy(t: np.ndarray) -> np.ndarray:
    """
    Sehr einfache, glatte BK-Proxy-Amplitude (keine Prime-Details):
    ~ (t + t0)^(-1) * exp(-t/t1). Reicht für Konsistenzplots.
    """
    t0, t1 = 0.02, 0.15
    y = np.exp(-t/t1) / (t + t0)
    y = (y - y.min()) / (y.max() - y.min() + 1e-12)
    return y


# -------------------- Universality Fit (Oberes Quartil) --------------------

def _model_N(E, a, b, c, d):
    E = np.asarray(E, float)
    return a*E*np.log(np.maximum(E,1.0)) + b*E + c*np.log(np.maximum(E,1.0)) + d

def fit_universality(E_sorted: np.ndarray, N_emp: np.ndarray,
                     frac: float = 0.75) -> Dict[str, float]:
    """
    Fit N(E) ≈ a E log E + b E + c log E + d
    nur auf oberem Quantil (default 75%..100%) für asymptotische Stabilität.
    """
    n = len(E_sorted)
    i0 = int(frac * n)
    E_fit = E_sorted[i0:]
    N_fit = N_emp[i0:]
    # Startwerte: grob
    p0 = [1e-3, 0.0, 0.0, 0.0]
    bounds = ([-np.inf, -np.inf, -np.inf, -np.inf],
              [ np.inf,  np.inf,  np.inf,  np.inf])
    if _HAVE_SCIPY:
        popt, _ = curve_fit(_model_N, E_fit, N_fit, p0=p0, bounds=bounds, maxfev=20000)
    else:
        # Fallback: lineare Regression auf Features
        X = np.vstack([E_fit*np.log(np.maximum(E_fit,1.0)),
                       E_fit,
                       np.log(np.maximum(E_fit,1.0)),
                       np.ones_like(E_fit)]).T
        popt, *_ = np.linalg.lstsq(X, N_fit, rcond=None)
    a, b, c, d = map(float, popt)
    return {"a": a, "b": b, "c": c, "d": d, "i0": i0}


def bk_k_from_linear_reg(E: np.ndarray, N_emp: np.ndarray) -> float:
    """
    Minimiert ∑ (N_emp - (E log E + k E)/(2π))^2  →  geschlossene Lösung für k.
    """
    two_pi = 2.0*np.pi
    y = two_pi*N_emp - E*np.log(np.maximum(E,1.0))
    denom = np.sum(E*E) + 1e-12
    k = np.sum(E*y)/denom
    return float(k)


# -------------------- FFT-Korrelation + Surrogates --------------------

def rigorous_correlation_analysis(quantum_trace: np.ndarray,
                                  bk_trace: np.ndarray,
                                  L: float,
                                  n_surrogates: int = 20,
                                  max_lag: int | None = None,
                                  seed: int = 123) -> Tuple[float, int, float]:
    rng = np.random.default_rng(seed)

    def zscore(x):
        x = np.asarray(x, float)
        x = x - np.nanmean(x)
        return x / (np.nanstd(x) + 1e-12)

    x = zscore(quantum_trace); y = zscore(bk_trace)
    x -= x.mean(); y -= y.mean()

    N = len(x)
    if max_lag is None:
        max_lag = min(300, N//4)

    # FFT-basiert:
    fx = np.fft.rfft(x); fy = np.fft.rfft(y)
    corr_full = np.fft.irfft(fx * np.conj(fy), n=N)
    corr_full /= (np.std(x)*np.std(y)*N + 1e-12)
    corr_full = np.roll(corr_full, N//2)
    lags = np.arange(-N//2, N//2)
    m = (lags >= -max_lag) & (lags <= max_lag)
    lags = lags[m]; correlations = corr_full[m]

    # Surrogates (Phasen randomisieren)
    null_vals = []
    for _ in range(n_surrogates):
        Y = np.fft.rfft(y)
        phases = rng.uniform(0, 2*np.pi, len(Y))
        Ys = np.abs(Y) * np.exp(1j*phases)
        ys = np.fft.irfft(Ys, n=N)
        ys = zscore(ys); ys -= ys.mean()

        fys = np.fft.rfft(ys)
        cs = np.fft.irfft(fx * np.conj(fys), n=N)
        cs /= (np.std(x)*np.std(ys)*N + 1e-12)
        cs = np.roll(cs, N//2)
        null_vals.append(np.max(np.abs(cs[m])))

    null_vals = np.asarray(null_vals)
    thresh = float(np.percentile(null_vals, 95))
    i = int(np.argmax(np.abs(correlations)))
    corr_max = float(np.abs(correlations[i])); lag_max = int(lags[i])

    # Plot
    plt.figure(figsize=(10, 4))
    plt.plot(lags, correlations, 'k-', lw=1.5, label='Pearson xcorr')
    plt.axhline(+thresh, ls='--', c='r', label=f'null max ≈ {thresh:.2f}')
    plt.axhline(-thresh, ls='--', c='r')
    plt.title(f"Proof 16 — Quantum vs BK correlation (L={L})\n"
              f"max={corr_max:.3f} at lag={lag_max}")
    plt.xlabel("Lag (samples)"); plt.ylabel("Correlation (normalized)")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f'proof16_rigorous_correlation_L{L}.png', dpi=150); plt.close()

    print(f"[corr] max={corr_max:.3f} @ lag={lag_max}  (null@95%≈{thresh:.2f}, surrogates={n_surrogates})")
    return corr_max, lag_max, thresh


# ------------------------------ Plots ------------------------------

def plot_universality(E_scaled: np.ndarray, N_emp: np.ndarray,
                      a: float, b: float, c: float, d: float,
                      L: float, lam: float, k_hat: float) -> None:
    two_pi = 2.0*np.pi
    # BK auf skaliertem E_phys':=E/λ (gleichwertig: N_BK(E_scaled) direkt)
    N_bk = (E_scaled*np.log(np.maximum(E_scaled,1.0)) + k_hat*E_scaled) / two_pi

    # Residuen
    R = N_emp - N_bk

    # Dichte (empirisch & BK)
    centers, rho_emp = density_histogram(E_scaled, bins=80)
    rho_bk = (np.log(np.maximum(centers,1.0)) + 1.0 + 0.0*centers) / two_pi  # derivative ~ (log E + 1)/(2π)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    ax = axes[0]
    ax.plot(E_scaled, N_emp, 'b.', ms=1.5, alpha=0.6, label="Empirisch N(E)")
    ax.plot(E_scaled, N_bk, 'r-', lw=1.4, label=f"BK-Fit: k̂={k_hat:.3f}")
    ax.set_xlabel("E (reskaliert)"); ax.set_ylabel("N(E)")
    ax.set_title(f"Skalenfaktor: {lam:.3g}")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(E_scaled, R, 'g-', lw=1.4)
    ax.axhline(0, color='r', ls='--', alpha=0.5)
    ax.set_xlabel("E"); ax.set_ylabel("N_emp - N_BK")
    ax.set_title("Residuen"); ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(centers, rho_emp, color='tab:blue', lw=1.4, label="Empirische Dichte")
    ax.plot(centers, rho_bk, color='tab:red', lw=1.4, label="BK-Dichte")
    ax.set_xlabel("E (reskaliert)"); ax.set_ylabel("ρ(E)")
    ax.set_title("Spektrale Dichte"); ax.legend(); ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"proof16_universal_scaling_L{L}.png", dpi=150)
    plt.close()


def plot_counting(E_scaled: np.ndarray, N_emp: np.ndarray,
                  k_hat: float, lam: float, L: float) -> None:
    two_pi = 2.0*np.pi
    N_bk = (E_scaled*np.log(np.maximum(E_scaled,1.0)) + k_hat*E_scaled) / two_pi

    plt.figure(figsize=(9, 6))
    plt.plot(E_scaled, N_emp, '-', lw=2, color='tab:blue', label="N_emp(E) (operator)")
    plt.plot(E_scaled, N_bk, '--', lw=2, color='tab:orange',
             label=f"N_BK(E'; k={k_hat:.3f}), λ={lam:g}")
    plt.xlabel("E"); plt.ylabel("Counting function N(E)")
    plt.title(f"Proof 16 — Counting comparison (L={L})")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"proof16_counting_L{L}.png", dpi=150); plt.close()


def plot_trace_overlay(t: np.ndarray, q_abs: np.ndarray, L: float) -> None:
    y_bk = trace_bk_proxy(t)
    q = (q_abs - q_abs.min())/(q_abs.max()-q_abs.min()+1e-12)

    plt.figure(figsize=(10, 6))
    plt.plot(t, q, color='tab:blue', lw=1.1, label='|Trace e^{-iHt}| (Quantum)')
    plt.plot(t, y_bk, color='tab:red', lw=1.5, label='|Trace_BK| (Semiclassical, proxy)')
    for p in [2,3,5,7,11,13,17,19,23,29]:
        plt.axvline(np.log(p), color='gray', ls='--', alpha=0.2)
    plt.xlabel("t"); plt.ylabel("Normalized amplitude")
    plt.title(f"Proof 16 — Trace overlay (L={L})")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"proof16_trace_overlay_L{L}.png", dpi=150); plt.close()


# ------------------------------ Main pipeline ------------------------------

def run_pipeline(eigs: np.ndarray, L: float, tmax: float, dt: float,
                 remove_dc: bool, n_surrogates: int, max_lag: int | None) -> Dict[str, float]:
    # Zähl-Funktion auf dichtem E-Gitter
    E_grid = np.linspace(eigs[0], eigs[-1], 1600)
    N_emp = empirical_counting(eigs, E_grid)

    # Universality-Fit (oberes Quartil)
    fit = fit_universality(E_grid, N_emp, frac=0.75)
    a, b, c, d = fit["a"], fit["b"], fit["c"], fit["d"]

    a_target = 1.0/(2.0*np.pi)  # theoretisch
    lam = a_target / a if a != 0.0 else 1.0  # E' = λ E

    # E skaliert + neuer BK-k (lineare Regression)
    E_scaled = lam * E_grid
    k_hat = bk_k_from_linear_reg(E_scaled, N_emp)

    # Plots (Counting, Universality)
    plot_universality(E_scaled, N_emp, a, b, c, d, L, lam, k_hat)
    plot_counting(E_scaled, N_emp, k_hat, lam, L)

    # Traces + Korrelation
    t, q_abs = trace_quantum(eigs, tmax, dt, remove_dc=remove_dc)
    plot_trace_overlay(t, q_abs, L)

    # Für Korrelation BK-Proxy auf *gleichem* t-Gitter
    y_bk = trace_bk_proxy(t)
    corr_max, lag_max, thresh = rigorous_correlation_analysis(
        q_abs, y_bk, L, n_surrogates=n_surrogates, max_lag=max_lag
    )

    return {
        "a": a, "b": b, "c": c, "d": d,
        "a_target": a_target, "lambda": lam, "k_hat": k_hat,
        "corr_max": corr_max, "lag_max": lag_max, "null95": thresh,
        "tmax": tmax, "dt": dt, "remove_dc": bool(remove_dc),
    }


def main():
    ap = argparse.ArgumentParser(description="Proof 16 — Berry–Keating Extension (optimized)")
    ap.add_argument("--eigenfile", type=str, default="eigenvalues_clean.txt")
    ap.add_argument("--L", type=float, default=100.0)
    ap.add_argument("--tmax", type=float, default=5.0)
    ap.add_argument("--dt", type=float, default=0.002)
    ap.add_argument("--remove_dc", action="store_true",
                    help="Dämpft den DC-Spike der Quantenspur")
    ap.add_argument("--n_surrogates", type=int, default=20,
                    help="Anzahl Surrogate für Null-Hypothese (95%-Schwelle)")
    ap.add_argument("--max_lag", type=int, default=None,
                    help="Maximaler Lag für Korrelation (Default=min(300,N/4))")
    args = ap.parse_args()

    print("=== Proof 16: Berry–Keating Extension (optimized) ===")
    eig_path = Path(args.eigenfile)
    if not eig_path.exists():
        raise FileNotFoundError(f"Eigenvalue file not found: {eig_path}")

    eigs = load_eigenvalues(eig_path)
    print(f"[proof16] Loaded {len(eigs)} eigenvalues from {eig_path.name}")
    print(f"[proof16] E-min={eigs[0]:.4g}, E-max={eigs[-1]:.4g}")

    res = run_pipeline(eigs, args.L, args.tmax, args.dt,
                       args.remove_dc, args.n_surrogates, args.max_lag)

    # Summary CSV
    df = pd.DataFrame([res])
    df.to_csv("proof16_summary.csv", index=False)
    print("[proof16] Summary → proof16_summary.csv")

    # Console key numbers
    print("\n=== Key numbers (optimized) ===")
    print(f"• a = {res['a']:.6e}  (target 1/(2π) = {res['a_target']:.6e})")
    print(f"• λ (scale from a) = {res['lambda']:.6g}")
    print(f"• k̂ (BK on E') = {res['k_hat']:.3f}")
    print(f"• corr_max = {res['corr_max']:.3f}  (null@95%≈{res['null95']:.2f}) @ lag={res['lag_max']}")
    print(f"• remove_dc={args.remove_dc}, L={args.L}, tmax={args.tmax}, dt={args.dt}")
    if not _HAVE_SCIPY:
        print("! SciPy nicht gefunden → Fallback-LSQ genutzt (ok, aber Fit kann etwas weniger stabil sein).")


if __name__ == "__main__":
    main()
