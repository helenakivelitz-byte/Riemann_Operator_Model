#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Proof 16 — Berry–Keating Extension (finetune, patched)

Features
--------
• robustes Fitten der BK-Skalierung: N(E) ~ (E/(2π λ)) * (log(E/(2π λ)) - c)
• automatische Reskalierung der Eigenwerte (--auto_rescale)
• Trace |sum exp(-i E t)| mit DC-Dämpfung (--remove_dc) und Hann-Fenster
• Kreuzkorrelation (Pearson, normiert) inkl. Null-Band (shuffle-basierte Schwelle)
• saubere Plots + Summary CSV

Start (Windows)
---------------
python proof16_bk_finetune.py --mode eigen --eigenfile eigenvalues_clean.txt ^
  --tmax 5 --dt 0.002 --remove_dc --auto_rescale --plot

"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import curve_fit
    _HAS_SCI = True
except Exception:
    _HAS_SCI = False

# --------------- utilities ----------------

TWOPI = 2.0 * np.pi
EPS = 1e-12


def load_eigs(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Eigenfile not found: {path}")
    eigs = np.loadtxt(path, dtype=float).ravel()
    return eigs


def counting_from_eigs(eigs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (E_sorted_positive, N_empirical)."""
    E = np.sort(eigs[eigs > 0.0])
    N = np.arange(1, len(E) + 1, dtype=float)
    return E, N


def N_BK(E: np.ndarray, lam: float, c: float) -> np.ndarray:
    """
    Berry–Keating counting model with scale λ and offset c:
    N_BK(E; λ, c) = (E / (2π λ)) * [ log(E / (2π λ)) - c ]
    """
    x = np.maximum(E / (TWOPI * max(lam, EPS)), EPS)
    return (E / (TWOPI * max(lam, EPS))) * (np.log(x) - c)


def fit_lambda_c(E: np.ndarray, N: np.ndarray) -> Tuple[float, float]:
    """
    Robust fit for (λ, c). Falls curve_fit fehlt/fehlschlägt,
    nutze kleine Gitter-Suche um die Lösung zu stabilisieren.
    """
    # sinnvolles Fit-Fenster (mittleres bis hohes E für Asymptotik)
    k0 = max(int(0.25 * len(E)), 10)
    Ex, Nx = E[k0:], N[k0:]

    lam0, c0 = 40.0, 1.0  # Heuristik nach bisherigen Resultaten
    bounds = ([1e-3, -20.0], [1e4, 20.0])

    if _HAS_SCI:
        try:
            popt, _ = curve_fit(N_BK, Ex, Nx, p0=[lam0, c0], bounds=bounds, maxfev=20000)
            return float(popt[0]), float(popt[1])
        except Exception:
            pass  # Fallback unten

    # Fallback: grobe Raster-Suche
    lam_grid = np.geomspace(0.5, 200.0, 80)
    c_grid = np.linspace(-10.0, 10.0, 81)
    best = (np.inf, lam0, c0)
    for lam in lam_grid:
        x = np.maximum(Ex / (TWOPI * lam), EPS)
        base = (Ex / (TWOPI * lam)) * np.log(x)
        # Nx ≈ base - (Ex/(2πλ)) * c  → lineare Regression in c
        A = (Ex / (TWOPI * lam))
        c_hat = np.sum((base - Nx) * A) / np.sum(A * A)
        if c_hat < c_grid.min() or c_hat > c_grid.max():
            # clamp, nur um RMSE zu bewerten
            c_hat = float(np.clip(c_hat, c_grid.min(), c_grid.max()))
        resid = Nx - (base - A * c_hat)
        rmse = np.sqrt(np.mean(resid * resid))
        if rmse < best[0]:
            best = (rmse, lam, c_hat)
    return float(best[1]), float(best[2])


def make_trace(eigs: np.ndarray, tmax: float, dt: float, remove_dc: bool) -> Tuple[np.ndarray, np.ndarray]:
    t = np.arange(0.0, tmax + 0.5 * dt, dt)
    tr = np.zeros_like(t, dtype=complex)
    for E in eigs:
        tr += np.exp(-1j * E * t)
    qabs = np.abs(tr).astype(float)
    if remove_dc and len(qabs) >= 6:
        k = min(100, len(qabs) - 1)
        base = np.median(qabs[1:k])
        qabs[0:3] = base
    return t, qabs


def bk_proxy_trace(t: np.ndarray) -> np.ndarray:
    """Eine monotone Proxy-Kurve (glatter Anteil / Dämpfer)."""
    tau = 0.07
    return 1.0 / (t + tau)


def pearson_xcorr(a: np.ndarray, b: np.ndarray, null_draws: int = 200) -> Tuple[np.ndarray, float]:
    """normierte xcorr + simple Null-Schwelle (95%-Quantil aus Phase-Shuffle)."""
    win = np.hanning(len(a))
    ax = (a - np.mean(a)) * win
    bx = (b - np.mean(b)) * win
    corr = np.correlate(ax, bx, mode="full")
    corr = corr / (np.max(np.abs(corr)) + 1e-12)

    # Null: phasen-shuffle von bx (Fourier-Amplitude konstant, Phase random)
    # -> braucht FFT; wenn nicht verfügbar, simple Permutation
    try:
        import numpy.fft as nfft
        B = nfft.rfft(bx)
        mag = np.abs(B)
        rng = np.random.default_rng(1234)
        maxima = []
        for _ in range(null_draws):
            phase = np.exp(1j * rng.uniform(0, 2 * np.pi, size=mag.shape))
            Bnull = mag * phase
            bnull = nfft.irfft(Bnull, n=len(bx))
            cnull = np.correlate(ax, bnull, mode="full")
            cnull = cnull / (np.max(np.abs(cnull)) + 1e-12)
            maxima.append(np.max(np.abs(cnull)))
        thr = float(np.quantile(maxima, 0.95))
    except Exception:
        # Fallback: einfache Zufallspermutation
        rng = np.random.default_rng(1234)
        maxima = []
        for _ in range(null_draws):
            bp = rng.permutation(bx)
            cnull = np.correlate(ax, bp, mode="full")
            cnull = cnull / (np.max(np.abs(cnull)) + 1e-12)
            maxima.append(np.max(np.abs(cnull)))
        thr = float(np.quantile(maxima, 0.95))
    return corr, thr


# --------------- plotting ----------------

def plot_trace_overlay(t: np.ndarray, qabs: np.ndarray, L: float, out: Path):
    proxy = bk_proxy_trace(t)
    plt.figure(figsize=(12, 7))
    qa = qabs / (np.max(qabs) + EPS)
    pb = proxy / (np.max(proxy) + EPS)
    plt.plot(t, qa, color="#1f77b4", lw=1.2, label=r"|Trace $e^{-iHt}$| (Quantum)")
    plt.plot(t, pb, color="#d62728", lw=1.2, label=r"|Trace$_{\rm BK}$| (proxy)")
    plt.xlabel("t")
    plt.ylabel("Normalized amplitude")
    plt.title(f"Proof 16 — Trace overlay (from eigenfile)")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def plot_correlation(corr: np.ndarray, null_thr: float, L: float, out: Path):
    lag = np.arange(-(len(corr) // 2), (len(corr) + 1) // 2)
    vmax = float(np.max(np.abs(corr)))
    arg = int(np.argmax(np.abs(corr)))
    lag_at = lag[arg]
    plt.figure(figsize=(12, 6))
    plt.plot(lag, corr, "k-", lw=1.6, label="Pearson xcorr")
    plt.axhline(+null_thr, ls="--", color="crimson", alpha=0.6, label=f"null max ≈ {null_thr:.2f}")
    plt.axhline(-null_thr, ls="--", color="crimson", alpha=0.6)
    plt.title(f"Proof 16 — Quantum vs BK correlation (from file)\nmax={vmax:.3f} at lag={lag_at}")
    plt.xlabel("Lag (samples)")
    plt.ylabel("Correlation (normalized)")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def plot_counting(E: np.ndarray, N: np.ndarray, lam: float, c: float, L: float, out: Path):
    plt.figure(figsize=(12, 7))
    plt.plot(E, N, ".", ms=2, color="#1f77b4", label="N_emp(E) (from file)")
    Egrid = np.linspace(max(E[5], 1.0), float(E[-1]), 600)
    Nfit = N_BK(Egrid, lam, c)
    plt.plot(Egrid, Nfit, "k--", lw=1.5, label=f"N_BK(E; λ={lam:.3f}, c={c:.2f})")
    plt.xlabel("E")
    plt.ylabel("Counting function N(E)")
    plt.title(f"Proof 16 — Counting (from eigenfile)")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


# --------------- main flow ----------------

def run_from_eigenfile(args: argparse.Namespace) -> Dict:
    eigs_raw = load_eigs(Path(args.eigenfile))
    # nur positive Energien für N(E)
    Epos, Nemp = counting_from_eigs(eigs_raw)

    # Fit λ, c
    lam, c = fit_lambda_c(Epos, Nemp)

    # optionales automatisches Reskalieren (nur für Trace/Korr/Plots sinnvoll)
    lam_corr = 1.0 / lam if args.auto_rescale else 1.0
    eigs = eigs_raw * lam_corr

    # Trace
    t, qabs = make_trace(eigs, args.tmax, args.dt, args.remove_dc)
    corr, null_thr = pearson_xcorr(qabs, bk_proxy_trace(t))

    # Plots
    out_trace = Path(f"proof16_trace_overlay_from_file.png")
    out_corr = Path(f"proof16_correlation_from_file.png")
    out_count = Path(f"proof16_counting_file.png")

    plot_trace_overlay(t, qabs, args.L, out_trace)
    plot_correlation(corr, null_thr, args.L, out_corr)
    plot_counting(Epos, Nemp, lam, c, args.L, out_count)

    vmax = float(np.max(np.abs(corr)))
    lag = np.arange(-(len(corr) // 2), (len(corr) + 1) // 2)
    lag_at = int(lag[int(np.argmax(np.abs(corr)))])

    # Summary
    summary = dict(
        mode="eigen",
        eigenfile=str(args.eigenfile),
        L=float(args.L),
        n_eigs=int(len(eigs_raw)),
        E_min=float(np.min(eigs_raw)),
        E_max=float(np.max(eigs_raw)),
        lambda_hat=float(lam),
        c_hat=float(c),
        auto_rescale=bool(args.auto_rescale),
        corr_max=float(vmax),
        corr_lag=int(lag_at),
        null_thr=float(null_thr),
        tmax=float(args.tmax),
        dt=float(args.dt),
        remove_dc=bool(args.remove_dc),
    )
    pd.DataFrame([summary]).to_csv("proof16_summary.csv", index=False)

    print("=== Proof 16: Berry–Keating Extension (patched) ===")
    print(f"[proof16] Loaded {len(eigs_raw)} eigenvalues from {args.eigenfile}")
    print(f"[proof16] E-min={np.min(eigs_raw):.4g}, E-max={np.max(eigs_raw):.4g}")
    print(f"[proof16] λ̂={lam:.3f}, ĉ={c:.3f}  "
          f"{'(auto-rescaled for trace)' if args.auto_rescale else '(no rescale for trace)'}")
    print(f"[proof16] corr_max={summary['corr_max']:.3f}  (null ≈ {null_thr:.2f})  @ lag={lag_at}")
    print(f"[proof16] Saved plots: {out_trace.name}, {out_corr.name}, {out_count.name}")
    print("[proof16] Summary → proof16_summary.csv")
    return summary


def main():
    p = argparse.ArgumentParser(description="Proof 16 — BK finetune (patched)")
    p.add_argument("--mode", choices=["eigen"], default="eigen",
                   help="Aktuell: nur 'eigen' (Eigenwerte aus Datei).")
    p.add_argument("--eigenfile", type=str, default="eigenvalues_clean.txt")
    p.add_argument("--L", type=float, default=100.0)

    # Trace / Korrelation
    p.add_argument("--tmax", type=float, default=5.0)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--remove_dc", action="store_true",
                   help="Dämpfe den DC-Spike in |Trace| (empfohlen).")

    # Scaling / Fit
    p.add_argument("--auto_rescale", action="store_true",
                   help="Verwende 1/λ̂ als globalen Faktor für die Trace-Analyse.")

    p.add_argument("--plot", action="store_true",
                   help="(Flag ohne Wirkung; Plots werden immer gespeichert.)")

    args = p.parse_args()

    if args.mode == "eigen":
        run_from_eigenfile(args)
    else:
        raise NotImplementedError("Nur --mode eigen ist implementiert.")


if __name__ == "__main__":
    main()
