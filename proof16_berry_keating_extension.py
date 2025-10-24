#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof16_berry_keating_extension.py
--------------------------------------------------------------
Proof 16 — Berry–Keating Extension (patched)

WAS IST NEU (gegenüber deiner Vorversion):
  1) Robuste Pearson-Kreuzkorrelation (z-score, FFT-frei) + Nulltest (Phase-Shuffle)
  2) Optionale DC-Unterdrückung der Quantenspur (--remove_dc)
  3) Energieskalen-Rescaling: E_phys = alpha * E / (L-1)^2   (alpha frei)
  4) Robuster Hochenergie-Fit:
        N(E) ≈ a·E log E + b·E + c·log E + d
     → aus a wird λ = a / (1/(2π)), dann BK-Fit mit E' = E_phys / λ
  5) Plots: Zählfunktion (emp vs. BK-rescaled), Trace-Overlay, Pearson-Korrelation mit Nullgrenze
  6) CSV-Summary mit a,b,c,d, λ, k̂, corr_max, lag_max, Parametern

Aufruf (Windows .bat-Beispiel):
  python proof16_berry_keating_extension.py ^
      --eigenfile eigenvalues_clean.txt ^
      --L 100 --tmax 5 --dt 0.002 ^
      --alpha 1.0 --remove_dc --plot

Autor: Privat/ChatGPT — 2025-10-23
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from typing import Tuple

# --------------------------- Hilfsfunktionen ---------------------------

def load_eigenvalues(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Eigenvalue file not found: {path}")
    vals = np.loadtxt(path, dtype=float).ravel()
    vals = np.sort(vals[~np.isnan(vals)])
    return vals

def counting_function_from_eigs(eigs: np.ndarray, E_grid: np.ndarray) -> np.ndarray:
    # N(E) = # {E_n ≤ E}
    return np.searchsorted(eigs, E_grid, side="right").astype(float)

def N_BK(E: np.ndarray, k: float) -> np.ndarray:
    # Berry–Keating counting (leading)
    E = np.asarray(E, float)
    Em = np.maximum(E, 1e-20)
    return (Em/(2*np.pi))*(np.log(Em) + k)

def fit_N_BK_k(E: np.ndarray, N: np.ndarray, E_min_fit: float = 1.0) -> float:
    mask = (E >= E_min_fit) & np.isfinite(N)
    if mask.sum() < 10:
        # fallback: fit auf allen Punkten
        mask = np.isfinite(N)
    Em = np.maximum(E[mask], 1e-20)
    y = N[mask]
    X = (Em/(2*np.pi))
    # y ≈ X * (log Em + k)  →  k ≈ (y/X) - log Em
    k_est = np.median(y/X - np.log(Em))
    return float(k_est)

def fit_general_model(E: np.ndarray, N: np.ndarray) -> Tuple[float,float,float,float]:
    """
    Fit: N(E) ≈ a*E log E + b*E + c*log E + d  (Least Squares)
    """
    Em = np.maximum(np.asarray(E, float), 1e-20)
    y = np.asarray(N, float)
    X = np.column_stack([Em*np.log(Em), Em, np.log(Em), np.ones_like(Em)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a,b,c,d = coef
    return float(a), float(b), float(c), float(d)

def compute_quantum_trace(eigs: np.ndarray, tmax: float, dt: float, remove_dc: bool) -> Tuple[np.ndarray, np.ndarray]:
    t = np.arange(0.0, tmax, dt)
    trace = np.zeros_like(t, dtype=complex)
    # naive time signal  ∑_n e^{-i E_n t}
    for En in eigs:
        trace += np.exp(-1j*En*t)
    q_abs = np.abs(trace).astype(float)
    if remove_dc:
        # DC-Spike (t≈0) glätten: ersetze die ersten Samples durch Median der Nachbarschaft
        k = min(100, max(10, len(q_abs)//200))
        base = float(np.median(q_abs[1:k]))
        q_abs[:max(3, k//10)] = base
    return t, q_abs

def bk_trace_proxy(t: np.ndarray, tau: float = 0.03, p: float = 2.0) -> np.ndarray:
    """
    Sehr einfacher glatter Proxy für eine semiklassische BK-Spur-Amplitude.
    Liefert |T_BK(t)| ~ (t+τ)^(-p) * exp(-t/τ_c), normiert auf [0,1].
    Das Ziel ist hier NUR eine robuste Form-Referenz für die Korrelation.
    """
    eps = 1e-6
    amp = (t + eps)**(-p) * np.exp(-t/max(tau, 1e-3))
    amp = amp / (np.max(amp) + 1e-12)
    return amp

def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    x = x - np.nanmean(x)
    s = np.nanstd(x) + 1e-12
    return x / s

def pearson_xcorr_with_null(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float, int, float]:
    """
    Standardisierte Kreuzkorrelation (Pearson) + Nulltest (phasen-geshuffeltes y).
    Gibt (lags, corr, cmax, lag_max, null_max) zurück.
    """
    xz = zscore(x)
    yz = zscore(y)

    corr_full = np.correlate(xz, yz, mode="full") / len(xz)
    lags = np.arange(-len(xz)+1, len(xz))
    i_max = int(np.argmax(np.abs(corr_full)))
    cmax = float(np.abs(corr_full[i_max]))
    lag_max = int(lags[i_max])

    # Nulltest: identisches Spektrum, zufällige Phase
    rng = np.random.default_rng(12345)
    Y = np.fft.rfft(yz)
    phases = rng.uniform(0, 2*np.pi, size=Y.size)
    y_surr = np.fft.irfft(np.abs(Y) * np.exp(1j*phases), n=len(yz))
    corr_null = np.correlate(xz, zscore(y_surr), mode="full") / len(xz)
    null_max = float(np.max(np.abs(corr_null)))

    return lags, corr_full, cmax, lag_max, null_max

# --------------------------- Plotter ---------------------------

def plot_counting(E_axis, N_emp, E_axis_bk, N_bk, L, k_hat, lambda_scale):
    plt.figure(figsize=(10, 6.8))
    plt.plot(E_axis, N_emp, label="N_emp(E) (operator)", lw=2)
    plt.plot(E_axis_bk, N_bk, "--", label=f"N_BK(E'; k={k_hat:.3f}),  λ={lambda_scale:.3g}", lw=2)
    plt.xlabel("E"); plt.ylabel("Counting function N(E)")
    plt.title(f"Proof 16 — Counting comparison (L={L})")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"proof16_counting_L{L}.png", dpi=150); plt.close()

def plot_trace_overlay(t, q_abs, bk_abs, L):
    qa = q_abs / (np.max(q_abs) + 1e-12)
    ba = bk_abs / (np.max(bk_abs) + 1e-12)
    plt.figure(figsize=(10, 6.8))
    plt.plot(t, qa, "b-", lw=1.2, label="|Trace e^{-iHt}| (Quantum)")
    plt.plot(t, ba, "r-", lw=1.5, alpha=0.85, label="|Trace_BK| (Semiclassical, proxy)")
    plt.xlabel("t"); plt.ylabel("Normalized amplitude")
    plt.title(f"Proof 16 — Trace overlay (L={L})")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"proof16_trace_overlay_L{L}.png", dpi=150); plt.close()

def plot_correlation(lags, corr, cmax, lag_max, null_max, L):
    plt.figure(figsize=(10.5, 5.6))
    plt.plot(lags, corr, "k-", lw=1.7, label="Pearson xcorr")
    plt.axhline(null_max, color="r", ls="--", lw=1.0, alpha=0.8, label=f"null max ≈ {null_max:.2f}")
    plt.axhline(-null_max, color="r", ls="--", lw=1.0, alpha=0.8)
    plt.title(f"Proof 16 — Quantum vs BK correlation (L={L})\nmax={cmax:.3f} at lag={lag_max}")
    plt.xlabel("Lag (samples)"); plt.ylabel("Correlation (normalized)")
    plt.grid(alpha=0.3); plt.legend(); plt.tight_layout()
    plt.savefig(f"proof16_correlation_L{L}.png", dpi=150); plt.close()

# --------------------------- Hauptanalyse ---------------------------

def main():
    p = argparse.ArgumentParser(description="Proof 16 — Berry–Keating Extension (patched)")
    p.add_argument("--eigenfile", type=str, default="eigenvalues_clean.txt")
    p.add_argument("--L", type=float, default=100.0)
    p.add_argument("--tmax", type=float, default=5.0)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--alpha", type=float, default=1.0, help="energy scale factor in E_phys=alpha*E/(L-1)^2")
    p.add_argument("--remove_dc", action="store_true")
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()

    print("=== Proof 16: Berry–Keating Extension (patched) ===")
    eig_path = Path(args.eigenfile)
    eigs = load_eigenvalues(eig_path)
    print(f"[proof16] Loaded {len(eigs)} eigenvalues from {eig_path}")
    print(f"[proof16] E-min={np.min(eigs):.4g}, E-max={np.max(eigs):.4g}")

    # ---------------- Zählfunktion + Skalen-Rescaling ----------------
    # rohes E-Gitter
    E_grid = np.linspace(max(1e-6, np.min(eigs)*0.9), np.max(eigs)*1.02, 800)
    N_emp = counting_function_from_eigs(eigs, E_grid)

    # physikalische Skala (1D-Länge): E_phys = alpha * E / (L-1)^2
    denom = max((args.L - 1.0)**2, 1e-12)
    E_phys = args.alpha * E_grid / denom

    # robuster allgemeiner Fit im Hochenergiefenster
    E_min_fit = max(1.0, 0.2*np.max(E_phys))
    mask_hi = E_phys >= E_min_fit
    a,b,c,d = fit_general_model(E_phys[mask_hi], N_emp[mask_hi])

    a_target = 1.0/(2*np.pi)
    lambda_scale = (a / a_target) if a > 0 else 1.0
    E_rescaled = E_phys / max(lambda_scale, 1e-12)

    # BK-Fit mit reskalierter Energie (nur k bleibt frei)
    k_hat = fit_N_BK_k(E_rescaled, N_emp, E_min_fit=E_min_fit/max(lambda_scale,1e-12))
    N_bk_rescaled = N_BK(E_rescaled, k_hat)

    print(f"[proof16] a={a:.4e} (target {a_target:.4e}), λ={lambda_scale:.4g}, k̂={k_hat:.3f}")

    # ---------------- Quantenspur & BK-Proxy-Trace ----------------
    t, q_abs = compute_quantum_trace(eigs, args.tmax, args.dt, remove_dc=args.remove_dc)
    bk_abs = bk_trace_proxy(t, tau=0.03, p=2.0)

    # robuste Pearson-Korrelation + Nulltest
    lags, corr, cmax, lag_max, null_max = pearson_xcorr_with_null(q_abs, bk_abs)

    # ---------------- Plots / Export ----------------
    if args.plot:
        plot_counting(E_grid, N_emp, E_grid, N_bk_rescaled, args.L, k_hat, lambda_scale)
        plot_trace_overlay(t, q_abs, bk_abs, args.L)
        plot_correlation(lags, corr, cmax, lag_max, null_max, args.L)

    # Summary CSV
    out = pd.DataFrame([{
        "L": args.L,
        "N_eigs": len(eigs),
        "alpha": args.alpha,
        "remove_dc": bool(args.remove_dc),
        "a": a, "b": b, "c": c, "d": d,
        "a_target": a_target,
        "lambda_scale": lambda_scale,
        "k_hat": k_hat,
        "corr_max": cmax,
        "lag_max": lag_max,
        "null_max": null_max,
        "tmax": args.tmax,
        "dt": args.dt
    }])
    out.to_csv("proof16_summary.csv", index=False)
    print("[proof16] Summary → proof16_summary.csv")

    # Key printout
    print("\n=== Key numbers (patched) ===")
    print(f"• a = {a:.4e}  (target 1/(2π) = {a_target:.4e})")
    print(f"• λ (scale from a) = {lambda_scale:.4g}")
    print(f"• k̂ (BK on E' = E_phys/λ) = {k_hat:.3f}")
    print(f"• corr_max = {cmax:.3f}  (null_max ≈ {null_max:.2f})  @ lag={lag_max}")
    print(f"• remove_dc={args.remove_dc}, alpha={args.alpha}, L={args.L}")

if __name__ == "__main__":
    main()
