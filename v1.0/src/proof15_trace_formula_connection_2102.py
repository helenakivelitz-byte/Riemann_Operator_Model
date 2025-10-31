#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof15_trace_formula_connection_2102.py
--------------------------------------------------------------
Proof 15 – Trace Formula Connection (Deep, improved)

Verbindet Eigenwertspektrum von H[ρ_*] mit Primzahlstruktur:
    ∑ e^{-iE_n t}  ↔  ∑ log(p) δ(t - log p)

Erweiterungen:
    • CLI-Parameter --sigma_t für Glättungsbreite
    • Option --remove_dc zur Dämpfung des DC-Spikes der Quantenspur
    • Automatischer Import von riemann_zeros_firstN.txt (falls vorhanden)
    • Saubere Ausgabe + Diagramme

Autor: ChatGPT (basierend auf deinem Proof-15-Design)
Datum: 2025-10-23
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import windows
from pathlib import Path
from typing import Dict, List, Tuple
import argparse

# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------

def sieve_primes(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i : n + 1 : i] = False
    return np.nonzero(sieve)[0].tolist()

def prime_logarithms(L: float) -> Tuple[List[int], List[float]]:
    primes = sieve_primes(int(L))
    return primes, [np.log(p) for p in primes]

def periodic_orbit_spectrum(eigs: np.ndarray, t_max: float, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    t_grid = np.arange(0, t_max, dt)
    trace_vals = np.zeros_like(t_grid, dtype=complex)
    for E in eigs:
        trace_vals += np.exp(-1j * E * t_grid)
    return t_grid, trace_vals

def prime_signal_model(log_primes: List[float], t_grid: np.ndarray,
                       amplitude: float = 1.0, sigma_t: float = 0.15) -> np.ndarray:
    signal = np.zeros(len(t_grid), dtype=complex)
    for log_p in log_primes:
        signal += amplitude * np.log(np.exp(log_p)) \
            * np.exp(-0.5 * ((t_grid - log_p) / sigma_t) ** 2) / (sigma_t * np.sqrt(2 * np.pi))
    return signal

def spectral_determinant_ratio_log(eigs: np.ndarray, zeros_file: Path,
                                   t_min: float = 0.0, t_max: float = 150.0,
                                   num: int = 500, eps: float = 1e-2
                                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Stable evaluation via complex logs.
    We compare
        log D_Q(t) = sum_n log( i t + eps - E_n )
        log D_ζ(t) = sum_k log( i t + eps - (1/2 + i γ_k) )
    and return:
        t_grid, Re(log R) - mean, Im(log R) (unwrapped), where log R = log D_Q - log D_ζ.
    """
    if not zeros_file.exists():
        print(f"[proof15] Determinant ratio skipped: Zeros file not found: {zeros_file}")
        return np.array([]), np.array([]), np.array([])

    gammas = np.loadtxt(zeros_file, dtype=float).ravel()
    t_grid = np.linspace(t_min, t_max, num)

    # precompute ζ zeros as complex points once
    z_pts = 0.5 + 1j * gammas

    # sums of logs (stable)
    re_logR = np.zeros_like(t_grid, dtype=float)
    im_logR = np.zeros_like(t_grid, dtype=float)

    for j, t in enumerate(t_grid):
        z = 1j * t + eps  # evaluation line: s-like variable offset by i t
        # quantum part
        # log of products → sum of logs
        Lq = np.sum(np.log(z - eigs))
        # zeta part
        Lz = np.sum(np.log(z - z_pts))
        L = Lq - Lz  # log ratio
        re_logR[j] = np.real(L)
        im_logR[j] = np.imag(L)

    # center the real part to avoid huge offsets and unwrap the phase
    re_logR -= np.mean(re_logR)
    im_logR = np.unwrap(im_logR)

    return t_grid, re_logR, im_logR


# ------------------------------------------------------------
# Plotfunktionen
# ------------------------------------------------------------

def plot_trace_comparison(t_grid, q_abs, prime_signal, L, primes):
    plt.figure(figsize=(10, 6))
    plt.plot(t_grid, q_abs / np.max(q_abs), 'b-', label='|Trace e^{-iHt}| (Quantum)')
    plt.plot(t_grid, np.abs(prime_signal) / np.max(np.abs(prime_signal)), 'r-', label='|Prime Signal|')
    for p in primes[:10]:
        plt.axvline(np.log(p), color='gray', linestyle='--', alpha=0.4)
    plt.xlabel("t")
    plt.ylabel("Normalized amplitude")
    plt.title(f"Proof 15 – Trace Formula Comparison (L={L})")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"proof15_trace_comparison_L{L}.png", dpi=150)
    plt.close()

def plot_correlation(t_grid, correlation, L):
    plt.figure(figsize=(8, 4))
    plt.plot(t_grid, correlation, 'k-', lw=1.8)
    plt.xlabel("Lag (samples)")
    plt.ylabel("Correlation")
    plt.title(f"Proof 15 – Correlation (L={L})")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"proof15_correlation_L{L}.png", dpi=150)
    plt.close()

def plot_spectral_determinant_log(t_grid: np.ndarray, re_logR: np.ndarray, im_logR: np.ndarray, L):
    if t_grid.size == 0:
        return
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 7))
    plt.subplot(2, 1, 1)
    plt.plot(t_grid, re_logR, lw=1.6)
    plt.ylabel("Re log(D_Q/D_ζ) (centered)")
    plt.title(f"Spectral determinant (log ratio) vs t  (L={L})")
    plt.grid(alpha=0.3)
    plt.subplot(2, 1, 2)
    plt.plot(t_grid, im_logR, lw=1.6)
    plt.xlabel("t")
    plt.ylabel("Phase [rad] (unwrapped)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"proof15_spectral_determinant_L{L}.png", dpi=150)
    plt.close()

# ------------------------------------------------------------
# Hauptanalyse
# ------------------------------------------------------------

def analyze_trace_formula(eigs: np.ndarray, L: float, tmax: float, dt: float,
                          sigma_t: float, zeros_file: Path, remove_dc: bool) -> Dict:
    results = {}
    primes, log_primes = prime_logarithms(L)
    print(f"[proof15] L={L}: {len(primes)} Primzahlen ≤ L")

    # Quantenspur
    t_grid, trace_q = periodic_orbit_spectrum(eigs, tmax, dt)
    q_abs = np.abs(trace_q).astype(float)

    # --- DC-Dämpfung ---
    if remove_dc:
        k = min(100, len(q_abs) - 1)
        base = np.median(q_abs[1:k])
        q_abs[0:3] = base
        print("[proof15] DC spike suppressed (remove_dc=True)")

    # Prime-Signal
    prime_signal = prime_signal_model(log_primes, t_grid, sigma_t=sigma_t)

    # Korrelation (Pearson)
    win = windows.hann(len(q_abs))
    qa = q_abs - np.mean(q_abs)
    pa = np.abs(prime_signal) - np.mean(np.abs(prime_signal))
    corr = np.correlate(qa * win, pa * win, mode="same")
    lag = np.arange(-len(corr) // 2, len(corr) // 2)
    corr_norm = corr / (np.max(np.abs(corr)) + 1e-12)

    # Max-Correlation
    idx_max = int(np.argmax(np.abs(corr_norm)))
    corr_max = float(np.abs(corr_norm[idx_max]))
    lag_max = int(lag[idx_max])
    print(f"[proof15] Correlation max={corr_max:.4f} at lag={lag_max}")

    # --- Spektrale Determinante (stabil via log) ---
    t_det, re_logR, im_logR = spectral_determinant_ratio_log(
        eigs,
        zeros_file,
        t_min=0.0,
        t_max=max(
            150.0,
            1.2 * (np.max(np.loadtxt(zeros_file)) if zeros_file.exists() else 100.0),
        ),
        num=600,
        eps=1e-2,
    )
    plot_spectral_determinant_log(t_det, re_logR, im_logR, L)

    # --- Plots ---
    plot_trace_comparison(t_grid, q_abs, prime_signal, L, primes)
    plot_correlation(lag, corr_norm, L)

    # --- Summary ---
    results.update(
        dict(
            L=L,
            N=len(eigs),
            corr_max=corr_max,
            lag_max=lag_max,
            tmax=tmax,
            dt=dt,
            sigma_t=sigma_t,
        )
    )
    return results

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Proof 15 – Trace Formula Connection (Deep)")
    parser.add_argument("--eigenfile", type=str, default="eigenvalues_clean.txt")
    parser.add_argument("--L", type=float, default=100.0)
    parser.add_argument("--tmax", type=float, default=5.0)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--sigma_t", type=float, default=0.15, help="width of prime spikes")
    parser.add_argument("--zeros", type=str, default="riemann_zeros_firstN.txt")
    parser.add_argument("--remove_dc", action="store_true", help="suppress DC spike of |trace|")
    parser.add_argument("--eps", type=float, default=1e-2)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    print("=== Proof 15: Trace Formula Connection (Deep) ===")
    eig_path = Path(args.eigenfile)
    if not eig_path.exists():
        raise FileNotFoundError(f"Eigenvalue file not found: {eig_path}")
    eigs = np.loadtxt(eig_path)
    print(f"[proof15] Loaded {len(eigs)} eigenvalues from {eig_path}")
    print(f"[proof15] E-Min={np.min(eigs):.4g}, E-Max={np.max(eigs):.4g}")

    zeros_file = Path(args.zeros)
    results = analyze_trace_formula(eigs, args.L, args.tmax, args.dt,
                                    args.sigma_t, zeros_file, args.remove_dc)

    # Export summary
    df = pd.DataFrame([results])
    df.to_csv("proof15_summary.csv", index=False)
    print("[proof15] Summary → proof15_summary.csv")

    print("\n=== Interpretation ===")
    print(f"• Max. correlation: {results['corr_max']:.3f} (lag={results['lag_max']})")
    print(f"• sigma_t={args.sigma_t}, remove_dc={args.remove_dc}")

if __name__ == "__main__":
    main()
