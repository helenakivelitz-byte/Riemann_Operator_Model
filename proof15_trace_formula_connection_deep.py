#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof15_trace_formula_robust.py
--------------------------------------------------------------
Proof 15 - Robust & Normalized Trace Formula Analysis

Korrekturen:
1. Pearson-Korrelation (z-score normalized)
2. Hann-Window für reduzierte Spectral Leakage
3. Log-Domain spektrale Determinante mit Dämpfung
4. Echte Riemann-Nullstellen Integration
5. Korrekte Primzahl-Gewichtung (log p / √p)
6. Nyquist-konforme Zeitabtastung
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Optional

# Robustere Primzahl-Funktionen
def sieve_of_eratosthenes(n: int) -> List[int]:
    """Return all primes <= n"""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i*i:n+1:i] = [False] * len(range(i*i, n+1, i))
    return [i for i, is_prime in enumerate(sieve) if is_prime]

def prime_logarithms(L_max: int) -> Tuple[List[int], List[float]]:
    """Return primes <= L_max and their natural logarithms"""
    primes = sieve_of_eratosthenes(int(L_max))
    log_primes = [np.log(p) for p in primes]
    return primes, log_primes

# KORREKTUR 1: Normalisierte Pearson-Korrelation + Hann-Window
def hann_window(x: np.ndarray) -> np.ndarray:
    """Hann window for reduced spectral leakage"""
    n = len(x)
    if n < 2:
        return x
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / (n - 1))
    return x * w

def xcorr_pearson(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Zero-mean, unit-variance cross-correlation (Pearson)"""
    # Normalisiere und fenstere
    xz = (x - np.mean(x)) / (np.std(x) + 1e-12)
    yz = (y - np.mean(y)) / (np.std(y) + 1e-12)
    xw, yw = hann_window(xz), hann_window(yz)
    
    # Kreuzkorrelation
    corr = np.correlate(xw, yw, mode='full') / len(xw)
    lags = np.arange(-len(xw) + 1, len(xw))
    
    return lags, corr

# KORREKTUR 2: Nyquist-konforme Zeitabtastung
def periodic_orbit_spectrum_nyquist(eigenvalues: np.ndarray, t_max: float = None, 
                                  dt: float = None) -> Tuple[np.ndarray, np.ndarray]:
    """Bandbegrenzte Spur-Berechnung mit Nyquist-Abtastung"""
    E_min, E_max = np.min(eigenvalues), np.max(eigenvalues)
    dE = max(E_max - E_min, 1e-9)
    
    # Nyquist-Abtastung
    if dt is None:
        dt = np.pi / dE  # Nyquist-artig
    
    if t_max is None:
        t_max = 2 * np.log(max(3, int(E_max)))  # Bis ~2*log(p_max)
    
    t_grid = np.arange(0.0, t_max + 1e-12, dt)
    
    # Stabil berechnete Spur
    trace_vals = np.zeros(len(t_grid), dtype=complex)
    for E_n in eigenvalues:
        trace_vals += np.exp(-1j * E_n * t_grid)
    
    return t_grid, trace_vals

# KORREKTUR 3: Korrekte Primzahl-Gewichtung
def prime_signal_model_weighted(log_primes: List[float], t_grid: np.ndarray, 
                              weight_type: str = "explicit", sigma_t: float = 0.2) -> np.ndarray:
    """
    Primzahl-Signal mit korrekter Gewichtung aus expliziter Formel
    """
    signal = np.zeros(len(t_grid), dtype=complex)
    
    for log_p in log_primes:
        p = np.exp(log_p)
        
        # Gewichtung entsprechend expliziter Formel
        if weight_type == "explicit":
            amp = np.log(p) / np.sqrt(p)  # log(p) * p^{-1/2}
        elif weight_type == "logonly":
            amp = np.log(p)
        else:  # "unit"
            amp = 1.0
        
        # Geglättete Delta-Funktion
        gaussian = np.exp(-0.5 * ((t_grid - log_p) / sigma_t)**2)
        gaussian /= (sigma_t * np.sqrt(2 * np.pi))
        
        signal += amp * gaussian
    
    return signal

# KORREKTUR 4: Log-Domain spektrale Determinante mit echten Nullstellen
def load_riemann_zeros(path: str = "riemann_zeros_first1000.txt", maxE: float = None) -> np.ndarray:
    """Lade echte Riemann-Nullstellen (nur positive imaginäre Teile)"""
    if Path(path).exists():
        zeros = np.loadtxt(path)
        if maxE is not None:
            zeros = zeros[zeros <= maxE]
        return zeros
    else:
        # Fallback: Erste paar Nullstellen (nur für Test)
        print(f"[proof15-robust] Warnung: {path} nicht gefunden, verwende Test-Nullstellen")
        return np.array([14.1347, 21.0220, 25.0108, 30.4249, 32.9351, 37.5862])

def log_spectral_determinant(eigenvalues: np.ndarray, E_points: np.ndarray, 
                           zeros: np.ndarray, eps: float = 1e-2) -> np.ndarray:
    """
    Stabile Berechnung von log(Det(E - H) / Det(E - H_ζ))
    mit Regularisierung E → E + iε
    """
    log_ratio = np.zeros_like(E_points, dtype=complex)
    
    for i, E in enumerate(E_points):
        z = E + 1j * eps  # Regularisierung
        
        # Log-Determinante für Quanten-Operator
        with np.errstate(divide='ignore', invalid='ignore'):
            log_det_H = np.sum(np.log(z - eigenvalues + 1e-15j))
            log_det_ζ = np.sum(np.log(z - zeros + 1e-15j))
        
        log_ratio[i] = log_det_H - log_det_ζ
    
    return log_ratio

# Hauptanalyse mit allen Korrekturen
def analyze_trace_formula_robust(eigenvalues: np.ndarray, L: float, 
                               make_plots: bool = True) -> Dict:
    """
    Robustere Trace-Formula Analyse mit allen Verbesserungen
    """
    results = {}
    
    # 1. Primzahl-Logarithmen
    primes, log_primes = prime_logarithms(int(L))
    results['primes'] = primes
    results['log_primes'] = log_primes
    results['n_primes'] = len(primes)
    
    print(f"[proof15-robust] L={L}, {len(primes)} Primzahlen")
    
    # 2. Nyquist-konforme Spur-Berechnung
    t_grid, trace_quantum = periodic_orbit_spectrum_nyquist(eigenvalues)
    results['t_grid'] = t_grid
    results['trace_quantum'] = trace_quantum
    
    # 3. Korrekt gewichtetes Primzahl-Signal
    prime_signal = prime_signal_model_weighted(log_primes, t_grid, weight_type="explicit")
    results['prime_signal'] = prime_signal
    
    # 4. Pearson-Korrelation mit Hann-Window
    lags, correlation = xcorr_pearson(np.abs(trace_quantum), np.abs(prime_signal))
    peak_idx = np.argmax(np.abs(correlation))
    peak_corr = correlation[peak_idx]
    peak_lag = lags[peak_idx] * (t_grid[1] - t_grid[0]) if len(t_grid) > 1 else 0
    
    results['correlation'] = correlation
    results['lags'] = lags
    results['peak_correlation'] = peak_corr
    results['peak_lag'] = peak_lag
    
    # 5. Stabile spektrale Determinante mit echten Nullstellen
    zeros = load_riemann_zeros("riemann_zeros_first1000.txt", maxE=100.0)
    E_test = np.linspace(0.5, 50.0, 200)
    log_det_ratio = log_spectral_determinant(eigenvalues, E_test, zeros)
    results['E_test'] = E_test
    results['log_det_ratio'] = log_det_ratio
    results['n_zeros_used'] = len(zeros)
    
    # 6. Zusätzliche Diagnostik
    results['trace_power'] = np.mean(np.abs(trace_quantum)**2)
    results['prime_power'] = np.mean(np.abs(prime_signal)**2)
    results['power_ratio'] = results['trace_power'] / (results['prime_power'] + 1e-15)
    
    # Plots
    if make_plots:
        plot_robust_analysis(results, L)
    
    return results

def plot_robust_analysis(results: Dict, L: float):
    """Verbesserte Visualisierung der robusten Analyse"""
    plt.figure(figsize=(15, 12))
    
    # Plot 1: Signalvergleich mit korrekter Gewichtung
    plt.subplot(3, 2, 1)
    plt.plot(results['t_grid'], np.abs(results['trace_quantum']), 'b-', alpha=0.8, 
             linewidth=1.5, label='|Spur e^{-iHt}| (Quanten)')
    plt.plot(results['t_grid'], np.abs(results['prime_signal']), 'r-', alpha=0.8, 
             linewidth=1.5, label='|Primzahl-Signal| (log p / √p)')
    
    # Primzahl-Positionen markieren
    for i, (p, log_p) in enumerate(zip(results['primes'], results['log_primes'])):
        if i < 8:  # Nur erste 8 beschriften
            plt.axvline(x=log_p, color='green', linestyle='--', alpha=0.6, 
                       label=f'log({p})' if i < 3 else "")
        else:
            plt.axvline(x=log_p, color='green', linestyle='--', alpha=0.3)
    
    plt.xlabel('Zeit t')
    plt.ylabel('Amplitude')
    plt.title(f'Proof 15 Robust - Gewichtete Signale (L={L})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 2: Pearson-Korrelation
    plt.subplot(3, 2, 2)
    plt.plot(results['lags'], results['correlation'], 'k-', linewidth=2)
    plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='Zero lag')
    plt.xlabel('Lag')
    plt.ylabel('Pearson-Korrelation')
    plt.title(f'Korrelation (Peak: {results["peak_correlation"]:.4f})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 3: Log-spektrale Determinante
    plt.subplot(3, 2, 3)
    plt.plot(results['E_test'], np.real(results['log_det_ratio']), 'g-', 
             label='Real(log Det-Ratio)')
    plt.plot(results['E_test'], np.imag(results['log_det_ratio']), 'purple', 
             label='Imag(log Det-Ratio)')
    plt.xlabel('Energie E')
    plt.ylabel('log(Det_H / Det_ζ)')
    plt.title(f'Log-Spektrale Determinante ({results["n_zeros_used"]} Nullstellen)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 4: Phasenvergleich
    plt.subplot(3, 2, 4)
    plt.plot(results['t_grid'], np.angle(results['trace_quantum']), 'b-', alpha=0.7, 
             label='Phase(Quanten)')
    plt.plot(results['t_grid'], np.angle(results['prime_signal']), 'r-', alpha=0.7, 
             label='Phase(Primzahlen)')
    plt.xlabel('Zeit t')
    plt.ylabel('Phase [rad]')
    plt.title('Phasenvergleich')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 5: Signal-Power Analyse
    plt.subplot(3, 2, 5)
    powers = [results['trace_power'], results['prime_power']]
    labels = ['Quanten-Signal', 'Primzahl-Signal']
    colors = ['blue', 'red']
    bars = plt.bar(labels, powers, alpha=0.7, color=colors)
    plt.ylabel('Mittlere Signal-Power')
    plt.title(f'Power-Verhältnis: {results["power_ratio"]:.3f}')
    
    # Plot 6: Korrelations-Histogramm
    plt.subplot(3, 2, 6)
    plt.hist(results['correlation'], bins=50, alpha=0.7, density=True)
    plt.axvline(x=results['peak_correlation'], color='red', linestyle='--', 
               label=f'Peak: {results["peak_correlation"]:.4f}')
    plt.xlabel('Korrelations-Wert')
    plt.ylabel('Dichte')
    plt.title('Korrelations-Verteilung')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'proof15_robust_analysis_L{L}.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"[proof15-robust] Plots gespeichert: proof15_robust_analysis_L{L}.png")

# Hauptprogramm
def main():
    parser = argparse.ArgumentParser(description='Proof 15 Robust - Verbesserte Trace Formula Analysis')
    parser.add_argument('--eigenfile', type=str, default='eigenvalues_H_rho_final.txt')
    parser.add_argument('--L', type=float, default=12.0)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    
    print("=== Proof 15 Robust: Verbesserte Trace Formula Connection ===")
    
    # Eigenwerte laden
    if not Path(args.eigenfile).exists():
        print(f"Fehler: {args.eigenfile} nicht gefunden!")
        return
    
    eigenvalues = np.loadtxt(args.eigenfile)
    print(f"[proof15-robust] {len(eigenvalues)} Eigenwerte geladen")
    
    # Robuste Analyse
    results = analyze_trace_formula_robust(
        eigenvalues, 
        L=args.L,
        make_plots=not args.no_plots
    )
    
    # Detaillierte Auswertung
    print("\n=== Proof 15 Robust - Quantitative Auswertung ===")
    print(f"Primzahlen: {results['n_primes']} ≤ {args.L}")
    print(f"Pearson-Korrelation (Peak): {results['peak_correlation']:.6f}")
    print(f"Peak-Lag: {results['peak_lag']:.6f}")
    print(f"Signal-Power Verhältnis: {results['power_ratio']:.6f}")
    print(f"Verwendete ζ-Nullstellen: {results['n_zeros_used']}")
    print(f"Primzahl-Gewichtung: log(p) / √p")
    
    # Interpretation
    print("\n=== Interpretation ===")
    if results['peak_correlation'] > 0.7:
        print("✅ AUSGEZEICHNETE Trace-Formula Connection!")
        print("   Starke Evidenz für Primzahl-Struktur in Quanten-Orbits.")
    elif results['peak_correlation'] > 0.4:
        print("✅ GUTE Trace-Formula Connection")
        print("   Klare Primzahl-Struktur erkennbar.")
    else:
        print("⚠️  Moderate Evidenz - weitere Optimierung möglich")
    
    # Ergebnisse speichern
    summary_df = pd.DataFrame({
        'L': [args.L],
        'n_primes': [results['n_primes']],
        'pearson_correlation': [results['peak_correlation']],
        'peak_lag': [results['peak_lag']],
        'power_ratio': [results['power_ratio']],
        'n_zeros_used': [results['n_zeros_used']],
        'n_eigenvalues': [len(eigenvalues)]
    })
    summary_df.to_csv('proof15_robust_summary.csv', index=False)
    print(f"\n[proof15-robust] Zusammenfassung → proof15_robust_summary.csv")

if __name__ == "__main__":
    main()
