#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Proof 18 – Asymptotic Density and Quantum–Statistical Coupling
Author: Helena Kivelitz (2025)

This module extends the operator-zeta framework (Proof 6–17)
to analyze the asymptotic spectral density ρ(E) and coupling
between the empirical spectrum of H[ρ] and the Berry–Keating (BK)
semiclassical prediction.

Input:
    proof17_enhanced_summary.json
    proof17_counting_rescaled.png  (for context)
    eigenvalues_clean.txt

Output:
    proof18_asymptotic_density.png
    proof18_stable_final.png
    proof18_energy_scale_optimization.png
"""

import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path
from scipy.stats import linregress
from scipy.optimize import curve_fit

TWOPI_INV = 1.0 / (2.0 * np.pi)

# ----------------------------------------------------------
# Utilities
# ----------------------------------------------------------

def load_eigs(path):
    """Load and clean positive eigenvalues."""
    E = np.loadtxt(path, dtype=float)
    E = E[np.isfinite(E)]
    return np.sort(E[E > 0.0])

def empirical_density(E, bins=200, smooth=5):
    """Compute normalized empirical density via histogram."""
    hist, edges = np.histogram(E, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # simple moving average smoothing
    if smooth > 1:
        kernel = np.ones(smooth) / smooth
        hist = np.convolve(hist, kernel, mode="same")
    return centers, hist / np.max(hist)

def bk_density(E):
    """Berry–Keating asymptotic density ρ_BK(E) ~ (1/(2π)) log(E/E0)."""
    E = np.asarray(E)
    return TWOPI_INV * np.log(np.maximum(E, 1e-8))

def coupled_density(E, alpha=1.0, beta=0.0):
    """Coupled model between empirical and BK density."""
    return alpha * bk_density(E) + beta * np.sqrt(E) / (1.0 + E)

def residual_metric(rho_emp, rho_bk):
    """L2 residual between empirical and model densities."""
    return np.sqrt(np.mean((rho_emp - rho_bk)**2))

# ----------------------------------------------------------
# Main Analysis
# ----------------------------------------------------------

def analyze_asymptotic_density(eigenfile, summaryfile, outdir="."):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading eigenvalues and summary...")
    E = load_eigs(eigenfile)
    summary = json.load(open(summaryfile))

    lam = summary["fit"]["lam"]
    mu = summary["fit"]["mu"]
    print(f"Using λ={lam:.6f}, μ={mu:.6f} from Proof 17")

    # Compute empirical density
    centers, rho_emp = empirical_density(E, bins=250, smooth=7)
    rho_bk = bk_density(centers)
    rho_coup = coupled_density(centers, alpha=1.0, beta=0.02)

    # Fit empirical density to BK form
    def fit_fn(E, a, b):
        return a * np.log(E) + b
    popt, _ = curve_fit(fit_fn, centers, rho_emp, p0=(TWOPI_INV, 0))
    a_fit, b_fit = popt
    print(f"Fitted asymptotic density: a={a_fit:.6f}, b={b_fit:.6f}")

    # Compute residual metrics
    res_bk = residual_metric(rho_emp, rho_bk)
    res_fit = residual_metric(rho_emp, fit_fn(centers, a_fit, b_fit))
    print(f"Residuals: BK={res_bk:.3e}, Fitted={res_fit:.3e}")

    # ----------------------------------------------------------
    # Figure 1 – Asymptotic Density Comparison
    # ----------------------------------------------------------
    plt.figure(figsize=(8,4))
    plt.plot(centers, rho_emp, label="Empirical ρ(E)", color="#1f77b4", lw=1.5)
    plt.plot(centers, rho_bk, label="BK asymptotic", color="r", lw=1.5)
    plt.plot(centers, fit_fn(centers, a_fit, b_fit), "--", color="k",
             label=f"Fit: a={a_fit:.5f}")
    plt.xlabel("E")
    plt.ylabel("ρ(E) (normalized)")
    plt.title("Proof 18 – Asymptotic Spectral Density")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir/"proof18_asymptotic_density.png", dpi=150)
    plt.close()

    # ----------------------------------------------------------
    # Figure 2 – Energy-Scale Optimization
    # ----------------------------------------------------------
    logE = np.log(np.maximum(centers, 1e-6))
    slope, intercept, r, _, _ = linregress(logE, rho_emp)
    plt.figure(figsize=(8,4))
    plt.scatter(logE, rho_emp, s=8, alpha=0.6, label="Empirical ρ(E)")
    plt.plot(logE, slope*logE+intercept, "r--",
             label=f"Linear fit: slope={slope:.4f}, R²={r**2:.3f}")
    plt.xlabel("log(E)")
    plt.ylabel("ρ(E)")
    plt.title("Proof 18 – Energy Scale Optimization")
    plt.legend(); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir/"proof18_energy_scale_optimization.png", dpi=150)
    plt.close()

    # ----------------------------------------------------------
    # Figure 3 – Stability of Coupling Term
    # ----------------------------------------------------------
    betas = np.linspace(-0.05, 0.05, 40)
    residuals = [residual_metric(rho_emp, coupled_density(centers, 1.0, b))
                 for b in betas]
    plt.figure(figsize=(8,4))
    plt.plot(betas, residuals, color="#2ca02c", lw=1.8)
    plt.axvline(0, color="k", ls="--", lw=1)
    plt.xlabel("β (coupling strength)")
    plt.ylabel("Residual L₂ error")
    plt.title("Proof 18 – Coupling Stability Analysis")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir/"proof18_stable_final.png", dpi=150)
    plt.close()

    # ----------------------------------------------------------
    # Output summary
    # ----------------------------------------------------------
    results = dict(
        lam_fit=float(a_fit),
        b_fit=float(b_fit),
        res_BK=float(res_bk),
        res_fit=float(res_fit),
        slope_logE=float(slope),
        r2=float(r**2),
    )

    outpath = outdir / "proof18_summary.json"
    json.dump(results, open(outpath, "w"), indent=2)
    print(f"Saved summary to {outpath}")
    return results


# ----------------------------------------------------------
# Run standalone
# ----------------------------------------------------------

if __name__ == "__main__":
    analyze_asymptotic_density(
        eigenfile="eigenvalues_clean.txt",
        summaryfile="proof17_enhanced_summary.json",
        outdir="."
    )
