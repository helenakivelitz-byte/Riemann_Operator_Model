#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Proof 17 - Enhanced with Bootstrap and Newton Convergence Analysis
"""

import argparse, math, json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# Constants and Utilities
# -----------------------------

TWOPI_INV = 1.0/(2.0*np.pi)

def load_eigs(path):
    E = np.loadtxt(path, dtype=float)
    E = E[np.isfinite(E)]
    return np.sort(E[E>0.0])

def empirical_counting(Epos):
    """Return sorted energies and cumulative counting N(E)."""
    E = np.sort(Epos)
    N = np.arange(1, len(E)+1, dtype=float)
    return E, N

def N_BK(Ep, lam, c0=-1.0):
    # standard smooth BK main term: lam * E' * (log(E') + c0)
    Ep = np.asarray(Ep, float)
    out = lam * Ep * (np.log(np.maximum(Ep,1e-12)) + c0)
    return out

def N_BK_inv(Ntarget, lam, c0=-1.0, iters=30):
    """
    Invert y = lam*E*(log E + c0) for E via Newton; vectorized & stable.
    Good initial guess with LambertW-like heuristic.
    """
    N = np.asarray(Ntarget, float)
    # initial guess: solve lam*E*logE ~ N  => E0 ~ N / (lam*log(N+e))
    g = np.maximum(N, 1e-9)
    E = g / (lam * np.maximum(np.log(g+math.e), 1.0))
    E = np.maximum(E, 1e-12)

    for _ in range(iters):
        f  = lam*E*(np.log(E)+c0) - N
        fp = lam*(np.log(E)+c0 + 1.0)
        step = f/np.maximum(fp, 1e-12)
        E = np.maximum(E - step, 1e-12)
        if np.max(np.abs(step)) < 1e-10:
            break
    return E

def kde_1d(x, grid, bw=None):
    """Simple Gaussian KDE on grid."""
    x = np.asarray(x)
    grid = np.asarray(grid)
    if bw is None:
        # Scott's rule
        std = np.std(x)
        bw = 1.06*std*np.power(len(x), -1/5) + 1e-12
    diffs = (grid[:,None]-x[None,:])/bw
    K = np.exp(-0.5*diffs*diffs) / (np.sqrt(2*np.pi)*bw)
    dens = K.mean(axis=1)
    return dens

def robust_slice(E, N, lo=0.10, hi=0.95):
    """Keep central quantiles [lo, hi] for robust fitting."""
    m = len(E)
    i0 = int(np.floor(lo*m))
    i1 = int(np.floor(hi*m))
    return E[i0:i1], N[i0:i1]

# -----------------------------
# OT/Quantile Rescale
# -----------------------------

def quantile_rescale(E_raw, lam=TWOPI_INV, c0=-1.0):
    """
    Build a monotone warping E' = f(E) such that
    N_emp(E) == N_BK(E', lam, c0).  (CDF matching)
    """
    E, N = empirical_counting(E_raw)
    Eprime = N_BK_inv(N, lam=lam, c0=c0)
    return Eprime, N, E  # E' aligned to E

# -----------------------------
# Mixed-model fit on rescaled spectrum
# -----------------------------

def design_matrix(Ep):
    Ep = np.asarray(Ep)
    return np.column_stack([
        Ep*np.log(np.maximum(Ep,1e-12)),   # E' log E'
        Ep,                                 # linear
        np.sqrt(Ep),                        # sqrt
        np.ones_like(Ep)                    # constant
    ])

def ridge_fit(Ep, Np, lam_target=TWOPI_INV, ridge=1e-2, weight_lambda=10.0):
    """
    Solve min ||Xθ - N||^2 + ridge*||θ||^2 + weight_lambda*(θ0 - lam_target)^2
    where θ = [λ, μ, ν, c].  (θ0 corresponds to column Ep*log Ep)
    """
    X = design_matrix(Ep)
    y = np.asarray(Np)

    # Normal equations with augmented rows for ridge + λ prior
    A = X.T @ X
    b = X.T @ y

    # ridge
    A += ridge*np.eye(A.shape[0])

    # prior on θ0 ≈ lam_target
    prior = np.zeros(4); prior[0] = lam_target
    P = np.zeros((4,4)); P[0,0] = weight_lambda
    A += P
    b += P @ prior

    theta = np.linalg.solve(A, b)
    return dict(lam=theta[0], mu=theta[1], nu=theta[2], c=theta[3], theta=theta)

def N_model(Ep, pars):
    lam, mu, nu, c = pars["lam"], pars["mu"], pars["nu"], pars["c"]
    return lam*Ep*np.log(np.maximum(Ep,1e-12)) + mu*Ep + nu*np.sqrt(Ep) + c

# -----------------------------
# Bootstrap Error Analysis
# -----------------------------

def bootstrap_error_analysis(Eprime, N, n_boot=2000, ridge=1e-2, lambda_weight=10.0):
    """Bootstrap estimates of parameter uncertainties"""
    n_samples = len(Eprime)
    indices = np.random.choice(n_samples, (n_boot, n_samples), replace=True)
    params_boot = []
    
    for i in range(n_boot):
        idx = indices[i]
        E_boot = Eprime[idx]
        N_boot = N[idx]
        try:
            pars = ridge_fit(E_boot, N_boot, lam_target=TWOPI_INV,
                           ridge=ridge, weight_lambda=lambda_weight)
            params_boot.append([pars['lam'], pars['mu'], pars['nu'], pars['c']])
        except np.linalg.LinAlgError:
            continue
    
    params_boot = np.array(params_boot)
    means = np.mean(params_boot, axis=0)
    errors = np.std(params_boot, axis=0)
    return errors, means, params_boot

def plot_bootstrap_distribution(params_boot, outdir):
    """Plot bootstrap distribution of parameters"""
    # Ensure output directory exists
    outdir.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    param_names = ['λ', 'μ', 'ν', 'c']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for i, (ax, name, color) in enumerate(zip(axes.flat, param_names, colors)):
        values = params_boot[:, i]
        ax.hist(values, bins=50, alpha=0.7, color=color, density=True)
        ax.axvline(np.mean(values), color='k', linestyle='--', linewidth=2, 
                  label=f'Mean: {np.mean(values):.6f}')
        if i == 0:
            ax.axvline(TWOPI_INV, color='r', linestyle=':', linewidth=2, label='Target 1/(2π)')
        else:
            ax.axvline(0, color='r', linestyle=':', linewidth=2, label='Zero')
        ax.set_xlabel(f'{name} value')
        ax.set_ylabel('Density')
        ax.legend()
        ax.set_title(f'Bootstrap: {name} = {np.mean(values):.6f} ± {np.std(values):.6f}')
    
    plt.tight_layout()
    plt.savefig(outdir / "proof17_bootstrap_dist.png", dpi=150)
    plt.close()

# -----------------------------
# Newton Convergence Analysis
# -----------------------------

def analyze_newton_convergence(E_target, N_target, lam=TWOPI_INV, c0=-1.0, max_iter=10):
    """Analyze Newton iteration convergence for N_BK inverse"""
    errors = []
    E_guesses = []
    
    # Initial guess
    E = N_target / (lam * np.maximum(np.log(N_target + np.e), 1.0))
    E = np.maximum(E, 1e-12)
    E_guesses.append(E)
    
    for iter in range(max_iter):
        f = lam * E * (np.log(E) + c0) - N_target
        fp = lam * (np.log(E) + c0 + 1.0)
        step = f / np.maximum(fp, 1e-12)
        E_new = np.maximum(E - step, 1e-12)
        
        error = np.abs(step)
        errors.append(error)
        E_guesses.append(E_new)
        
        E = E_new
        if error < 1e-12:
            break
    
    return np.array(errors), np.array(E_guesses)

def plot_newton_convergence(Eprime, N, E_raw, outdir, n_samples=5):
    """Plot Newton convergence for representative samples"""
    # Ensure output directory exists
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Select representative points across the spectrum
    indices = np.linspace(0, len(Eprime)-1, n_samples, dtype=int)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot convergence for selected points
    colors = plt.cm.viridis(np.linspace(0, 1, n_samples))
    
    for idx, color in zip(indices, colors):
        E_val = E_raw[idx]  # Original eigenvalue
        N_val = N[idx]      # Empirical counting
        
        errors, E_guesses = analyze_newton_convergence(E_val, N_val)
        iterations = np.arange(len(errors))
        
        ax1.semilogy(iterations, errors, 'o-', color=color, 
                    label=f'E={E_val:.1f}', alpha=0.8, markersize=4)
        ax2.plot(iterations, E_guesses[:-1], 'o--', color=color,
                label=f'E={E_val:.1f}', alpha=0.6, markersize=4)
    
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Step Size (log scale)')
    ax1.set_title('Newton Convergence: Error vs Iteration')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('E guess')
    ax2.set_title('Newton Convergence: E guess evolution')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(outdir / "proof17_newton_convergence.png", dpi=150)
    plt.close()

# -----------------------------
# Original Plots (slightly enhanced)
# -----------------------------

def plot_all(E_raw, Eprime, N, pars, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists

    # 1) Residuals N_emp - N_model(E')
    res = N - N_model(Eprime, pars)
    plt.figure(figsize=(10,3.6))
    plt.plot(Eprime, res, lw=1.5, color='#2e86ab')
    plt.axhline(0, ls='--', c='k', lw=1)
    plt.xlabel("E' (rescaled)")
    plt.ylabel("Residual  N_emp - N_model")
    plt.title("Proof 17 — Residuals after OT rescaling")
    plt.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(outdir/"proof17_residuals.png"); plt.close()

    # 2) Density + local slope a(E)
    g = np.linspace(max(Eprime.min()*0.9, 1e-6), Eprime.max()*1.05, 400)
    rho_emp = kde_1d(Eprime, g)
    rho_emp /= rho_emp.max() + 1e-12
    rho_bk = (TWOPI_INV*(np.log(np.maximum(g,1e-12))+1.0))
    rho_bk = rho_bk / (rho_bk.max()+1e-12)

    fig,ax = plt.subplots(1,2, figsize=(11,4))
    ax[0].plot(g, rho_bk, 'r-', label='BK density (norm.)', alpha=0.8)
    ax[0].plot(g, rho_emp, color='#3b83bd', label='Empirical KDE (norm.)', alpha=0.8)
    ax[0].set_xlabel("E' (rescaled)"); ax[0].set_ylabel("ρ(E') (norm.)")
    ax[0].legend(); ax[0].set_title("Spectral density after rescaling")
    ax[0].grid(True, alpha=0.3)

    # local slope a(E) := dN/dE' / log term approx
    dN_dE = np.gradient(N_model(g, pars), g)
    a_loc = dN_dE / np.maximum(np.log(np.maximum(g,1e-12))+1.0, 1e-12)
    ax[1].plot(g, a_loc, 'k-', label='local a(E)', alpha=0.8)
    ax[1].axhline(TWOPI_INV, ls='--', c='orange', label='1/(2π)')
    ax[1].set_xlabel("E'"); ax[1].set_ylabel('a(E)')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    ax[1].set_title("Local slope diagnostics")
    plt.tight_layout(); plt.savefig(outdir/"proof17_density_rescaled.png"); plt.close(fig)

    # 3) Counting overlay
    g_fine = np.linspace(Eprime.min()*0.95, Eprime.max()*1.05, 500)
    plt.figure(figsize=(10,4))
    plt.scatter(Eprime, N, s=8, alpha=0.6, color='#1f77b4', label="N_emp(E') (rescaled)")
    plt.plot(g_fine, N_model(g_fine, pars), 'k-', lw=2, alpha=0.8,
             label=f"N_model(E')")
    plt.xlabel("E' (rescaled energy)"); plt.ylabel("Counting function N(E')")
    plt.title(f"Proof 17 — Counting after OT rescaling (λ={pars['lam']:.6f})")
    plt.legend(); plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(outdir/"proof17_counting_rescaled.png"); plt.close()

# -----------------------------
# Enhanced Main Function
# -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eigenfile", default="eigenvalues_clean.txt")
    ap.add_argument("--lam0", type=float, default=TWOPI_INV, help="target λ ≈ 1/(2π)")
    ap.add_argument("--c0", type=float, default=-1.0, help="BK constant in E logE + c0 E")
    ap.add_argument("--trim_lo", type=float, default=0.10)
    ap.add_argument("--trim_hi", type=float, default=0.95)
    ap.add_argument("--ridge", type=float, default=1e-2)
    ap.add_argument("--lambda_weight", type=float, default=10.0)
    ap.add_argument("--n_bootstrap", type=int, default=2000)
    ap.add_argument("--outdir", default="results_proof17")
    args = ap.parse_args()

    # Create output directory at the very beginning
    outdir_path = Path(args.outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {outdir_path.absolute()}")

    print("Loading eigenvalues...")
    Eraw = load_eigs(args.eigenfile)
    if len(Eraw) < 50:
        raise SystemExit("Not enough positive eigenvalues.")
    print(f"Loaded {len(Eraw)} eigenvalues")

    # 1) OT/Quantile rescale to match counting
    print("Performing OT quantile rescale...")
    Eprime_all, N_all, E_sorted = quantile_rescale(Eraw, lam=args.lam0, c0=args.c0)

    # 2) Robust central slice for stable fit
    Eprime, N = robust_slice(Eprime_all, N_all, args.trim_lo, args.trim_hi)
    print(f"Using {len(Eprime)} points for fitting ({args.trim_lo*100:.0f}%-{args.trim_hi*100:.0f}%)")

    # 3) Mixed model (with prior on λ → 1/(2π))
    print("Fitting mixed model...")
    pars = ridge_fit(Eprime, N, lam_target=args.lam0,
                     ridge=args.ridge, weight_lambda=args.lambda_weight)

    # 4) Bootstrap analysis
    print(f"Running bootstrap analysis (n={args.n_bootstrap})...")
    errors, means, params_boot = bootstrap_error_analysis(
        Eprime, N, n_boot=args.n_bootstrap, ridge=args.ridge, 
        lambda_weight=args.lambda_weight)

    # 5) Newton convergence analysis
    print("Analyzing Newton convergence...")
    plot_newton_convergence(Eprime, N_all, E_sorted, outdir=outdir_path)

    # 6) Plot bootstrap distributions
    plot_bootstrap_distribution(params_boot, outdir=outdir_path)

    # 7) Original plots
    print("Generating plots...")
    plot_all(Eraw, Eprime_all, N_all, pars, outdir=args.outdir)

    # 8) Enhanced summary with bootstrap results
    summary = {
        "lam_target": args.lam0,
        "fit": {k: float(v) for k,v in pars.items() if k != "theta"},
        "bootstrap": {
            "lam": {"mean": float(means[0]), "std": float(errors[0])},
            "mu": {"mean": float(means[1]), "std": float(errors[1])},
            "nu": {"mean": float(means[2]), "std": float(errors[2])},
            "c": {"mean": float(means[3]), "std": float(errors[3])},
            "n_boot": args.n_bootstrap
        },
        "n_eigs": int(len(Eraw)),
        "trim": [args.trim_lo, args.trim_hi],
        "c0": args.c0
    }
    
    summary_file = outdir_path / "proof17_enhanced_summary.json"
    summary_file.write_text(json.dumps(summary, indent=2))
    
    print("\n=== Proof 17 (Enhanced with Bootstrap) ===")
    print(f"λ = {means[0]:.6f} ± {errors[0]:.6f}  (target: {TWOPI_INV:.6f})")
    print(f"μ = {means[1]:.6f} ± {errors[1]:.6f}")
    print(f"ν = {means[2]:.6f} ± {errors[2]:.6f}")
    print(f"c = {means[3]:.6f} ± {errors[3]:.6f}")
    print(f"\nBootstrap completed with {args.n_bootstrap} replicates")
    print(f"Results saved to {outdir_path.absolute()}")

if __name__ == "__main__":
    main()
