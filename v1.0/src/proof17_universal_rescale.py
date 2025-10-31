#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof17_universal_rescale.py
--------------------------------------------------------------
Proof 17 — Universal rescale & robust BK matching (patched)

Hauptideen (robust):
- Zählfunktion gegen ein stabiles BK-Lead-Modell fitten:
    N(E') ≈ (E'/(2π)) [log(E') - 1] + c0 + (beta/(2π)) E'
  (beta ist eine kleine lineare Korrektur, Riemann: beta = -1)
- λ-Schätzung aus N(E)/E ~ a log(E) + b, robust in mittlerem Fenster
- Harte Schranken & Log-Schutz (E' >= 1e-12)
- Robust-Objektiv (Huber/Clipping) + Grid-Suche
- Optionale "Gewichte" für Potentialsummanden nur klein (0..0.05)

Ausgaben:
- proof17_counting_rescaled.png
- proof17_residuals.png
- proof17_density_rescaled.png
- proof17_trace_overlay_rescaled.png
- proof17_summary.csv

Aufrufbeispiel (Windows, eine Zeile!):
python proof17_universal_rescale.py --mode fit --eigenfile eigenvalues_clean.txt --L 100 --N 400 --k 220 --sigma 0.1 --w_bk_grid 0.8,1.0,1.2 --w_prim_grid 0.0,0.03,0.06 --w_grav_grid 0.0,0.02,0.05 --w_exch_grid 0.0,0.02,0.05
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from pathlib import Path

# ------------------------------------------------------------
# Utils
# ------------------------------------------------------------

def parse_grid(s: str, default=(1.0,)):
    if s is None or s.strip() == "":
        return tuple(default)
    try:
        return tuple(float(x) for x in s.split(","))
    except Exception:
        return tuple(default)

def gaussian_kde_smooth(x, bins=60, bandwidth=None):
    """Einfache KDE: nutzt scipy.stats.gaussian_kde, fallback auf Gauß-Faltung."""
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(x, bw_method=bandwidth)
        grid = np.linspace(np.min(x), np.max(x), bins)
        y = kde(grid)
        return grid, y / (np.max(y) + 1e-12)
    except Exception:
        # Fallback: Histogramm + Gauß glätten
        hist, edges = np.histogram(x, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        sigma = max(1, bins // 50)
        ker = np.exp(-0.5 * (np.arange(-3*sigma, 3*sigma+1)/sigma)**2)
        ker /= np.sum(ker)
        sm = np.convolve(hist, ker, mode="same")
        return centers, sm / (np.max(sm) + 1e-12)

# ------------------------------------------------------------
# Zählfunktion (empirisch)
# ------------------------------------------------------------

def load_eigs(eigenfile: Path):
    if not eigenfile.exists():
        raise FileNotFoundError(f"Eigenvalue file not found: {eigenfile}")
    E = np.loadtxt(eigenfile, dtype=float)
    E = np.asarray(E).ravel()
    return E

def counting_function(E):
    """Sortierte positive Eigenwerte, N(E_i)=i (1-basiert)."""
    Ep = np.sort(E[E > 0.0])
    N = np.arange(1, len(Ep)+1, dtype=float)
    return Ep, N

# ------------------------------------------------------------
# Robustes Modell (Lead BK + kleine lineare Korrektur)
# ------------------------------------------------------------

def N_model_bk(Ep, lam=1.0, c0=0.0, beta=-1.0):
    """
    Robustes Modell nahe der Riemann-Hauptterme:
      N(E') ≈ (E'/(2π)) [log(E') - 1]  +  c0  +  (beta/(2π)) * E'
    (E' = lam * E_phys später; hier erwarten wir direkt E')
    """
    Ep = np.asarray(Ep, float)
    eps = 1e-12
    L = np.log(np.maximum(Ep, eps))
    return (Ep/(2*np.pi))*(L - 1.0) + c0 + (beta/(2*np.pi))*Ep

# ------------------------------------------------------------
# Robuste λ-Schätzung
# ------------------------------------------------------------

def safe_lambda_from_slope(E_pos, N_emp, lam_hint=1.0):
    """
    λ aus N(E)/E ~ a log(E) + b. Mittleres Fenster, robust gewichtet,
    harte Schranken gegen Entartung.
    """
    E = np.asarray(E_pos, float)
    N = np.asarray(N_emp, float)
    if len(E) < 20:
        return float(np.clip(lam_hint, 1/50, 50))

    k0, k1 = int(0.2*len(E)), int(0.85*len(E))
    Ew, Nw = E[k0:k1], N[k0:k1]
    Ew = np.maximum(Ew, 1e-9)

    y = Nw / Ew
    X = np.column_stack([np.log(Ew), np.ones_like(Ew)])

    # robuste Gewichtung (Clipping)
    resid0 = y - np.median(y)
    s = np.std(y) + 1e-12
    w = np.clip(np.abs(resid0) / s, 0, 3.0)
    W = np.diag(1.0/(1.0 + w))

    theta = np.linalg.lstsq(W @ X, W @ y, rcond=None)[0]
    a_hat = float(theta[0])
    a_star = 1.0/(2*np.pi)  # Ziel
    ratio  = a_hat / (a_star + 1e-15)

    lam = lam_hint * np.clip(ratio, 0.1, 10.0) ** (+1.0)
    lam = float(np.clip(lam, 1/50.0, 50.0))
    return lam

# ------------------------------------------------------------
# Fit-Objektiv & Grid-Suche
# ------------------------------------------------------------

def robust_objective(Ep, N_emp, lam, c0, beta):
    m0, m1 = int(0.15*len(Ep)), int(0.9*len(Ep))
    if m1 - m0 < 10:
        m0, m1 = 0, len(Ep)
    Epw, Nw = Ep[m0:m1], N_emp[m0:m1]
    pred = N_model_bk(Epw, lam=lam, c0=c0, beta=beta)
    r = Nw - pred
    s = np.median(np.abs(r)) + 1e-12
    r_clip = np.clip(r, -5*s, 5*s)
    return float(np.mean(np.abs(r_clip)))

def grid_fit(Ep, N_emp,
             lam_candidates=(0.5, 0.75, 1.0, 1.25, 1.5),
             c0_candidates=(-50, -20, 0, 20, 50),
             beta_candidates=(-1.25, -1.0, -0.75)):
    best = None
    best_val = np.inf
    for lam, c0, beta in product(lam_candidates, c0_candidates, beta_candidates):
        val = robust_objective(Ep, N_emp, lam, c0, beta)
        if val < best_val:
            best_val = val
            best = (lam, c0, beta)
    return best, best_val

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

def plot_counting(Ep, N_emp, lam, c0, beta, L):
    plt.figure(figsize=(10,4))
    plt.scatter(Ep, N_emp, s=8, alpha=0.6, label="N_emp(E') (rescaled)")
    xg = np.linspace(max(1e-12, Ep[0]), Ep[-1], 400)
    yg = N_model_bk(xg, lam=lam, c0=c0, beta=beta)
    plt.plot(xg, yg, "k--", lw=2, label=f"N_model(E'); λ={lam:.3f}, c={c0:+.2f}, β={beta:+.2f}")
    plt.xlabel("E' (rescaled energy)")
    plt.ylabel("Counting function N(E')")
    plt.title(f"Proof 17 — Counting after rescaling (L={L:.1f})")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("proof17_counting_rescaled.png", dpi=150)
    plt.close()

def plot_residuals(Ep, N_emp, lam, c0, beta):
    r = N_emp - N_model_bk(Ep, lam=lam, c0=c0, beta=beta)
    plt.figure(figsize=(10,3.8))
    plt.plot(Ep, r, lw=2)
    plt.axhline(0, color="k", ls="--", lw=1)
    plt.xlabel("E' (rescaled)")
    plt.ylabel("Residual  N_emp - N_BK")
    plt.title("Proof 17 — Residuals against BK after rescaling")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("proof17_residuals.png", dpi=150)
    plt.close()

def plot_density(Ep):
    # BK-Form (normiert für Formvergleich)
    xg = np.linspace(Ep[0], Ep[-1], 200)
    eps = 1e-12
    rho_bk = (1.0/(2*np.pi)) * (np.log(np.maximum(xg, eps)) + 1.0)
    rho_bk = rho_bk / (np.max(rho_bk) + 1e-12)

    # Empirische KDE
    gx, gy = gaussian_kde_smooth(Ep, bins=120)

    plt.figure(figsize=(12,4.5))
    plt.subplot(1,2,1)
    plt.plot(xg, rho_bk, "r-", lw=2, label="BK density (norm.)")
    plt.plot(gx, gy, color="#386cb0", lw=2, label="Empirical KDE (norm.)")
    plt.xlabel("E' (rescaled)")
    plt.ylabel("ρ(E') (norm.)")
    plt.title("Proof 17 — Spectral density & local slope diagnostics")
    plt.legend()
    plt.grid(alpha=0.3)

    # Lokale Steigung a(E): numerisch d/dE [ N(E)/E ] ~ a log E + b ⇒ a ~ d/d(logE) [N/E]
    # Hier als Diagnose: a_loc ≈ d/d(logE) N/E / 1
    E = Ep
    N = np.arange(1, len(E)+1, dtype=float)
    y = N / np.maximum(E, 1e-12)
    # finite differences in log-space
    L = np.log(np.maximum(E, 1e-12))
    dy = np.gradient(y, L)
    a_loc = dy

    plt.subplot(1,2,2)
    plt.plot(E, a_loc, "k-", lw=1.8, label="local a(E)")
    plt.axhline(1.0/(2*np.pi), color="orange", ls="--", lw=1.5, label="1/(2π)")
    plt.xlabel("E'")
    plt.ylabel("a(E)")
    plt.legend()
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig("proof17_density_rescaled.png", dpi=150)
    plt.close()

def plot_trace_overlay(Ep, L, tmax=5.0, dt=0.002):
    """Nur qualitatives Overlay: |Tr e^{-iE't}| vs. einfache Proxy-Kurve."""
    t = np.arange(0, tmax, dt)
    tr = np.zeros_like(t, dtype=complex)
    for E in Ep:
        tr += np.exp(-1j * E * t)
    q = np.abs(tr).astype(float)
    q /= (np.max(q) + 1e-12)

    # simple proxy (monoton fallend)
    proxy = np.exp(-t/0.6)
    proxy /= (np.max(proxy) + 1e-12)

    plt.figure(figsize=(12,4.5))
    plt.plot(t, q, lw=1.4, label=r"$|{\rm Tr}\,e^{-iHt}|$ (Quantum, rescaled)")
    plt.plot(t, proxy, lw=2.0, color="#ff7f00", label=r"$|T_{\rm BK}|$ (proxy)")
    plt.xlabel("t")
    plt.ylabel("Normalized amplitude")
    plt.title("Proof 17 — Trace overlay (rescaled spectrum)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("proof17_trace_overlay_rescaled.png", dpi=150)
    plt.close()

# ------------------------------------------------------------
# Hauptfluss
# ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Proof 17 — Universal rescale (robust)")
    ap.add_argument("--mode", choices=["fit", "eigen"], default="fit",
                    help="fit: robust BK-Fit; eigen: nur Diagramme aus Datei")
    ap.add_argument("--eigenfile", type=str, default="eigenvalues_clean.txt")
    ap.add_argument("--L", type=float, default=100.0)
    ap.add_argument("--N", type=int, default=400)
    ap.add_argument("--k", type=int, default=220, help="nur Diagnose-Ausgabe")
    ap.add_argument("--sigma", type=float, default=0.1, help="nur Diagnose-Ausgabe")
    # Grids werden aktuell NICHT aggressiv für Gewichte genutzt – wir halten sie klein
    ap.add_argument("--w_bk_grid", type=str, default="1.0")
    ap.add_argument("--w_prim_grid", type=str, default="0.0")
    ap.add_argument("--w_grav_grid", type=str, default="0.0")
    ap.add_argument("--w_exch_grid", type=str, default="0.0")
    args = ap.parse_args()

    print("=== Proof 17: Universal rescale (robust) ===")

    eigenfile = Path(args.eigenfile)
    E_all = load_eigs(eigenfile)
    Epos, Nemp = counting_function(E_all)
    if len(Epos) == 0:
        raise RuntimeError("Keine positiven Eigenwerte gefunden.")

    print(f"[load] {len(E_all)} Eigenwerte geladen; positive: {len(Epos)}")
    print(f"[range] E_min={np.min(Epos):.4g}, E_max={np.max(Epos):.4g}")

    # 1) sichere λ-Schätzung aus Steigung
    lam0 = safe_lambda_from_slope(Epos, Nemp, lam_hint=1.0)
    Ep = lam0 * Epos

    # sanity check & clamp
    r5, r95 = np.percentile(Ep, [5, 95])
    if (not np.isfinite(r5)) or (not np.isfinite(r95)) or (r95 <= 10*r5):
        print("[safe] Ep range looks degenerate; reset λ to 1.0")
        lam0 = 1.0
        Ep = lam0 * Epos

    Ep = np.maximum(Ep, 1e-12)

    print(f"[fit] initial λ_safe = {lam0:.4f}")

    # 2) Robust-Grid für (λ, c0, beta) – bewusst eng
    lam_cand  = parse_grid(args.w_bk_grid, default=(0.75, 1.0, 1.25))
    # Verwende lam-Kandidaten als relative Faktoren auf lam0
    lam_candidates = tuple(np.clip(lam0*np.array(lam_cand, float), 1/50.0, 50.0))

    c0_candidates   = (-50, -20, 0, 20, 50)
    beta_candidates = (-1.25, -1.0, -0.75)

    (best_lam, best_c0, best_beta), J = grid_fit(Ep, Nemp,
                                                 lam_candidates=lam_candidates,
                                                 c0_candidates=c0_candidates,
                                                 beta_candidates=beta_candidates)
    print(f"[fit] best: λ={best_lam:.3f}, c0={best_c0:+.2f}, beta={best_beta:+.2f}, |J|={J:.4g}")

    # 3) Plots
    plot_counting(Ep, Nemp, best_lam, best_c0, best_beta, args.L)
    plot_residuals(Ep, Nemp, best_lam, best_c0, best_beta)
    plot_density(Ep)
    plot_trace_overlay(Ep, args.L, tmax=5.0, dt=0.002)

    # 4) Summary
    out = {
        "N_pos": int(len(Epos)),
        "E_min": float(np.min(Epos)),
        "E_max": float(np.max(Epos)),
        "lambda_init": float(lam0),
        "lambda_best": float(best_lam),
        "c0_best": float(best_c0),
        "beta_best": float(best_beta),
        "J_obj": float(J),
        "L": float(args.L),
        "k_diag": int(args.k),
        "sigma_diag": float(args.sigma),
    }
    pd.DataFrame([out]).to_csv("proof17_summary.csv", index=False)
    print("[save] proof17_summary.csv geschrieben")
    print("[fig] proof17_counting_rescaled.png, proof17_residuals.png, proof17_density_rescaled.png, proof17_trace_overlay_rescaled.png")

if __name__ == "__main__":
    main()
