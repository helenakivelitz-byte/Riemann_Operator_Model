# src/compute_diagnostics.py

# Robust intra-package imports: prefer relative; fall back if run directly.
try:
    from .diagnostics_utils import (
        load_results, eigenvalues_from_results, ensure_dirs,
        dense_s_grid, bridge_ratio_absdev, weyl_count_leading,
        quantile_rescale_BK
    )
except ImportError:
    # Fallback for accidental "python src/compute_diagnostics.py" execution
    import os, sys
    sys.path.append(os.path.dirname(__file__))
    from diagnostics_utils import (
        load_results, eigenvalues_from_results, ensure_dirs,
        dense_s_grid, bridge_ratio_absdev, weyl_count_leading,
        quantile_rescale_BK
    )

import os, json, csv
import numpy as np
import matplotlib.pyplot as plt

def compute_scf_convergence():
    path = "data/diagnostics/scf_history.json"
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        hist = json.load(f)
    it = [h["iter"] for h in hist]
    nm = [h["norm"] for h in hist]
    return {"iters": it, "norms": nm}

def compute_bridge_dense(RAW, RES):
    # Use more eigenvalues for better zeta convergence
    E = eigenvalues_from_results(RAW)
    
    # Extend s-grid to 1.05–1.6 for broader diagnostic range
    sgrid = dense_s_grid(RAW, 1.05, 1.60, 0.01)
    
    # Use optimal K truncation (use all available eigenvalues for better accuracy)
    K_optimal = min(len(E), 200)  # Use up to 200 eigenvalues for better convergence
    
    dense_raw = {}
    for s in sgrid:
        val = bridge_ratio_absdev(E, s, K=K_optimal)
        if val is not None:
            dense_raw[str(s)] = val
            
    # rescaled with lambda from res file
    lam = float(RES.get("lambda_used", 1/(2*np.pi)))
    Eprime = quantile_rescale_BK(E, lam)
    
    dense_res = {}
    for s in sgrid:
        val = bridge_ratio_absdev(Eprime, s, K=K_optimal)
        if val is not None:
            dense_res[str(s)] = val
            
    return sgrid, dense_raw, dense_res, lam, K_optimal

def estimate_lambda_with_error(RAW):
    # Use more eigenvalues for better lambda estimation
    E = eigenvalues_from_results(RAW)
    center = 1/(2*np.pi)
    grid = np.linspace(0.95*center, 1.05*center, 31)  # Finer grid: 31 points
    
    # Extended s-grid for more robust optimization
    sgrid = dense_s_grid(RAW, 1.10, 1.50, 0.02)
    
    # Use optimal K truncation
    K_optimal = min(len(E), 200)
    
    scores = []
    for lam in grid:
        Eprime = quantile_rescale_BK(E, lam)
        # objective: median(|Psi-1|) across extended s-grid (more robust)
        vals = [bridge_ratio_absdev(Eprime, s, K=K_optimal) for s in sgrid]
        vals = [v for v in vals if v is not None]
        scores.append(np.median(vals) if vals else np.inf)
        
    scores = np.array(scores)
    j = int(np.argmin(scores))
    lam_hat = grid[j]

    # Improved quadratic fit with more points around minimum
    fit_range = 3  # Use 3 points on each side for better fit
    j_start = max(0, j - fit_range)
    j_end = min(len(grid), j + fit_range + 1)
    
    if j_end - j_start >= 3:  # Need at least 3 points for quadratic fit
        x = grid[j_start:j_end]
        y = scores[j_start:j_end]
        # fit y = a x^2 + b x + c
        A = np.vstack([x*x, x, np.ones_like(x)]).T
        a, b, c = np.linalg.lstsq(A, y, rcond=None)[0]
        
        # Better error estimation using curvature and residual
        if a > 0:
            # Use residual standard error for better error bars
            y_pred = a*x*x + b*x + c
            residuals = y - y_pred
            rmse = np.sqrt(np.mean(residuals**2))
            half_width = rmse / (a * len(x)**0.5) if a > 0 else 0.001
        else:
            half_width = 0.001
    else:
        half_width = 0.001

    return float(lam_hat), float(half_width), list(map(float, grid)), list(map(float, scores)), K_optimal

def compute_weyl(RAW):
    L = RAW.get("precheck", {}).get("domain", {}).get("L", 12.0)
    E = eigenvalues_from_results(RAW)
    if len(E) == 0:
        return {}
    Emax = float(np.max(E))
    N_emp = int(len(E))
    N_weyl = weyl_count_leading(L, Emax)
    return {"L": float(L), "Emax": Emax, "N_emp": N_emp, "N_weyl_leading": float(N_weyl)}

def save_json_csv(obj, json_path, csv_rows=None):
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    if csv_rows:
        csv_path = os.path.splitext(json_path)[0] + ".csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerows(csv_rows)

def plot_scf_convergence(conv):
    if not conv: return
    it = conv["iters"]; nm = conv["norms"]
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7.5,5))
    plt.semilogy(it, nm, "o-")
    plt.xlabel("SCF iteration"); plt.ylabel(r"$\|\rho_{k+1}-\rho_k\|_2$")
    plt.title("SCF Convergence")
    plt.grid(True, ls="--", alpha=0.4)
    os.makedirs("data/figures", exist_ok=True)
    plt.tight_layout()
    plt.savefig("data/figures/scf_convergence.png", dpi=220)
    plt.close()

def plot_bridge_dense(sgrid, dense_raw, dense_res, K_used):
    import matplotlib.pyplot as plt
    x = np.array([float(s) for s in dense_raw.keys()])
    y1 = np.array([dense_raw[str(s)] for s in x])
    y2 = np.array([dense_res[str(s)] for s in x])
    plt.figure(figsize=(8,5))
    plt.plot(x, y1, "o-", label="raw")
    plt.plot(x, y2, "o-", label="rescaled")
    plt.yscale("log")
    plt.xlabel("s"); plt.ylabel(r"$|\Psi(s;L)-1|$")
    plt.title(f"Bridge Ratio (K={K_used}, dense grid, log-y)")
    plt.grid(True, ls="--", alpha=0.35)
    plt.axhline(1e-2, ls=":", color="gray", lw=1, label=r"$10^{-2}$")
    plt.axhline(1e-1, ls=":", color="red", lw=1, label=r"$10^{-1}$")
    plt.legend()
    os.makedirs("data/figures", exist_ok=True)
    plt.tight_layout()
    plt.savefig("data/figures/psi_compare_dense_log.png", dpi=220)
    plt.close()

def plot_lambda_sweep(grid, scores, lam_hat, err, K_used):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7.5,5))
    plt.plot(grid, scores, "o-")
    plt.axvline(lam_hat, color="tab:red", ls="--", label=fr"$\hat\lambda={lam_hat:.6f}\pm{err:.4f}$")
    plt.axvline(1/(2*np.pi), color="tab:green", ls=":", label=r"$1/(2\pi)$")
    plt.xlabel(r"$\lambda$")
    plt.ylabel("median |Psi_rescaled-1| over s-grid")
    plt.title(f"Lambda sweep (K={K_used}, BK rescaling)")
    plt.grid(True, ls="--", alpha=0.35)
    plt.legend()
    os.makedirs("data/figures", exist_ok=True)
    plt.tight_layout()
    plt.savefig("data/figures/lambda_sweep.png", dpi=220)
    plt.close()

def main():
    ensure_dirs()
    RAW, RES = load_results()

    # 1) SCF convergence
    conv = compute_scf_convergence()
    save_json_csv(conv or {}, "data/diagnostics/scf_convergence.json")
    plot_scf_convergence(conv)

    # 2) Dense bridge ratio (raw & rescaled) with extended s-grid
    sgrid, dense_raw, dense_res, lam_used, K_optimal = compute_bridge_dense(RAW, RES)
    dense = {
        "s_grid": list(map(float, sgrid)),
        "absdev_raw": dense_raw,
        "absdev_rescaled": dense_res,
        "lambda_used": float(lam_used),
        "K_used": K_optimal
    }
    # also export as a flat CSV for LaTeX tables
    rows = [["s","absdev_raw","absdev_rescaled"]]
    for s in sgrid:
        sr = str(float(s))
        rows.append([float(s), dense_raw.get(sr,None), dense_res.get(sr,None)])
    save_json_csv(dense, "data/diagnostics/bridge_dense.json", rows)
    plot_bridge_dense(sgrid, dense_raw, dense_res, K_optimal)

    # 3) Weyl count comparison
    w = compute_weyl(RAW)
    rows = [["L","Emax","N_emp","N_weyl_leading"],
            [w.get("L"), w.get("Emax"), w.get("N_emp"), w.get("N_weyl_leading")]]
    save_json_csv(w, "data/diagnostics/weyl_comparison.json", rows)

    # 4) Lambda estimate with error bar and improved optimization
    lam_hat, lam_err, grid, scores, K_lambda = estimate_lambda_with_error(RAW)
    lam_out = {
        "lambda_hat": lam_hat,
        "lambda_err": lam_err,
        "grid": grid,
        "scores_median_absdev": scores,
        "K_used": K_lambda
    }
    rows = [["lambda","median_abs_deviation"]] + [[g,s] for g,s in zip(grid,scores)]
    save_json_csv(lam_out, "data/diagnostics/lambda_fit.json", rows)
    plot_lambda_sweep(grid, scores, lam_hat, lam_err, K_lambda)

    print("✅ Diagnostics complete.")
    print(f"   K used for zeta sums: {K_optimal} eigenvalues")
    print(f"   Lambda estimate: {lam_hat:.6f} ± {lam_err:.4f}")
    print("   data/diagnostics/: scf_convergence.json, bridge_dense.json, weyl_comparison.json, lambda_fit.json")
    print("   data/figures/: scf_convergence.png, psi_compare_dense_log.png, lambda_sweep.png")

if __name__ == "__main__":
    main()
