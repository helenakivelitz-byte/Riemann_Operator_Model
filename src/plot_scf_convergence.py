# src/plot_scf_convergence.py
import json, os
import numpy as np
import matplotlib.pyplot as plt

def load_scf_history(path="data/diagnostics/scf_history.json"):
    if not os.path.exists(path):
        raise FileNotFoundError(f"SCF history not found at {path}. "
                                "Run the main pipeline first so it writes scf_history.json.")
    with open(path, "r", encoding="utf-8") as f:
        hist = json.load(f)
    if not hist:
        raise ValueError("SCF history is empty.")
    iters = np.array([h["iter"] for h in hist], dtype=int)
    norms = np.array([h["norm"] for h in hist], dtype=float)
    return iters, norms, hist

def plot_scf_convergence(iters, norms, out_path="data/figures/scf_convergence.png"):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.figure(figsize=(7.5, 5))
    plt.semilogy(iters, norms, "o-", linewidth=1.5, markersize=4)
    plt.xlabel("SCF iteration")
    plt.ylabel(r"$\|\rho_{k+1}-\rho_k\|_2$")
    plt.title("SCF Convergence")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()
    return out_path

def main():
    iters, norms, hist = load_scf_history()
    out = plot_scf_convergence(iters, norms)
    print(f"✅ SCF convergence plot saved to: {out}")

if __name__ == "__main__":
    main()
