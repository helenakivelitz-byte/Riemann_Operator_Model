# src/asymptotic_rescaling.py
import os, json
import numpy as np
import mpmath as mp

DEFAULTS = {
    "lambda_guess": 1.0 / (2.0 * np.pi),
    "newton_tol": 1e-10,
    "newton_maxit": 50
}

def _load_eigenvalues():
    """
    Robust loader:
    - Prefer data/eigenvalues_raw.json (list OR dict-with-key)
    - Fallback to data/results.json (expects dict with 'eigenvalues')
    """
    import json, os
    paths = [
        "data/eigenvalues_raw.json",
        "data/results.json",
    ]
    for p in paths:
        if not os.path.exists(p):
            continue
        with open(p, "r", encoding="utf-8") as f:
            obj = json.load(f)
        # If it's already a list of numbers
        if isinstance(obj, list):
            return np.array(obj, dtype=float)
        # If it's a dict, check for common keys
        if isinstance(obj, dict):
            for key in ("eigenvalues", "evals", "E", "E_raw"):
                if key in obj:
                    return np.array(obj[key], dtype=float)
    raise FileNotFoundError(
        "Could not load eigenvalues. Expected data/eigenvalues_raw.json (list) "
        "or data/results.json with key 'eigenvalues'."
    )

def _to_float(x):
    return float(mp.re(x)) if hasattr(x, "real") else float(x)

def _to_complex_dict(z):
    return {"re": float(mp.re(z)), "im": float(mp.im(z))}

def N_BK(Ep, lam):
    """Berry–Keating counting: N_BK(E') = λ E' (log E' − 1),  E'>0."""
    return lam * Ep * (mp.log(Ep) - 1)

def inv_N_BK(Ntarget, lam, x0=None, tol=1e-10, maxit=50):
    """Invert N_BK via Newton; Ntarget>0, lam>0."""
    if Ntarget <= 0:
        return mp.mpf("1.0")
    if x0 is None:
        x0 = max(2.0, float(Ntarget / lam + 1.0))
    Ep = mp.mpf(x0)
    for _ in range(maxit):
        f  = lam * Ep * (mp.log(Ep) - 1) - Ntarget
        df = lam * (mp.log(Ep))              # d/dEp [λ Ep (log Ep − 1)] = λ log Ep
        step = f / df
        Ep -= step
        if abs(step) < tol * (1 + abs(Ep)):
            break
    return float(Ep)

class AsymptoticRescaling:
    def __init__(self, eigenvalues, lambda_param=None):
        self.evals_raw = np.sort(np.asarray(eigenvalues, dtype=float))
        self.lam = float(DEFAULTS["lambda_guess"] if lambda_param is None else lambda_param)
        self.evals_rescaled = None
        self.stats = {}

    def N_emp(self, E):
        """Empirical counting: #{E_n ≤ E}."""
        return int(np.searchsorted(self.evals_raw, E, side="right"))

    def perform(self):
        """E'_n = N_BK^{-1}(N_emp(E_n))."""
        print("🔄 Performing Optimal-Transport (BK) rescaling …")
        k = np.arange(1, len(self.evals_raw) + 1, dtype=int)  # N_emp(E_n)=n for sorted spectrum
        Ep = [inv_N_BK(mp.mpf(int(ki)), self.lam,
                       x0=float(ki / self.lam + 1.0),
                       tol=DEFAULTS["newton_tol"],
                       maxit=DEFAULTS["newton_maxit"])
              for ki in k]
        self.evals_rescaled = np.array(Ep, dtype=float)
        self._compute_stats()
        return self.evals_rescaled

    def _compute_stats(self):
        s = self.evals_rescaled / np.maximum(self.evals_raw, 1e-12)
        self.stats = {
            "lambda_parameter": float(self.lam),
            "n_eigenvalues": int(len(self.evals_raw)),
            "mean_scaling_factor": float(np.mean(s)),
            "std_scaling_factor": float(np.std(s)),
            "raw_eigenvalues_range": [float(self.evals_raw.min()), float(self.evals_raw.max())],
            "rescaled_eigenvalues_range": [float(self.evals_rescaled.min()), float(self.evals_rescaled.max())],
            "newton": {
                "tol": DEFAULTS["newton_tol"],
                "maxit": DEFAULTS["newton_maxit"]
            }
        }

    def zeta_rescaled(self, s, K=None):
        """ζ_H_rescaled(s) with first K levels (default: all)."""
        if self.evals_rescaled is None:
            raise RuntimeError("Call perform() first.")
        E = self.evals_rescaled if K is None else self.evals_rescaled[:int(K)]
        Epos = [mp.mpf(e) for e in E if e > 0]
        s_c = mp.mpc(s)  # supports complex s
        return mp.fsum(e ** (-s_c) for e in Epos)

    def bridge_ratio_rescaled(self, s_values, K=None):
        """Ψ_rescaled(s) = ζ_H_rescaled(s)/ζ(s), only for Re(s)>1; JSON-friendly output."""
        mp.mp.dps = 50
        out = {}
        for s in s_values:
            s_c = s if isinstance(s, complex) else complex(s)
            if s_c.real <= 1.0:
                print(f"⚠️ Skipping s={s}: Re(s) must be > 1.")
                continue
            if abs(s_c - 1.0) < 1e-6:
                print(f"⚠️ Skipping s={s}: too close to pole at 1.")
                continue
            zH = self.zeta_rescaled(s_c, K=K)
            zR = mp.zeta(mp.mpc(s_c))
            if abs(zR) < mp.mpf("1e-30"):
                print(f"⚠️ Skipping s={s}: |ζ(s)| too small.")
                continue
            psi = zH / zR
            out[str(s)] = {
                "truncation_K": int(len(self.evals_rescaled) if K is None else int(K)),
                "zeta_H_rescaled": _to_complex_dict(zH),
                "zeta_R": _to_complex_dict(zR),
                "psi_rescaled": _to_complex_dict(psi),
                "abs_deviation_rescaled": float(abs(psi - 1))
            }
        return out

    # --- add inside class AsymptoticRescaling ------------------------------------
    def run_complete_analysis(self, s_test_points=None):
        """
        1) Perform BK/OT rescaling, 2) compute Ψ_rescaled(s) for given s,
        3) return a JSON-serializable dict.
        """
        import mpmath as mp
        mp.mp.dps = 50

        if s_test_points is None:
#            s_test_points = [1.10, 1.15, 1.20]
            s_test_points = np.arange(1.05, 1.51, 0.01)  # dense real grid
             
        # 1) Perform rescaling
        self.perform()

        # 2) Compute bridge ratios
        bridge = self.bridge_ratio_rescaled(s_test_points)

        # 3) Assemble result (JSON safe)
        result = {
            "lambda_used": float(self.lam),
            "lambda_theoretical": float(1.0 / (2.0 * np.pi)),
            "K": int(len(self.evals_rescaled)),
            "eigenvalues_rescaled": self.evals_rescaled.tolist(),
            "rescaling_parameters": self.stats,
            "zeta_rescaled": bridge,
            "s_values": s_test_points
        }

        # Pretty summary
        print("\n" + "=" * 60)
        print("🎯 RESCALING SUMMARY")
        print("=" * 60)
        print(f"• λ used: {result['lambda_used']:.6f}  (theoretical 1/(2π) ≈ {result['lambda_theoretical']:.6f})")
        print(f"• mean scaling factor: {self.stats['mean_scaling_factor']:.3f}  ± {self.stats['std_scaling_factor']:.3f}")
        for s, dat in bridge.items():
            print(f"  s={s:>6} → |Ψ_rescaled−1| = {dat['abs_deviation_rescaled']:8.3e}")
        print("✅ Rescaling analysis complete.")
        
        return result
    # --- end add -----------------------------------------------------------------

def main():
    """Main function - load eigenvalues and run rescaling"""
    try:
        # 1) Load eigenvalues robustly
        eigenvalues = _load_eigenvalues()
        print(f"📊 Loaded {len(eigenvalues)} eigenvalues")

        # 2) s-grid: reuse from results.json if present, else default
        s_values = [1.10, 1.15, 1.20]
        res0 = os.path.join("data", "results.json")
        if os.path.exists(res0):
            try:
                R = json.load(open(res0, "r", encoding="utf-8"))
                if "zeta_s_values" in R:
                    s_values = [float(x) for x in R["zeta_s_values"]]
            except Exception:
                pass

        print("🚀 STARTING ASYMPTOTIC RESCALING ANALYSIS")
        print("=" * 60)
        print(f"• Loaded {len(eigenvalues)} eigenvalues")
        print(f"• s-grid: {s_values}")

        # 2) Run rescaling
        rescaling = AsymptoticRescaling(eigenvalues)
        results = rescaling.run_complete_analysis(s_test_points=s_values)

        # 3) Save results
        os.makedirs('data', exist_ok=True)
        with open('data/results_rescaled.json', 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2)
        print("💾 Rescaled results saved to: data/results_rescaled.json")

        return results

    except Exception as e:
        print(f"❌ Rescaling failed: {e}")
        raise

if __name__ == "__main__":
    main()
