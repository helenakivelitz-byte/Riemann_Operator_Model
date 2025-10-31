# src/main_pipeline.py
import os
import json
import numpy as np
from .operator import OperatorBuilder
from .scf_solver import SCFSolver
from .spectral_analyzer import SpectralAnalyzer
from .zeta_comparison import ZetaComparison
from .analytic_checks import AnalyticPrecheck

DEFAULT_CONFIG = {
    "L": 12.0,
    "n_points": 1201,
    "epsilon": 0.1,
    "gamma": 0.35,
    "eta": 0.25,
    "sigma": 0.5,
    "rho_min": 1e-6,
    "N_states": 32,
    "beta": 0.5,
    "tol": 1e-6,
    "max_iter": 200,
    "temperature": None,
    "s_values": [1.1, 1.15, 1.2],
    "out_dir": "data"
}

class RiemannPipeline:
    def __init__(self, cfg):
        self.cfg = {**DEFAULT_CONFIG, **(cfg or {})}
        self.op = OperatorBuilder(
            L=self.cfg["L"],
            n_points=self.cfg["n_points"],
            epsilon=self.cfg["epsilon"],
            gamma=self.cfg["gamma"],
            eta=self.cfg["eta"],
            sigma=self.cfg["sigma"],
            rho_min=self.cfg["rho_min"]
        )
        self.scf_solver = SCFSolver(
            self.op,
            N_states=self.cfg["N_states"],
            beta=self.cfg["beta"],
            tol=self.cfg["tol"],
            max_iter=self.cfg["max_iter"],
            temperature=self.cfg["temperature"]
        )
        os.makedirs(self.cfg["out_dir"], exist_ok=True)

    def run_complete_analysis(self):
        print("🚀 START RIEMANN-OPERATOR PIPELINE")
        print("=" * 60)

        # 0) Analytic pre-check
        pre = AnalyticPrecheck(self.op).run()

        # 1) SCF
        print("\n1) Computing self-consistent density ...")
        rho_star, eigenvalues, eigenvectors, scf_hist = self.scf_solver.solve()

        # 2) Spectral sanity
        print("\n2) Spectral sanity checks ...")
        sa = SpectralAnalyzer(self.op.x, eigenvalues)
        weyl = sa.weyl_sanity()
        stats = sa.spacing_statistics()

        # 3) Zeta comparison (safe domain)
        print("\n3) Zeta bridge ratio (Re(s) > 1) ...")
        s_values = self.cfg.get("s_values", [1.10, 1.15, 1.20])
        zc = ZetaComparison(eigenvalues, s_values)
        psi = zc.bridge_ratio()

        # ensure output directory exists
        out_dir = self.cfg.get("out_dir", "data")
        os.makedirs(out_dir, exist_ok=True)

        # save eigenvalues (first K used in zeta)
        K = self.cfg["N_states"]
        with open(os.path.join(out_dir, "eigenvalues_raw.json"), "w", encoding="utf-8") as f:
            json.dump(list(map(float, eigenvalues[:K])), f, indent=2)

        # prepare and save results bundle
        results = {
            "precheck": pre,
            "spectral": {"weyl": weyl, "gaps": stats},
            "eigenvalues": list(map(float, eigenvalues[:K])),
            "zeta": psi,               # your existing Ψ-raw structure
            "zeta_s_values": self.cfg.get("s_values", [1.10, 1.15, 1.20]),
            "scf_history": scf_hist,             # from scf_solver.solve() return
            "K": int(K),
        }

        with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)

        # save extra copies for diagnostics
        os.makedirs("data/diagnostics", exist_ok=True)
        with open("data/diagnostics/scf_history.json", "w", encoding="utf-8") as f:
            json.dump(scf_hist, f, indent=2)

        print("💾 Eigenvalues saved to", os.path.join(out_dir, "eigenvalues_raw.json"))
        print("✅ Done. Saved to", os.path.join(out_dir, "results.json"))
        
        return results


def load_yaml(path):
    try:
        import yaml
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}

def main(config_path="config/parameters.yaml"):
    print("🎯 Starte Riemann Operator Pipeline...")
    cfg = load_yaml(config_path)
    pipe = RiemannPipeline(cfg)
    return pipe.run_complete_analysis()
