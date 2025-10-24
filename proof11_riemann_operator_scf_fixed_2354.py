#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof11_riemann_operator_scf_fixed_2334.py

Dirichlet-Innenraum-Version des selbstkonsistenten Riemann-Operators H[ρ]
- Diskretisierung auf [1, L] mit N Gitterpunkten
- Dirichlet-Randbedingungen: ψ(1)=ψ(L)=0 ⇒ Innenraum N_in=N-2
- Tridiagonaler Laplace-Operator, Diagonal-Potential V_eff nur im Innenraum
- SCF-Iteration mit Dichte ρ aus Boltzmann-gemittelten Eigenzuständen
- Energie-Shift in Boltzmann-Gewichten für numerische Stabilität
- Konsistente Normierung ∫ ρ dx = 1 (numerisch)
- Konvergenzkriterium: (||Δρ||₂ < eps) ODER (|ΔF| < eps_FE)
- Eigenwerte werden gespeichert (eigenvalues_H_rho_final.txt)
- Bridge-Fallback: versucht proof10_equivalence_rigorous_bridge,
  ansonsten _1513, dann _1414; falls nichts passt → interne (trunkierte) Spektral-ζ
- Export: ζ_H-Tabelle (zetaH_table.csv), relative Symmetrie-Heuristik (skaleninvariant)
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional
from pathlib import Path
import importlib
import csv

try:
    # Stabil und schnell für reelle sym. Tridiagonale
    from scipy.linalg import eigh_tridiagonal
except Exception as e:
    raise RuntimeError("Benötigt scipy.linalg.eigh_tridiagonal.") from e


# =========================
# Konfig
# =========================
@dataclass
class SCFConfig:
    L: float = 8.0          # Intervall [1, L]
    N: int = 800             # Gesamtgitter inkl. Ränder
    beta: float = 0.30       # Mixing-Faktor für ρ (0..1)
    eps: float = 1e-3        # Konvergenzschwelle für ||Δρ||
    eps_FE: float = 1e-5     # Konvergenzschwelle für ΔF
    max_iter: int = 200
    n_eigs: int = 400        # Anzahl Eigenpaare für Dichte (≤ N_in)
    epsilon_floor: float = 1e-12
    save_dir: str = "./"
    seed: int = 42
    T: float = 0.75          # „Temperatur“ für Boltzmann-Gewichte


# =========================
# Operator-Klasse
# =========================
class RiemannOperatorSCF:
    def __init__(self, cfg: SCFConfig):
        self.cfg = cfg
        rng = np.random.default_rng(cfg.seed)

        # Gitter & Innenraum
        self.x = np.linspace(1.0, cfg.L, cfg.N)
        self.dx = (cfg.L - 1.0) / (cfg.N - 1)
        self.idx_inner = slice(1, cfg.N - 1)
        self.x_inner = self.x[self.idx_inner]
        self.N_in = cfg.N - 2

        # Initiale Dichte im Innenraum: leicht verrauschte Uniforme
        rho0 = np.ones(self.N_in)
        rho0 += 0.01 * rng.standard_normal(self.N_in)
        rho0 = np.clip(rho0, 0.0, None)
        rho0 /= (np.sum(rho0) * self.dx)
        self.rho = rho0

        self.evals: Optional[np.ndarray] = None
        self.evecs: Optional[np.ndarray] = None

    # --------- Modell-Potentiale (glatte Platzhalter) ---------
    def V_prim(self, x: np.ndarray) -> np.ndarray:
        # Glatte, schwache Modulation – genügt als Testbett
        return 0.02 * np.log(np.log(np.maximum(x, 1.000001)))

    def V_grav(self, rho: np.ndarray) -> np.ndarray:
        # Einfacher lokaler „Hartree“-Term
        kappa = 0.5
        return kappa * rho

    def V_exch(self, rho: np.ndarray) -> np.ndarray:
        # Schwacher Austausch-Term
        gamma = 0.3
        return -gamma * np.sqrt(np.maximum(rho, 0.0) + 1e-16)

    def V_eff(self, rho: np.ndarray) -> np.ndarray:
        return self.V_prim(self.x_inner) + self.V_grav(rho) + self.V_exch(rho)

    # --------- Diskretisierung: -d²/dx² (Dirichlet im Innenraum) ---------
    def _laplacian_tridiagonal(self) -> Tuple[np.ndarray, np.ndarray]:
        inv_dx2 = 1.0 / (self.dx ** 2)
        diag = 2.0 * inv_dx2 * np.ones(self.N_in)
        off = -1.0 * inv_dx2 * np.ones(self.N_in - 1)
        return diag, off

    def _hamiltonian_tridiagonal(self, V_eff: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        diag, off = self._laplacian_tridiagonal()
        diag = diag + V_eff  # H = -d²/dx² + V
        return diag, off

    # --------- Eigenproblem ---------
    def solve_eigen(self, V_eff: np.ndarray, n_eigs: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        if n_eigs is None:
            n_eigs = self.cfg.n_eigs
        n_eigs = int(min(n_eigs, self.N_in))
        diag, off = self._hamiltonian_tridiagonal(V_eff)
        w, Z = eigh_tridiagonal(diag, off, select='i', select_range=(0, n_eigs - 1))
        # Normierung der Eigenvektoren bzgl. ∑ |ψ|² dx = 1
        norm = np.sqrt(np.sum(Z * Z, axis=0) * self.dx)
        Z /= norm
        return w, Z

    # --------- Dichte aus Spektrum (mit Energie-Shift) ---------
    def density_from_spectrum(self, evals: np.ndarray, evecs: np.ndarray) -> np.ndarray:
        betaT = 1.0 / max(self.cfg.T, 1e-12)
        E0 = float(np.min(evals))
        weights = np.exp(-betaT * (evals - E0))
        weights /= np.sum(weights)
        rho_new = np.sum((evecs * evecs) * weights[None, :], axis=1)
        rho_new = np.maximum(rho_new, self.cfg.epsilon_floor)
        rho_new /= (np.sum(rho_new) * self.dx)
        return rho_new

    def free_energy(self, evals: np.ndarray) -> float:
        betaT = 1.0 / max(self.cfg.T, 1e-12)
        E0 = float(np.min(evals))
        Z = np.sum(np.exp(-betaT * (evals - E0)))
        F = E0 - (1.0 / betaT) * np.log(Z)
        return float(F)

    # --------- SCF ---------
    def run_scf(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        rho = self.rho.copy()
        F_prev: Optional[float] = None

        for k in range(1, self.cfg.max_iter + 1):
            V = self.V_eff(rho)
            evals, evecs = self.solve_eigen(V, self.cfg.n_eigs)
            rho_new = self.density_from_spectrum(evals, evecs)

            # Linear Mixing
            rho_mixed = (1.0 - self.cfg.beta) * rho + self.cfg.beta * rho_new
            rho_mixed = np.maximum(rho_mixed, self.cfg.epsilon_floor)
            rho_mixed /= (np.sum(rho_mixed) * self.dx)

            delta = np.linalg.norm(rho_mixed - rho)
            F = self.free_energy(evals)
            dF = np.inf if F_prev is None else abs(F - F_prev)

            print(f"[Iter {k:03d}] ||Δρ||₂ = {delta:.4e} | ΔF = {dF:.4e} | Emin = {evals[0]:.6f}")

            rho = rho_mixed
            F_prev = F

            # ODER statt UND → schnelleres, robusteres Stopkriterium
            if (delta < self.cfg.eps) or (dF < self.cfg.eps_FE):
                print("Konvergenz erreicht.")
                break
        else:
            print("Warnung: Max. Iterationen erreicht (evtl. nicht voll konvergiert).")

        V_final = self.V_eff(rho)
        evals_full, evecs_full = self.solve_eigen(V_final, min(self.cfg.n_eigs, self.N_in))

        self.rho = rho
        self.evals = evals_full
        self.evecs = evecs_full
        return rho, evals_full, evecs_full

    def a1_seeley_dewitt(self) -> float:
        # a1 = (1/√(4π)) * ∫ V_eff(x) dx (heuristisch für dieses Modell)
        V = self.V_eff(self.rho)
        integral = np.sum(V) * self.dx
        return float(integral / np.sqrt(4.0 * np.pi))


# =========================
# Bridge / ζ-Auswertung
# =========================
def import_bridge_module() -> Optional[object]:
    """
    Versucht die Bridge in dieser Reihenfolge zu laden:
    1) proof10_equivalence_rigorous_bridge
    2) proof10_equivalence_rigorous_bridge_1513
    3) proof10_equivalence_rigorous_bridge_1414
    """
    names = [
        "proof10_equivalence_rigorous_bridge",
        "proof10_equivalence_rigorous_bridge_1513",
        "proof10_equivalence_rigorous_bridge_1414",
    ]
    for name in names:
        try:
            return importlib.import_module(name)
        except Exception:
            continue
    return None


def try_bridge_and_export(evals: np.ndarray, save_dir: Path,
                          L: float, a1_value: float,
                          s_min=0.6, s_max=3.0, s_step=0.1) -> Tuple[Path, dict]:
    """
    Versucht die Bridge zu nutzen; ansonsten Fallback: trunkierte ζ_H(s)=∑ λ^{-s} (λ>0).
    Exportiert zetaH_table.csv und liefert Summary (Bridge genutzt?, Symmetrie-Heuristik).
    """
    save_dir.mkdir(parents=True, exist_ok=True)
    table_path = save_dir / "zetaH_table.csv"

    # --- s-Gitter: feines Fenster um 1/2 + grobes Raster 0.6..3.0 ---
    # Ziel:
    #  • Viele Punkte in (0.46, 0.56) (ohne s=0.5) für den Residuen-Fit
    #  • Genug Punkte in [0.6, 3.0] für die globale Kurve
    #  • Keine exakt problematischen Punkte (s=0.5, optional s=1.0)
    s_coarse = np.linspace(0.60, 3.00, 28)           # globales Raster
    s_fine   = np.linspace(0.46, 0.56, 21)           # dichtes Fenster
    s_fine   = s_fine[np.abs(s_fine - 0.50) > 1e-10] # s=0.5 entfernen
    s_coarse = s_coarse[np.abs(s_coarse - 1.00) > 1e-12]  # optional s=1 meiden

    s_values = np.unique(np.concatenate([s_coarse, s_fine]))
    s_values.sort()

    n_window = np.count_nonzero((s_values > 0.46) & (s_values < 0.56))
    print(f"[proof11] s-grid: {len(s_values)} Punkte, davon {n_window} im Fenster (0.46, 0.56) um 1/2.")

    # --- Summary & Container vorbereiten ---
    summary = {"used_bridge": False, "symmetry_error": None, "residue_s_half": None}
    zeta_vals = None

    # --- Bridge laden und versuchen ---
    bridge = import_bridge_module()
    if bridge is not None:
        try:
            if hasattr(bridge, "compute_zeta_from_eigs"):
                try:
                    # neue Signatur (mit L, a1)
                    zeta_vals = bridge.compute_zeta_from_eigs(evals, s_values, L=L, a1_value=a1_value)
                except TypeError:
                    # ältere Signatur (ohne L, a1)
                    zeta_vals = bridge.compute_zeta_from_eigs(evals, s_values)
            elif hasattr(bridge, "main"):
                zeta_vals = bridge.main(evals=evals, s_values=s_values)

            if zeta_vals is not None:
                summary["used_bridge"] = True
                if hasattr(bridge, "symmetry_error"):
                    summary["symmetry_error"] = getattr(bridge, "symmetry_error")
                if hasattr(bridge, "residue_s_half"):
                    summary["residue_s_half"] = getattr(bridge, "residue_s_half")
        except Exception:
            zeta_vals = None

    # --- Fallback: truncierte Spektral-ζ aus Eigenwerten ---
    if zeta_vals is None:
        print("[Info] Bridge nicht verfügbar/inkompatibel – verwende interne ζ-Schätzung.")
        evals_pos = evals[evals > 1e-14]
        zeta_vals = np.array([np.sum(evals_pos ** (-s)) for s in s_values], dtype=float)

        # relative, skaleninvariante Symmetrie-Heuristik (nur Diagnose)
        s_probe = np.array([0.75, 0.85, 0.90])
        if (s_values.min() <= s_probe.min()) and (s_values.max() >= (1.0 - s_probe).max()):
            z_left  = np.interp(s_probe,        s_values, zeta_vals)
            z_right = np.interp(1.0 - s_probe,  s_values, zeta_vals)
            denom = np.maximum(np.mean(np.abs(z_left)), 1e-12)
            summary["symmetry_error"] = float(np.mean(np.abs(z_right - z_left)) / denom)

    # --- CSV-Export ---
    with open(table_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["s", "zeta_H(s)"])
        for s, z in zip(s_values, zeta_vals):
            writer.writerow([f"{s:.3f}", f"{float(z):.12e}"])

    return table_path, summary


# =========================
# Main
# =========================
def main():
    cfg = SCFConfig()
    op = RiemannOperatorSCF(cfg)

    # SCF
    rho, evals, evecs = op.run_scf()

    save_dir = Path(cfg.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # Eigenwerte speichern
    eig_path = save_dir / "eigenvalues_H_rho_final.txt"
    np.savetxt(eig_path, evals)
    print(f"Eigenwerte gespeichert: {eig_path}")

    # a1 (Seeley–DeWitt) ausgeben
    a1 = op.a1_seeley_dewitt()
    print(f"a1 (Seeley–DeWitt) = {a1:.6e}")

    # Bridge / ζ_H – mit L und a1_value
    table_path, summary = try_bridge_and_export(evals, save_dir, L=cfg.L, a1_value=a1)
    print(f"ζ_H-Tabelle exportiert: {table_path}")

    # Zusammenfassung
    print("\n===== ZUSAMMENFASSUNG (proof11_fixed_v2) =====")
    print(f"N_in = {op.N_in}, dx = {op.dx:.4e}, T = {cfg.T}, beta (mix) = {cfg.beta}")
    print(f"Min/Max Eigenwert: {float(np.min(evals)):.6f} / {float(np.max(evals)):.6f}")
    print(f"∫ρ dx (numerisch) = {float(np.sum(rho)*op.dx):.6f}")
    print(f"a1 = {a1:.6e}")
    print(f"Bridge verwendet: {summary['used_bridge']}")
    print(f"Symmetrie-Fehler (heur.): {summary['symmetry_error']}")
    print(f"Residuum bei s=1/2 (falls von Bridge geliefert): {summary['residue_s_half']}")


if __name__ == "__main__":
    main()
