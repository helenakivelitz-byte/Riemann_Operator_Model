#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof13_limitcheck.py
Konvergenztest L→∞:
  Ψ_L(s) = ζ_H,L(s) / ζ(s),  Ψ_{L,rel}(s) := Ψ_L(s) / Ψ_L(s_ref)

Unterstützte Modi:
  (A) --mode csv  : lies zetaH-Tabellen von Disk (ein File pro L)
      -> --Ls 8 10 12 14  --csv-pattern "zetaH_L{L}.csv"
  (B) --mode evals: erzeuge zetaH-Tabellen aus Eigenwertlisten via Bridge
      -> --Ls 8 10 12 14  --evals-pattern "eigs_L{L}.txt"

Ausgaben (mit --out-prefix proof13):
  - proof13_curves.png            (Ψ_rel-Kurven)
  - proof13_dev_curves.png        (|Ψ_rel-1|)
  - proof13_sup_vs_invL.png       (sup_s |Ψ_rel-1| vs. 1/L + Linearfit)
  - proof13_summary.csv/.txt/.md  (Kennzahlen je L)
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Optional


# ---------- robuste CSV-Leser für zetaH ----------
def load_zetaH_table(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path)
    # s-Spalte
    s_col = None
    for c in df.columns:
        if c.strip().lower() == "s":
            s_col = c; break
    if s_col is None:
        s_col = df.columns[0]
    # zeta-Spalte
    z_cols = [c for c in df.columns if "zeta" in c.lower()]
    z_col = z_cols[0] if z_cols else df.columns[1]
    s = df[s_col].to_numpy(float)
    z = df[z_col].to_numpy(float)
    return s, z


# ---------- Bridge (falls wir aus Eigenwerten rechnen) ----------
def try_import_bridge():
    for name in [
        "proof10_equivalence_rigorous_bridge",
        "proof10_equivalence_rigorous_bridge_1513",
        "proof10_equivalence_rigorous_bridge_1414",
    ]:
        try:
            return __import__(name)
        except Exception:
            continue
    return None


def compute_zetaH_from_evals(evals: np.ndarray, s: np.ndarray, L: float, a1: Optional[float] = None) -> np.ndarray:
    bridge = try_import_bridge()
    if bridge is None:
        # Fallback: reine Dirichletserie (NICHT für Residuen-Fit, aber hier ok)
        epos = evals[evals > 1e-14]
        return np.array([np.sum(epos**(-si)) for si in s], dtype=float)
    # bevorzugte neue Signatur
    try:
        return np.asarray(bridge.compute_zeta_from_eigs(evals, s, L=float(L), a1_value=(None if a1 is None else float(a1))), float)
    except TypeError:
        return np.asarray(bridge.compute_zeta_from_eigs(evals, s), float)


# ---------- Ψ-Bildung ----------
def safe_zeta_real(s: np.ndarray) -> np.ndarray:
    out = np.empty_like(s, dtype=float)
    for i, si in enumerate(s):
        if abs(si - 1.0) < 1e-12:
            out[i] = np.nan  # Pol
        else:
            out[i] = float(mp.zeta(si))
    return out


def build_psi_rel(s: np.ndarray, zetaH: np.ndarray, s_ref: float) -> Tuple[np.ndarray, np.ndarray, float]:
    zeta = safe_zeta_real(s)
    mask = np.isfinite(zeta) & np.isfinite(zetaH) & (np.abs(zeta) > 0)
    s_use, zH_use, z_use = s[mask], zetaH[mask], zeta[mask]

    # scale bei s_ref
    if not (s_use.min() <= s_ref <= s_use.max()):
        # nimm nächstliegenden Punkt
        s_ref_eff = float(s_use[np.argmin(np.abs(s_use - s_ref))])
    else:
        s_ref_eff = float(s_ref)
    zH_ref = float(np.interp(s_ref_eff, s_use, zH_use))
    z_ref  = float(np.interp(s_ref_eff, s_use, z_use))
    scale  = zH_ref / z_ref

    psi     = zH_use / z_use
    psi_rel = psi / scale
    return s_use, psi_rel, s_ref_eff


# ---------- Plot-Helfer ----------
def plot_curves(s_grid: List[np.ndarray], y_grid: List[np.ndarray], labels: List[str],
                ylabel: str, title: str, fname: Path, yref: Optional[float] = None):
    plt.figure(figsize=(7.0, 4.8))
    for s, y, lab in zip(s_grid, y_grid, labels):
        plt.plot(s, y, label=lab)
    if yref is not None:
        plt.axhline(yref, color="k", lw=1.0, ls="--", alpha=0.6)
    plt.xlabel("s")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=160)
    plt.close()


def linear_fit(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    # y ≈ a*(1/L) + b;  R^2 zurückgeben
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = a*x + b
    ss_res = float(np.sum((y - yhat)**2))
    ss_tot = float(np.sum((y - np.mean(y))**2) + 1e-18)
    r2 = 1.0 - ss_res/ss_tot
    return float(a), float(b), float(r2)


# ---------- Hauptlogik ----------
def run_mode_csv(Ls: List[float], csv_pattern: str, s_ref: float, out_prefix: str) -> Dict:
    s_list, psi_rel_list, labels = [], [], []
    sup_dev, mean_dev, med_dev = [], [], []

    for L in Ls:
        path = Path(csv_pattern.format(L=int(L)))
        if not path.exists():
            raise FileNotFoundError(f"CSV nicht gefunden: {path}")
        s, zH = load_zetaH_table(path)
        s_use, psi_rel, s_ref_eff = build_psi_rel(s, zH, s_ref=s_ref)

        s_list.append(s_use)
        psi_rel_list.append(psi_rel)
        labels.append(f"L={int(L)}")

        dev = np.abs(psi_rel - 1.0)
        sup_dev.append(float(np.max(dev)))
        mean_dev.append(float(np.mean(dev)))
        med_dev.append(float(np.median(dev)))

    # Plots
    plot_curves(s_list, psi_rel_list, labels,
                ylabel="Ψ_rel(s)", title="Ψ_rel(s) je L", fname=Path(f"{out_prefix}_curves.png"), yref=1.0)

    dev_lists = [np.abs(y - 1.0) for y in psi_rel_list]
    plot_curves(s_list, dev_lists, labels,
                ylabel="|Ψ_rel(s) - 1|", title="Abweichungen |Ψ_rel-1|", fname=Path(f"{out_prefix}_dev_curves.png"))

    invL = 1.0/np.asarray(Ls, float)
    sup_dev = np.asarray(sup_dev, float)
    a, b, r2 = linear_fit(invL, sup_dev)

    # sup vs 1/L Plot
    xx = np.linspace(invL.min()*0.9, invL.max()*1.05, 100)
    yy = a*xx + b
    plt.figure(figsize=(6.8, 4.6))
    plt.plot(invL, sup_dev, "o", label="Daten")
    plt.plot(xx, yy, "-", label=f"Fit: a*(1/L)+b\n a={a:.4g}, b={b:.4g}, R²={r2:.4f}")
    plt.xlabel("1/L")
    plt.ylabel("sup_s |Ψ_rel(s)-1|")
    plt.title("Sup-Abweichung vs. 1/L")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_sup_vs_invL.png", dpi=160)
    plt.close()

    # Summary
    summary_df = pd.DataFrame({
        "L": Ls,
        "invL": invL,
        "sup_dev": sup_dev,
        "mean_dev": mean_dev,
        "median_dev": med_dev,
        "s_ref_eff": [float(s_ref) for _ in Ls],
    }).sort_values("L")
    summary_df.to_csv(f"{out_prefix}_summary.csv", index=False)

    txt = [
        "=== Proof 13 – L→∞ Limitcheck (CSV) ===",
        f"Files: {csv_pattern} für L in {Ls}",
        f"s_ref = {s_ref}",
        f"Linearfit: sup ≈ a*(1/L) + b  mit a={a:.6g}, b={b:.6g}, R²={r2:.6f}",
        "Tabelle:",
        summary_df.to_string(index=False),
    ]
    Path(f"{out_prefix}_summary.txt").write_text("\n".join(txt), encoding="utf-8")

    md = [
        "# Proof 13 – L→∞ Limitcheck (CSV)",
        f"- CSV-Muster: `{csv_pattern}`",
        f"- L-Werte: `{Ls}`",
        f"- s_ref: `{s_ref}`",
        f"- **Fit:** `sup ≈ a*(1/L) + b`  mit  `a={a:.6g}`, `b={b:.6g}`, `R²={r2:.6f}`",
        "",
        "## Kennzahlen",
        summary_df.to_markdown(index=False),
        "",
        "## Bewertung",
        "- **Konvergenz bestätigt**, wenn `b ≈ 0` (klein) und `R²` nahe 1, sowie die Kurven `Ψ_rel(s)` an 1 kleben.",
    ]
    Path(f"{out_prefix}_summary.md").write_text("\n".join(md), encoding="utf-8")

    return {"a": a, "b": b, "r2": r2}


def run_mode_evals(Ls: List[float], evals_pattern: str, smin: float, smax: float, Ns: int,
                   s_ref: float, out_prefix: str) -> Dict:
    # s-Gitter inkl. Fenster um 1/2; s=1 meiden
    s_coarse = np.linspace(max(0.60, smin), smax, Ns-12)
    s_fine   = np.linspace(0.46, 0.56, 13)
    s_fine   = s_fine[np.abs(s_fine - 0.50) > 1e-12]
    s_coarse = s_coarse[np.abs(s_coarse - 1.00) > 1e-12]
    s_all    = np.unique(np.concatenate([s_coarse, s_fine]))
    s_all.sort()

    s_list, psi_rel_list, labels = [], [], []
    sup_dev, mean_dev, med_dev = [], [], []

    for L in Ls:
        p = Path(evals_pattern.format(L=int(L)))
        if not p.exists():
            raise FileNotFoundError(f"Eigenwerte nicht gefunden: {p}")
        evals = np.loadtxt(p, dtype=float)
        evals = np.sort(evals[evals > 1e-14])
        zH = compute_zetaH_from_evals(evals, s_all, L=float(L), a1=None)

        s_use, psi_rel, s_ref_eff = build_psi_rel(s_all, zH, s_ref=s_ref)
        s_list.append(s_use)
        psi_rel_list.append(psi_rel)
        labels.append(f"L={int(L)}")

        dev = np.abs(psi_rel - 1.0)
        sup_dev.append(float(np.max(dev)))
        mean_dev.append(float(np.mean(dev)))
        med_dev.append(float(np.median(dev)))

        # optional Zwischen-CSV pro L
        pd.DataFrame({"s": s_use, "Psi_rel": psi_rel}).to_csv(f"{out_prefix}_PsiRel_L{int(L)}.csv", index=False)

    # Plots analog zu CSV-Modus
    plot_curves(s_list, psi_rel_list, labels, "Ψ_rel(s)", "Ψ_rel(s) je L", Path(f"{out_prefix}_curves.png"), yref=1.0)
    dev_lists = [np.abs(y - 1.0) for y in psi_rel_list]
    plot_curves(s_list, dev_lists, labels, "|Ψ_rel-1|", "Abweichungen |Ψ_rel-1|", Path(f"{out_prefix}_dev_curves.png"))

    invL = 1.0/np.asarray(Ls, float)
    sup_dev = np.asarray(sup_dev, float)
    a, b, r2 = linear_fit(invL, sup_dev)

    xx = np.linspace(invL.min()*0.9, invL.max()*1.05, 100)
    yy = a*xx + b
    plt.figure(figsize=(6.8, 4.6))
    plt.plot(invL, sup_dev, "o", label="Daten")
    plt.plot(xx, yy, "-", label=f"Fit: a*(1/L)+b\n a={a:.4g}, b={b:.4g}, R²={r2:.4f}")
    plt.xlabel("1/L")
    plt.ylabel("sup_s |Ψ_rel(s)-1|")
    plt.title("Sup-Abweichung vs. 1/L")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_sup_vs_invL.png", dpi=160)
    plt.close()

    summary_df = pd.DataFrame({
        "L": Ls,
        "invL": invL,
        "sup_dev": sup_dev,
        "mean_dev": mean_dev,
        "median_dev": med_dev,
        "s_ref_eff": [float(s_ref) for _ in Ls],
    }).sort_values("L")
    summary_df.to_csv(f"{out_prefix}_summary.csv", index=False)

    Path(f"{out_prefix}_summary.txt").write_text(
        "=== Proof 13 – L→∞ Limitcheck (EVALS) ===\n"
        f"Evals-Pattern: {evals_pattern}\n"
        f"L-Werte: {Ls}\n"
        f"s-Grid: [{smin}, {smax}], Ns={Ns}\n"
        f"s_ref = {s_ref}\n"
        f"Linearfit: a={a:.6g}, b={b:.6g}, R²={r2:.6f}\n\n"
        + summary_df.to_string(index=False),
        encoding="utf-8"
    )

    return {"a": a, "b": b, "r2": r2}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["csv", "evals"], default="csv")
    ap.add_argument("--Ls", type=float, nargs="+", required=True, help="Liste der L-Werte, z.B. 8 10 12 14")
    ap.add_argument("--csv-pattern", type=str, default="zetaH_L{L}.csv", help="Pfadmuster für CSVs (nur mode=csv)")
    ap.add_argument("--evals-pattern", type=str, default="eigs_L{L}.txt", help="Pfadmuster für Eigenwerte (nur mode=evals)")
    ap.add_argument("--smin", type=float, default=0.46, help="nur mode=evals")
    ap.add_argument("--smax", type=float, default=3.0,  help="nur mode=evals")
    ap.add_argument("--Ns",   type=int,   default=60,   help="nur mode=evals")
    ap.add_argument("--sref", type=float, default=3.0, help="Referenzpunkt für Ψ_rel")
    ap.add_argument("--out-prefix", type=str, default="proof13")
    args = ap.parse_args()

    if args.mode == "csv":
        run_mode_csv(args.Ls, args.csv_pattern, args.sref, args.out_prefix)
    else:
        run_mode_evals(args.Ls, args.evals_pattern, args.smin, args.smax, args.Ns, args.sref, args.out_prefix)


if __name__ == "__main__":
    main()
