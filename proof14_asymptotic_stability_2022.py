#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proof14_asymptotic_stability.py — Komfortversion
--------------------------------------------------------------
Beweis 14 — Asymptotic Stability & Uniform Convergence of Ψ(s;L)

Features:
  • Gemeinsames s-Grid (Shared Grid) über alle L.
  • Referenz-Normierung mit Auto-Snap (kein Absturz, wenn s_ref außerhalb).
  • Wahl der Fehlermetrik: sup(|Ψ_rel-1|) ODER log(|log|Ψ_rel||) → skaleninvariant.
  • Zwei-Punkt-Normierung (entfernt globalen Faktor + linearen s-Trend im log-Betrag).
  • Smarter CSV-Import (Trennzeichen ; , Tab, Whitespace; Dezimal-Komma; Re/Im o. komplexe Spalte).
  • Plots + CSV-Summary.

CLI-Beispiele:
  python proof14_asymptotic_stability.py --Ls 8 10 12 14 --smin 1.30 --smax 2.50 --sref 2.40 --plot_psirel
  python proof14_asymptotic_stability.py --Ls 8 10 12 14 --metric log --two_point 2.1 2.4
"""

from __future__ import annotations

import argparse
import importlib
from pathlib import Path
from typing import Dict, Tuple, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.typing import ArrayLike
from mpmath import zeta as riemann_zeta


# ------------------------------------------------------------
# Optional: rigoroser Brückenmodul (falls vorhanden)
# ------------------------------------------------------------
try:
    bridge = importlib.import_module("proof10_equivalence_rigorous_bridge")
    print("[proof14] Bridge module loaded successfully.")
except Exception as e:
    print("[proof14] Info: no bridge module found (optional):", e)
    bridge = None


# ============================================================
# 1) CSV-Loader: robust gegen Format-Varianz
# ============================================================

def _parse_complex_cell(val) -> complex:
    """String 'a+bj' o. ä. sicher in complex umwandeln (oder NaN+NaNj)."""
    if isinstance(val, complex):
        return val
    try:
        return complex(val)
    except Exception:
        return np.nan + 1j*np.nan

def _canon(name: str) -> str:
    """Spaltenname kanonisieren für Matching."""
    return "".join(ch for ch in str(name).lower() if ch.isalnum())

def _read_csv_smart(path: Path) -> pd.DataFrame:
    """
    Versucht mehrere (sep, decimal)-Kombinationen; wählt die beste (mehr
    numerische Zellen + genügend Spalten).
    """
    candidates = [
        (None, "."),     # sniff
        (None, ","),     # sniff + decimal comma
        (",",  "."),     # Komma
        (";",  "."),     # Semikolon (Excel DE)
        (";",  ","),     # Semikolon + Dezimal-Komma
        ("\t", "."),     # TSV
        ("\s+", "."),    # beliebiges Whitespace
    ]
    best_df, best_score = None, -1
    for sep, dec in candidates:
        try:
            df = pd.read_csv(path, sep=sep, engine="python", decimal=dec)
            num_cells = df.select_dtypes(include=["number"]).size
            score = (df.shape[1] * 10) + num_cells
            if score > best_score and df.shape[1] >= 2:
                best_df, best_score = df, score
        except Exception:
            continue
    if best_df is None:
        best_df = pd.read_csv(path)
    return best_df

def load_zetaH_table(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Akzeptierte Formate:
      • s , zetaH                     (zetaH reell oder komplexer String)
      • s , Re , Im                   (oder zH_real/zH_imag)
    """
    df = _read_csv_smart(path)
    original_cols = list(df.columns)
    cols_canon = {col: _canon(col) for col in df.columns}

    def find_col(*tokens, must_have=None, must_not=None):
        toks = [t.lower() for t in tokens if t]
        mh  = [m.lower() for m in (must_have or [])]
        mn  = [m.lower() for m in (must_not or [])]
        best = None
        for col, c in cols_canon.items():
            if all(t in c for t in toks) and all(m in c for m in mh) and all(m not in c for m in mn):
                best = col
                break
        return best

    # s-Spalte
    s_col = find_col("s") or find_col("sigma") or find_col("svalue", "s")
    if s_col is None:
        first = df.columns[0]
        try:
            _ = df[first].astype(float)
            s_col = first
        except Exception:
            raise ValueError(f"{path.name}: keine 's'-Spalte erkannt. Gefunden: {original_cols}")

    # Einzelspalte zetaH?
    z_single = (find_col("zetah", must_not=["re","im"]) or
                find_col("zh", "val", must_not=["re","im"]) or
                find_col("zeta", "h", must_not=["re","im"]))
    if z_single is not None:
        s = df[s_col].astype(float).to_numpy()
        z_series = df[z_single]
        if np.issubdtype(z_series.dtype, np.number):
            zH = z_series.to_numpy(dtype=float).astype(complex)
        else:
            zH = np.array([_parse_complex_cell(v) for v in z_series.to_numpy()], dtype=complex)
        print(f"[proof14] {path.name}: mapped columns → s='{s_col}', zetaH='{z_single}'")
        return s, zH

    # Re/Im?
    re_col = (find_col("zetah","re") or find_col("zh","re") or
              find_col("re", must_have=["z"]) or find_col("real", must_have=["z"]))
    im_col = (find_col("zetah","im") or find_col("zh","im") or
              find_col("im", must_have=["z"]) or find_col("imag", must_have=["z"]))
    if re_col is None or im_col is None:
        re_col = re_col or find_col("real")
        im_col = im_col or find_col("imag") or find_col("im")
    if re_col is not None and im_col is not None:
        s = df[s_col].astype(float).to_numpy()
        re = df[re_col].astype(float).to_numpy()
        im = df[im_col].astype(float).to_numpy()
        zH = re + 1j*im
        print(f"[proof14] {path.name}: mapped columns → s='{s_col}', Re='{re_col}', Im='{im_col}'")
        return s, zH

    # Fallback: numerische Spalten nach s
    df_num = df.drop(columns=[s_col], errors="ignore")
    numeric_cols = [c for c in df_num.columns if pd.api.types.is_numeric_dtype(df_num[c])]
    if len(numeric_cols) >= 2:
        s = df[s_col].astype(float).to_numpy()
        re = df_num[numeric_cols[0]].astype(float).to_numpy()
        im = df_num[numeric_cols[1]].astype(float).to_numpy()
        zH = re + 1j*im
        print(f"[proof14] {path.name}: fallback mapped → s='{s_col}', Re='{numeric_cols[0]}', Im='{numeric_cols[1]}'")
        return s, zH
    if len(numeric_cols) == 1:
        s = df[s_col].astype(float).to_numpy()
        zH = df_num[numeric_cols[0]].astype(float).to_numpy().astype(complex)
        print(f"[proof14] {path.name}: fallback mapped → s='{s_col}', zetaH='{numeric_cols[0]}' (real-only)")
        return s, zH

    raise ValueError(
        f"{path.name}: Spaltenformat nicht erkannt.\n"
        f"Gefunden: {original_cols}\n"
        f"Erwarte: (s,zetaH) ODER (s,Re,Im)."
    )


# ============================================================
# 2) Ψ(s)=ζ_H/ζ und Shared-Grid Utilities
# ============================================================

def compute_Psi(s: ArrayLike, zH: ArrayLike) -> np.ndarray:
    """Elementweise Ψ(s)=ζ_H(s)/ζ(s) (s reell)."""
    s = np.asarray(s, float)
    zH = np.asarray(zH, complex)
    out: List[complex] = []
    for si, zi in zip(s, zH):
        try:
            zeta_val = complex(riemann_zeta(si))
            out.append(zi / zeta_val)
        except Exception:
            out.append(np.nan + 1j*np.nan)
    return np.array(out, dtype=complex)

def load_all_L_tables(L_values: Tuple[int, ...],
                      pattern: str = "zetaH_L{L}.csv") -> Dict[int, Dict[str, np.ndarray]]:
    data: Dict[int, Dict[str, np.ndarray]] = {}
    for L in L_values:
        path = Path(pattern.format(L=L))
        if not path.exists():
            print(f"[proof14] Missing file: {path}")
            continue
        s, zH = load_zetaH_table(path)
        psi = compute_Psi(s, zH)
        data[L] = {"s": s, "zH": zH, "psi": psi}
        print(f"[proof14] Loaded {path.name} ({len(s)} points)")
    if not data:
        raise RuntimeError("Keine Daten geladen. Prüfe Pfade und Dateinamen.")
    return data

def shared_grid_from_data(data: Dict[int, Dict[str, np.ndarray]],
                          s_min=0.60, s_max=3.00, ds=0.05) -> np.ndarray:
    s_common = np.arange(s_min, s_max + 1e-12, ds)
    ok = np.ones_like(s_common, dtype=bool)
    for _, d in data.items():
        sL = d["s"]
        if len(sL) < 3:
            ok &= False
            continue
        ok &= (s_common >= sL.min()) & (s_common <= sL.max())
    return s_common[ok]

def interpolate_to_grid(s_src: ArrayLike, y_src: ArrayLike, s_tgt: ArrayLike) -> np.ndarray:
    s_src = np.asarray(s_src, float); s_tgt = np.asarray(s_tgt, float)
    y_src = np.asarray(y_src)
    if np.iscomplexobj(y_src):
        y_re = np.interp(s_tgt, s_src, np.real(y_src))
        y_im = np.interp(s_tgt, s_src, np.imag(y_src))
        return y_re + 1j*y_im
    return np.interp(s_tgt, s_src, y_src)

def compute_Psi_on_grid(d: Dict[str, np.ndarray], s_grid: np.ndarray) -> np.ndarray:
    zH_grid = interpolate_to_grid(d["s"], d["zH"], s_grid)
    return compute_Psi(s_grid, zH_grid)


# ============================================================
# 3) Hauptanalyse mit Komfort-Patches
# ============================================================

def analyse_convergence_refnorm(
    data: Dict[int, Dict[str, np.ndarray]],
    s_min=0.60, s_max=3.00, ds=0.05,
    s_ref=3.00,
    metric: str = "sup",           # "sup" oder "log"
    two_point: Optional[Tuple[float, float]] = None,
    robust: bool = False
) -> Tuple[pd.DataFrame, np.ndarray, Dict[int, np.ndarray]]:
    """
    Shared-Grid + Normierung + Δ_L-Auswertung.
      metric="sup":  err = |Ψ_rel-1|
      metric="log":  err = |log |Ψ_rel||  (skaleninvariant)
      two_point=(s1,s2): Zwei-Punkt-Normierung im log-Betrag (Trend-Korrektur)
    """
    s_grid = shared_grid_from_data(data, s_min, s_max, ds)
    if s_grid.size < 5:
        raise RuntimeError("Gemeinsames s-Gitter zu klein—prüfe s_min/s_max/ds und CSV-Abdeckung.")

    # Auto-Snap/Clamp von s_ref
    if not (s_grid[0] <= s_ref <= s_grid[-1]):
        s_ref_old = s_ref
        s_ref = float(np.clip(s_ref, s_grid[0], s_grid[-1]))
        print(f"[proof14] s_ref={s_ref_old} outside [{s_grid[0]}, {s_grid[-1]}] → using s_ref={s_ref} (auto-snap).")

    s_ref_idx = int(np.argmin(np.abs(s_grid - s_ref)))
    s_ref_eff = float(s_grid[s_ref_idx])

    rows: List[Dict[str, float]] = []
    psi_rel_store: Dict[int, np.ndarray] = {}

    for L in sorted(data.keys()):
        psi = compute_Psi_on_grid(data[L], s_grid)
        mask = np.isfinite(np.real(psi)) & np.isfinite(np.imag(psi))
        if mask.sum() < 5:
            print(f"[proof14] L={L}: zu wenige gültige Punkte, übersprungen.")
            continue

        # Normierung
        if two_point is not None:
            s1, s2 = two_point
            i1 = int(np.argmin(np.abs(s_grid - s1)))
            i2 = int(np.argmin(np.abs(s_grid - s2)))
            if i1 == i2:
                i2 = min(i1 + 1, len(s_grid) - 1)
            # log-Betrag linear approx: log|Ψ(s)| ≈ A + B*s
            y1 = np.log(np.maximum(1e-300, np.abs(psi[i1])))
            y2 = np.log(np.maximum(1e-300, np.abs(psi[i2])))
            B = (y2 - y1) / (s_grid[i2] - s_grid[i1] + 1e-15)
            A = y1 - B * s_grid[i1]
            psi_rel = psi * np.exp(-(A + B * s_grid))  # entfernt Skala + Trend
            ref_info = f"two-point({s_grid[i1]:.2f},{s_grid[i2]:.2f})"
        else:
            # Ein-Punkt-Normierung
            psi_ref = psi[s_ref_idx]
            if not np.isfinite(np.real(psi_ref)) or abs(psi_ref) < 1e-300:
                neigh = np.where(mask & (np.abs(s_grid - s_ref) <= 0.2))[0]
                if neigh.size == 0:
                    print(f"[proof14] L={L}: kein brauchbarer Referenzpunkt, übersprungen.")
                    continue
                psi_ref = psi[neigh[np.argmin(np.abs(s_grid[neigh] - s_ref))]]
            psi_rel = psi / psi_ref
            ref_info = f"s_ref={s_ref_eff:.2f}"

        # Fehlermetrik wählen
        if metric == "log":
            eps = 1e-15
            err_vec = np.abs(np.log(np.maximum(eps, np.abs(psi_rel[mask]))))
        else:
            err_vec = np.abs(psi_rel - 1.0)[mask]

        # Auswertung
        if robust:
            med = float(np.median(err_vec))
            mad = float(np.median(np.abs(err_vec - med)))
            delta = med + 1.4826 * mad
            delta_sup_print = float(np.max(err_vec))
            print(f"[proof14] L={L:3d} → robust Δ={delta:.4e} (sup={delta_sup_print:.4e}) @ {ref_info}")
        else:
            delta = float(np.max(err_vec))
            print(f"[proof14] L={L:3d} → Δ={delta:.4e} @ {ref_info}")

        psi_rel_store[L] = psi_rel
        rows.append({"L": float(L), "Delta_sup": delta})

    df = pd.DataFrame(rows).sort_values("L")
    return df, s_grid, psi_rel_store


# ============================================================
# 4) Fit & Plots
# ============================================================

def fit_asymptotic_decay(df: pd.DataFrame) -> Dict[str, float]:
    Ls = df["L"].to_numpy(dtype=float)
    invL = 1.0 / Ls
    Δ = df["Delta_sup"].to_numpy(dtype=float)
    if Δ.size < 2:
        raise RuntimeError("Zu wenige Punkte für Fit.")
    coeffs = np.polyfit(invL, Δ, 1)
    a, b = float(coeffs[0]), float(coeffs[1])
    pred = np.polyval(coeffs, invL)
    ss_res = float(np.sum((Δ - pred)**2))
    ss_tot = float(np.sum((Δ - np.mean(Δ))**2))
    R2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 1.0
    print(f"[proof14] Fit: Δ_L ≈ {a:.4f}/L + {b:.4f}, R²={R2:.3f}")
    return {"a": a, "b": b, "R2": R2}

def plot_convergence(df: pd.DataFrame, fit_params: Dict[str, float], out_prefix="proof14"):
    invL = 1.0 / df["L"].to_numpy(dtype=float)
    Δ = df["Delta_sup"].to_numpy(dtype=float)
    a, b = fit_params["a"], fit_params["b"]
    x_fit = np.linspace(0, max(invL)*1.1, 200)
    y_fit = a * x_fit + b

    plt.figure(figsize=(6.0, 4.2))
    plt.plot(invL, Δ, "o-", label="Δ_L (sup |Ψ_rel-1|)")
    plt.plot(x_fit, y_fit, "--", label=f"Fit: {a:.3g}/L + {b:.3g}")
    plt.xlabel("1 / L"); plt.ylabel("Δ_L")
    plt.title("Proof 14 — Convergence of Ψ(s;L) on shared grid")
    plt.grid(True, alpha=0.35); plt.legend(); plt.tight_layout()
    plt.savefig(f"{out_prefix}_sup_vs_invL.png", dpi=160)
    plt.close()
    print(f"[proof14] Saved plot: {out_prefix}_sup_vs_invL.png")

def plot_psirel_grid(s_grid: np.ndarray,
                     psi_rel_store: Dict[int, np.ndarray],
                     out_prefix="proof14"):
    plt.figure(figsize=(7.2, 4.2))
    for L in sorted(psi_rel_store.keys()):
        psi_rel = psi_rel_store[L]
        plt.plot(s_grid, np.abs(psi_rel), label=f"L={L}")
    plt.xlabel("s"); plt.ylabel("|Ψ_rel(s;L)|")
    plt.title("Proof 14 — |Ψ_rel(s;L)| across the shared s-grid")
    plt.grid(True, alpha=0.35); plt.legend(ncol=2, fontsize=9)
    plt.tight_layout(); plt.savefig(f"{out_prefix}_psirel_grid.png", dpi=160)
    plt.close()
    print(f"[proof14] Saved plot: {out_prefix}_psirel_grid.png")


# ============================================================
# 5) Analytischer Platzhalter
# ============================================================

def resolvent_stability_placeholder(L_values: Tuple[int, ...]):
    for L in L_values:
        print(f"[proof14] Heuristik: ‖ΔR(λ;L)‖₁ ~ O(1/{L}) für λ außerhalb [0,∞).")


# ============================================================
# 6) CLI & Main
# ============================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Proof 14 — Asymptotic Stability of Ψ(s;L)")
    p.add_argument("--Ls", nargs="+", type=int, default=[8, 10, 12, 14],
                   help="Liste der L-Werte (z.B. --Ls 8 10 12 14)")
    p.add_argument("--pattern", type=str, default="zetaH_L{L}.csv",
                   help="Dateimuster der CSVs (nutzt {L} als Platzhalter)")
    p.add_argument("--smin", type=float, default=0.60, help="untere s-Grenze")
    p.add_argument("--smax", type=float, default=3.00, help="obere s-Grenze")
    p.add_argument("--ds", type=float, default=0.05, help="s-Schrittweite für Shared Grid")
    p.add_argument("--sref", type=float, default=3.00, help="Referenz s_ref für Ein-Punkt-Normierung")
    p.add_argument("--metric", choices=["sup", "log"], default="sup",
                   help="Fehlermetrik: 'sup' (|Ψ_rel-1|) oder 'log' (|log|Ψ_rel||).")
    p.add_argument("--two_point", nargs=2, type=float, metavar=("S1","S2"),
                   help="Zwei-Punkt-Normierung im log-Betrag bei s=S1,S2 (z.B. --two_point 2.1 2.4).")
    p.add_argument("--robust", action="store_true",
                   help="robustes Δ (Median + 1.4826*MAD) statt reiner sup.")
    p.add_argument("--plot_psirel", action="store_true",
                   help="Plot von |Ψ_rel(s;L)| erzeugen.")
    return p.parse_args()

def main():
    print("=== Proof 14: Asymptotic Stability of Ψ(s;L) ===")
    args = parse_args()

    data = load_all_L_tables(tuple(args.Ls), pattern=args.pattern)

    df, s_grid, psi_rel_store = analyse_convergence_refnorm(
        data,
        s_min=args.smin, s_max=args.smax, ds=args.ds,
        s_ref=args.sref,
        metric=args.metric,
        two_point=tuple(args.two_point) if args.two_point is not None else None,
        robust=args.robust
    )
    if len(df) < 2:
        print("[proof14] Zu wenige L-Punkte für Regression. Ausgabe nur Rohwerte:")
        print(df)
        return

    fit_params = fit_asymptotic_decay(df)
    plot_convergence(df, fit_params, out_prefix="proof14")
    if args.plot_psirel:
        plot_psirel_grid(s_grid, psi_rel_store, out_prefix="proof14")

    resolvent_stability_placeholder(tuple(sorted(psi_rel_store.keys())))

    df_out = df.copy()
    for k, v in fit_params.items():
        df_out[k] = v
    df_out["s_ref"] = args.sref
    df_out["s_min"] = args.smin
    df_out["s_max"] = args.smax
    df_out["ds"] = args.ds
    df_out["metric"] = args.metric
    df_out["two_point"] = "" if args.two_point is None else f"{args.two_point[0]};{args.two_point[1]}"
    df_out["robust"] = args.robust
    df_out.to_csv("proof14_summary.csv", index=False)
    print("[proof14] Exported summary → proof14_summary.csv")

    print("\n=== Summary ===")
    print(df_out.to_string(index=False))
    print("\nInterpretation:")
    print("  • Uniforme Konvergenz plausibel, wenn b ≈ 0 und R² hoch (≳ 0.9).")
    print("  • 'log'-Metrik: robust gg. globale Skalenfehler.")
    print("  • 'two_point S1 S2': entfernt zusätzlich linearen s-Trend im log-Betrag.")
    print("  • Gleiche s-Gitter/Parameter bei der Erzeugung der zetaH_L*.csv verwenden.")

if __name__ == "__main__":
    main()
