# src/paper_plots.py
import os, sys, json, argparse
import numpy as np
import matplotlib.pyplot as plt

# --- robust relative import of diagnostics_utils ---
try:
    from .diagnostics_utils import (
        load_results, ensure_dirs, dense_s_grid,
        eigenvalues_from_results, bridge_ratio_absdev
    )
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    from diagnostics_utils import (
        load_results, ensure_dirs, dense_s_grid,
        eigenvalues_from_results, bridge_ratio_absdev
    )

def _load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def _save_fig(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.savefig(path, dpi=200, bbox_inches="tight")
    print(f"✅ Saved figure: {path}")

# ---------- helpers to coerce various JSON shapes to arrays ----------
def _to_float_array(x):
    """Coerce list-like of numbers to np.array(float)."""
    return np.array(list(x), dtype=float)

def _coerce_series_from_dict_with_sgrid(sgrid, series_dict):
    """series_dict: mapping 's'->value (keys may be str/float)."""
    # normalize keys to strings with same repr used in JSON
    out = []
    for s in sgrid:
        # try exact float key, str key, and formatted str key
        candidates = [s, str(s), f"{s:.12g}"]
        val = None
        for k in candidates:
            if k in series_dict:
                val = series_dict[k]
                break
        if val is None:
            raise KeyError(f"Missing value for s={s} in diagnostics dict.")
        out.append(float(val))
    return np.array(out, dtype=float)

def _coerce_from_points(points):
    """points: list of {'s':..., 'raw':..., 'rescaled':...}"""
    s = [float(p["s"]) for p in points]
    idx = np.argsort(s)
    s = np.array(s, dtype=float)[idx]
    raw = np.array([float(points[i]["raw"]) for i in idx], dtype=float)
    res = np.array([float(points[i]["rescaled"]) for i in idx], dtype=float)
    return s, raw, res

# ---------- primary loaders ----------
def load_dense_from_diagnostics():
    """
    Preferred: data/diagnostics/bridge_dense.json
    Accepts multiple shapes:
      1) {'sgrid':[...], 'absdev_raw':[...], 'absdev_rescaled':[...]}
      2) same but with dicts for absdev_*
      3) {'points':[{'s':..,'raw':..,'rescaled':..}, ...]}
    Returns: (s, absdev_raw, absdev_res) as np.float arrays, or None
    """
    path = os.path.join("data", "diagnostics", "bridge_dense.json")
    if not os.path.exists(path):
        return None
    D = _load_json(path)

    # shape (3): points list
    if "points" in D and isinstance(D["points"], list):
        try:
            return _coerce_from_points(D["points"])
        except Exception:
            return None

    # shapes (1) / (2): sgrid + two series
    sgrid = D.get("sgrid")
    if sgrid is None:
        return None
    s = _to_float_array(sgrid)

    raw = D.get("absdev_raw")
    res = D.get("absdev_rescaled")
    if raw is None or res is None:
        return None

    # case: arrays already
    if isinstance(raw, (list, tuple)) and isinstance(res, (list, tuple)):
        try:
            return s, _to_float_array(raw), _to_float_array(res)
        except Exception:
            return None

    # case: dicts keyed by s
    if isinstance(raw, dict) and isinstance(res, dict):
        try:
            raw_arr = _coerce_series_from_dict_with_sgrid(s, raw)
            res_arr = _coerce_series_from_dict_with_sgrid(s, res)
            return s, raw_arr, res_arr
        except Exception:
            return None

    return None

def load_sparse_from_summary():
    """Fallback: parse sparse (3 points) from data/paper_summary.json."""
    path = os.path.join("data", "paper_summary.json")
    if not os.path.exists(path):
        return None
    S = _load_json(path)
    zr, zres = S.get("zeta_raw", {}), S.get("zeta_rescaled", {})
    s_raw, a_raw = zr.get("grid", []), zr.get("absdev", [])
    s_res, a_res = zres.get("grid", []), zres.get("absdev", [])
    if not s_raw or not a_raw or not s_res or not a_res:
        return None
    s_common = sorted(set(map(float, s_raw)).intersection(map(float, s_res)))
    if not s_common:
        return None

    def pick(vals, grid_src, grid_tgt):
        idx = {float(ss): i for i, ss in enumerate(grid_src)}
        return [float(vals[idx[float(ss)]]) for ss in grid_tgt]

    s = np.array(s_common, dtype=float)
    absdev_raw = np.array(pick(a_raw, list(map(float, s_raw)), s_common))
    absdev_res = np.array(pick(a_res, list(map(float, s_res)), s_common))
    return s, absdev_raw, absdev_res

# ---------- plotting ----------
def plot_bridge_curve(s, y, title, ylabel, outpath, overwrite=False):
    if (not overwrite) and os.path.exists(outpath):
        print(f"↪︎ Skipping (exists): {outpath}")
        return
    plt.figure(figsize=(9, 6))
    plt.plot(s, y, marker="o", linewidth=2)
    plt.xlabel("s")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.35)
    _save_fig(outpath)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Create publication plots for bridge ratio.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing images.")
    args = parser.parse_args()

    ensure_dirs()
    fig_dir = os.path.join("data", "figures")
    os.makedirs(fig_dir, exist_ok=True)

    # 1) Dense diagnostics (preferred)
    loaded = load_dense_from_diagnostics()
    if loaded is not None:
        s, absdev_raw, absdev_res = loaded
        print(f"• Using diagnostics grid with {len(s)} s-values.")
    else:
        # 2) Sparse summary
        loaded = load_sparse_from_summary()
        if loaded is not None:
            s, absdev_raw, absdev_res = loaded
            print(f"• Using sparse summary with {len(s)} s-values.")
        else:
            # 3) Recompute from results (guaranteed fallback)
            print("• Diagnostics not found. Recomputing dense grid from results…")
            RAW, RES = load_results()
            E = eigenvalues_from_results(RAW)
            s = dense_s_grid(RAW, start=1.05, stop=1.50, step=0.01)
            # truncation if available
            K_res = None
            try:
                K_res = int(RES.get("bridge_ratio_rescaled", {}).get("truncation_K", 0))
                if K_res <= 0:
                    K_res = None
            except Exception:
                K_res = None
            absdev_raw = np.array([bridge_ratio_absdev(E, float(ss), k=None) for ss in s])
            absdev_res  = np.array([bridge_ratio_absdev(E, float(ss), k=K_res) for ss in s])

    # --- make plots ---
    plot_bridge_curve(
        s, absdev_raw,
        title="Bridge Ratio Before Rescaling",
        ylabel=r"$|\Psi_{\mathrm{raw}}(s;L) - 1|$",
        outpath=os.path.join(fig_dir, "psi_raw_vs_s.png"),
        overwrite=args.overwrite
    )
    plot_bridge_curve(
        s, absdev_res,
        title="Bridge Ratio After Rescaling (Berry–Keating / OT)",
        ylabel=r"$|\Psi_{\mathrm{rescaled}}(s;L) - 1|$",
        outpath=os.path.join(fig_dir, "psi_rescaled_vs_s.png"),
        overwrite=args.overwrite
    )
    print("🎨 Plot generation finished.")

if __name__ == "__main__":
    main()
