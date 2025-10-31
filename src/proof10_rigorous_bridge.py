# src/proof10_rigorous_bridge.py
# Rigorous heat-trace based spectral zeta bridge:
#   zeta_H(s) = (1/Gamma(s)) * ∫_0^∞ t^{s-1} Tr(e^{-tH}) dt
# with small-t subtraction of Weyl terms (a0 t^{-1/2} + a1).

from __future__ import annotations
import os, json, math, argparse
from typing import Iterable, Dict, Any, Tuple

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt


# ---------- Utility: IO ----------

def _load_eigenvalues() -> np.ndarray:
    """
    Load eigenvalues from data/eigenvalues_raw.json,
    fallback to data/results.json (key 'eigenvalues').
    """
    paths = ["data/eigenvalues_raw.json", "data/results.json"]
    for p in paths:
        if os.path.exists(p):
            with open(p, "r", encoding="utf-8") as f:
                obj = json.load(f)
            # eigenvalues_raw.json -> {"eigenvalues": [...]}
            # results.json -> may contain "eigenvalues" as well
            if isinstance(obj, dict) and "eigenvalues" in obj:
                return np.array(obj["eigenvalues"], dtype=float)
            # Or the file *is* just a list already:
            if isinstance(obj, list):
                return np.array(obj, dtype=float)
    raise FileNotFoundError("Could not find eigenvalues in data/.")


def _ensure_dirs():
    os.makedirs("data", exist_ok=True)
    os.makedirs("data/diagnostics", exist_ok=True)
    os.makedirs("data/figures", exist_ok=True)


# ---------- Heat-trace, Weyl coefficients ----------

def a0_coefficient(L: float) -> mp.mpf:
    """
    1D Dirichlet Laplacian leading heat-trace coefficient:
      Tr(e^{-tH0}) ~ a0 t^{-1/2} + a1 + ...
    with a0 = (L-1)/(2√π)
    (Here the domain is [1, L] ⇒ length = L-1)
    """
    return mp.mpf(L - 1) / (2 * mp.sqrt(mp.pi))


def heat_trace_from_eigs(E: np.ndarray, t: mp.mpf, K: int | None = None) -> mp.mpf:
    """
    Compute Tr(e^{-tH}) from discrete spectrum (positive eigenvalues).
    """
    if K is None or K > len(E):
        K = len(E)
    total = mp.mpf("0")
    tt = mp.mpf(t)
    for lam in E[:K]:
        if lam > 0:
            total += mp.e ** (-tt * mp.mpf(lam))
    return total


# ---------- a1 estimation from eigenvalues (small-t fit) ----------

def estimate_a1_from_eigs(E: np.ndarray, tmin: float = 1e-4, tmax: float = 5e-3, m: int = 40, K: int | None = None) -> float:
    """
    Fit small-t expansion:
        Tr(e^{-tH}) ~ a0 * t^{-1/2} + a1
    on a log-spaced grid t ∈ [tmin, tmax] to obtain a1.
    """
    ts = np.geomspace(tmin, tmax, m)
    # Design matrix columns: t^{-1/2}, 1
    X = np.column_stack([ts ** (-0.5), np.ones_like(ts)])

    y = np.zeros_like(ts, dtype=float)
    for i, t in enumerate(ts):
        y[i] = float(heat_trace_from_eigs(E, mp.mpf(t), K=K))

    # Least squares fit
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)  # [a0_fit, a1_fit]
    a1_fit = float(coef[1])
    return a1_fit


# ---------- Quadratures for the zeta integral ----------

def small_t_integral_subtracted(
    E: np.ndarray,
    s: mp.mpf,
    L: float,
    a1: mp.mpf,
    t0: float,
    n_small: int,
    K: int | None = None,
) -> mp.mpf:
    """
    ∫_0^{t0} t^{s-1} [Tr(e^{-tH}) - a0 t^{-1/2} - a1] dt
    (Weyl-subtracted small-t integral)
    """
    a0 = a0_coefficient(L)
    t0_mp = mp.mpf(t0)
    # Gauss-Legendre on [0, t0] via substitution u∈[-1,1] -> t(u)
    xs, ws = np.polynomial.legendre.leggauss(n_small)
    total = mp.mpf("0")
    for x, w in zip(xs, ws):
        # affine map [-1,1] -> [0, t0]
        t = 0.5 * (x + 1.0) * t0_mp
        if t <= 0:
            continue
        tr = heat_trace_from_eigs(E, t, K=K)
        integrand = mp.power(t, s - 1) * (tr - a0 * mp.power(t, -mp.mpf("0.5")) - a1)
        total += mp.mpf(w) * integrand
    total *= 0.5 * t0_mp
    return total


def large_t_integral(
    E: np.ndarray,
    s: mp.mpf,
    t0: float,
    n_large: int,
    u_max: float,
    K: int | None = None,
) -> mp.mpf:
    """
    ∫_{t0}^{∞} t^{s-1} Tr(e^{-tH}) dt
    with substitution t = t0 * exp(u), u ∈ [0, u_max].
    """
    t0_mp = mp.mpf(t0)

    # Gauss-Legendre on [0, u_max]
    xs, ws = np.polynomial.legendre.leggauss(n_large)
    # affine map [-1,1] -> [0, u_max]
    total = mp.mpf("0")
    for x, w in zip(xs, ws):
        u = 0.5 * (x + 1.0) * u_max
        t = t0_mp * mp.e ** (mp.mpf(u))
        tr = heat_trace_from_eigs(E, t, K=K)
        integrand = mp.power(t, s - 1) * tr * (t)  # dt = t du
        total += mp.mpf(w) * integrand
    total *= 0.5 * mp.mpf(u_max)
    return total


def small_t_closedform(s: mp.mpf, L: float, a1: mp.mpf, t0: float) -> mp.mpf:
    """
    Closed form for the subtracted part ∫_0^{t0} t^{s-1}(a0 t^{-1/2} + a1) dt.
    """
    a0 = a0_coefficient(L)
    t0_mp = mp.mpf(t0)
    term_a0 = a0 * mp.power(t0_mp, s - mp.mpf("0.5")) / (s - mp.mpf("0.5"))
    term_a1 = a1 * mp.power(t0_mp, s) / s
    return term_a0 + term_a1


def zeta_from_heattrace(
    E: np.ndarray,
    L: float,
    s_values: Iterable[complex],
    a1_value: float | None = None,
    t0: float = 1e-2,
    n_small: int = 300,
    n_large: int = 800,
    u_max: float = 40.0,
    K: int | None = None,
) -> Dict[str, complex]:
    """
    Compute ζ_H(s) via heat-trace splitting + Weyl subtraction.
    Returns dict mapping s_str -> complex ζ_H(s).
    """
    mp.mp.dps = 50
    if K is not None:
        E = np.sort(np.asarray(E, dtype=float))[:K]

    if a1_value is None:
        # auto-fit a1 from eigenvalues on small t interval
        try:
            a1_value = estimate_a1_from_eigs(E, K=K)
            print(f"• a1 auto-fit from eigenvalues: {a1_value:.6e}")
        except Exception as exc:
            print(f"• a1 auto-fit failed ({exc}), using 0.0")
            a1_value = 0.0

    a1 = mp.mpf(a1_value)

    res: Dict[str, complex] = {}
    for s in s_values:
        s_mp = mp.mpf(s) if isinstance(s, (int, float)) else mp.mpc(s)
        pref = 1.0 / mp.gamma(s_mp)
        I_small = small_t_integral_subtracted(E, s_mp, L, a1, t0, n_small, K=K)
        I_large = large_t_integral(E, s_mp, t0, n_large, u_max, K=K)
        I_closed = small_t_closedform(s_mp, L, a1, t0)
        zH = pref * (I_small + I_large + I_closed)
        res[str(s)] = complex(zH)
    return res


# ---------- Residue & symmetry diagnostics ----------

def residue_at_half(zh: Dict[str, complex], s0: float = 0.5, eps: float = 1e-6) -> complex:
    """
    Estimate Res_{s=1/2} ζ_H(s) via symmetric difference:
        Res ≈ (ζ(1/2+ε) - ζ(1/2-ε)) / (2ε)
    """
    sp = str(s0 + eps)
    sm = str(s0 - eps)
    if sp in zh and sm in zh:
        return (zh[sp] - zh[sm]) / (2.0 * eps)
    return complex("nan")


def symmetry_proxy(zh: Dict[str, complex], s0: float = 0.5, eps: float = 1e-3) -> float:
    """
    Check | ζ(1/2+ε) - conj(ζ(1/2-ε)) | as a crude symmetry proxy.
    """
    sp = str(s0 + eps)
    sm = str(s0 - eps)
    if sp in zh and sm in zh:
        return abs(zh[sp] - np.conj(zh[sm]))
    return float("nan")


# ---------- Main driver ----------

def main():
    parser = argparse.ArgumentParser(description="Proof10: rigorous heat-trace spectral zeta bridge")
    parser.add_argument("--L", type=float, default=12.0, help="Domain end (interval [1,L])")
    parser.add_argument("--a1", type=str, default="auto", help="a1 value (float) or 'auto'")
    parser.add_argument("--t0", type=float, default=1e-2, help="split point for heat-trace integrals")
    parser.add_argument("--smin", type=float, default=1.05)
    parser.add_argument("--smax", type=float, default=1.60)
    parser.add_argument("--sstep", type=float, default=0.01)
    parser.add_argument("--nsmall", type=int, default=300, help="GL nodes small-t")
    parser.add_argument("--nlarge", type=int, default=800, help="GL nodes large-t")
    parser.add_argument("--umax", type=float, default=40.0, help="u-max for large-t substitution")
    parser.add_argument("--K", type=int, default=200, help="truncation of eigenvalues for sums/traces")
    args = parser.parse_args()

    _ensure_dirs()
    E_all = _load_eigenvalues()
    if args.K is not None and args.K > 0:
        E = np.sort(np.asarray(E_all, dtype=float))[: args.K]
    else:
        E = np.sort(np.asarray(E_all, dtype=float))

    if args.a1.strip().lower() == "auto":
        a1_val: float | None = None
        a1_mode = "auto-fit from eigenvalues"
    else:
        a1_val = float(args.a1)
        a1_mode = f"fixed={a1_val}"

    print("🚀 Proof10 — Rigorous heat-trace zeta bridge")
    print("============================================================")
    print(f"• L = {args.L}, K = {len(E)} eigenvalues")
    print(f"• a1 mode: {a1_mode}")
    print(f"• split t0 = {args.t0}, n_small = {args.nsmall}, n_large = {args.nlarge}, u_max = {args.umax}")

    # s-grid for |Psi(s)-1| in Re(s) > 1
    s_grid = np.round(np.arange(args.smin, args.smax + 1e-12, args.sstep), 12)
    # add the pair around 1/2 for residue/symmetry
    s_diag = [0.5 - 1e-6, 0.5 + 1e-6]

    # compute zeta on union of sets
    s_values_for_eval = list(map(float, s_grid)) + s_diag
    zH = zeta_from_heattrace(
        E,
        args.L,
        s_values_for_eval,
        a1_value=a1_val,
        t0=args.t0,
        n_small=args.nsmall,
        n_large=args.nlarge,
        u_max=args.umax,
        K=len(E),
    )

    # Bridge ratio and |Psi-1|
    mp.mp.dps = 50
    absdev = []
    for s in s_grid:
        zH_s = zH[str(s)]
        zR_s = mp.zeta(s)
        psi = zH_s / complex(zR_s)
        absdev.append(abs(psi - 1.0))

    # Residue & symmetry
    res_est = residue_at_half(zH, s0=0.5, eps=1e-6)
    res_theory = (args.L - 1.0) / (2.0 * math.pi)
    symm_err = symmetry_proxy(zH, s0=0.5, eps=1e-6)

    # Save JSON
    out_json = {
        "L": args.L,
        "K": int(len(E)),
        "a1_mode": a1_mode,
        "t0": args.t0,
        "n_small": args.nsmall,
        "n_large": args.nlarge,
        "u_max": args.umax,
        "s_grid": list(map(float, s_grid)),
        "absdev": list(map(float, absdev)),
        "residue_half_est": [float(np.real(res_est)), float(np.imag(res_est))],
        "residue_half_theory": float(res_theory),
        "symmetry_proxy_eps1e6": float(symm_err),
        "zetaH_samples": {k: [float(np.real(v)), float(np.imag(v))] for k, v in zH.items()},
    }
    with open("data/diagnostics/proof10_bridge.json", "w", encoding="utf-8") as f:
        json.dump(out_json, f, indent=2)
    print("💾 Saved: data/diagnostics/proof10_bridge.json")

    # ---- Plots ----
    # 1) |Psi-1| vs s (linear y)
    plt.figure(figsize=(9, 6))
    plt.plot(s_grid, absdev, marker="o")
    plt.xlabel("s")
    plt.ylabel(r"$|\Psi(s;L)-1|$")
    plt.title(r"Bridge Ratio Deviation $|\Psi(s;L)-1|$ (heat-trace, truncated)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("data/figures/proof10_zeta_absdev.png", dpi=140)
    plt.close()

    # 2) Residue check at s=1/2
    plt.figure(figsize=(8, 5))
    est_real = np.real(res_est)
    plt.axhline(res_theory, color="g", linestyle="--", label=r"theory $(L-1)/(2\pi)$")
    plt.scatter([0.5], [est_real], label="estimate (real part)", zorder=5)
    plt.xlabel("s = 1/2")
    plt.ylabel("Residue (real part)")
    plt.title("Residue check at s=1/2")
    plt.legend()
    plt.tight_layout()
    plt.savefig("data/figures/proof10_residue_check.png", dpi=140)
    plt.close()

    # 3) Symmetry proxy magnitude
    plt.figure(figsize=(8, 5))
    plt.bar([0], [symm_err])
    plt.xticks([0], [r"$|\zeta_H(1/2+\varepsilon) - \overline{\zeta_H(1/2-\varepsilon)}|$"])
    plt.ylabel("magnitude")
    plt.title("Symmetry proxy near the critical line")
    plt.tight_layout()
    plt.savefig("data/figures/proof10_symmetry_proxy.png", dpi=140)
    plt.close()

    print("🖼  Plots:")
    print("   - data/figures/proof10_zeta_absdev.png")
    print("   - data/figures/proof10_residue_check.png")
    print("   - data/figures/proof10_symmetry_proxy.png")
    print("✅ Proof10 completed.")


if __name__ == "__main__":
    main()
