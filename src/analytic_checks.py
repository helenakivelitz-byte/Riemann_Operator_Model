# src/analytic_checks.py
import numpy as np

class AnalyticPrecheck:
    """
    Lightweight boundedness & setup sanity check:
      - grid + step size
      - sup-norms of V_prim, V_grav[rho0], V_exch[rho0]
      - Kato–Rellich precondition note
    """
    def __init__(self, op):
        self.op = op

    def run(self):
        x = self.op.x
        L = self.op.L
        h = self.op.h

        # uniform density on [1,L], normalized by trapezoid rule
        rho0 = np.ones_like(x, dtype=float)
        rho0 /= np.trapz(rho0, x)

        # potentials (note: V_prim() takes NO args; V_grav/V_exch take rho)
        Vp = self.op.V_prim()             # shape (n,)
        Vg = self.op.V_grav(rho0)         # shape (n,)
        Ve = self.op.V_exch(rho0)         # shape (n,)

        info = {
            "L": float(L),
            "n_points": int(len(x)),
            "h": float(h),
            "V_prim_sup": float(np.max(np.abs(Vp))),
            "V_grav_sup_r0": float(np.max(np.abs(Vg))),
            "V_exch_sup_r0": float(np.max(np.abs(Ve))),
            "kato_rellich": "all V are bounded multiplications ⇒ relative bound 0",
        }

        # pretty print (matches your existing console style)
        print("0) Analytic pre-check")
        print(f"   • Domain [1,L] with L={L}, n_points={len(x)}, h={h:.4e}")
        print(f"   • ||V_prim||_∞  ≈ {info['V_prim_sup']: .3e}")
        print(f"   • ||V_grav[r0]||_∞ ≈ {info['V_grav_sup_r0']: .3e}")
        print(f"   • ||V_exch[r0]||_∞ ≈ {info['V_exch_sup_r0']: .3e}")
        print("   • Kato–Rellich precondition: all V are bounded multiplications ⇒ relative bound 0")

        return info
