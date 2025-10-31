# src/zeta_comparison.py
import numpy as np
import mpmath as mp

def _to_float(x):
    return float(mp.re(x)) if hasattr(x, 'real') else float(x)

class ZetaComparison:
    def __init__(self, eigenvalues, s_values, K=None):
        self.E = np.asarray(eigenvalues, dtype=float)
        self.s_values = s_values
        # --- NEW: record how many eigenvalues to include ---
        self.K = len(self.E) if K is None else min(int(K), len(self.E))

    def spectral_zeta(self, s):
        Epos = [mp.mpf(e) for e in self.E[:self.K] if e > 0]
        return mp.fsum(e**(-s) for e in Epos)


    def bridge_ratio(self):
        results = {}
        for s in self.s_values:
            s_c = complex(s) if not isinstance(s, complex) else s
            if s_c.real <= 1.0:
                print(f"⚠️ Skipping s={s}: outside safe domain Re(s)>1.")
                continue
            if abs(s_c - 1.0) < 1e-3:
                print(f"⚠️ Skipping s={s}: too close to the pole at 1.")
                continue

            zH = self.spectral_zeta(mp.mpc(s_c))
            zR = mp.zeta(mp.mpc(s_c))
            if abs(zR) < mp.mpf('1e-30'):
                print(f"⚠️ Skipping s={s}: |ζ(s)| too small.")
                continue

            psi = zH / zR
            results[str(s)] = {                     # keys as strings for JSON
                "truncation_K": int(self.K),             # <--- NEW
                "zeta_H": {"re": float(mp.re(zH)), "im": float(mp.im(zH))},
                "zeta_R": {"re": float(mp.re(zR)), "im": float(mp.im(zR))},
                "psi":    {"re": float(mp.re(psi)), "im": float(mp.im(psi))},
                "abs_deviation": float(abs(psi - 1))
            }
        return results
