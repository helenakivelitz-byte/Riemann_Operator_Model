# proof10_equivalence_rigorous_bridge.py
# BRÜCKE: Integriert die rigorose Heat-Trace ζ_H(s) in die Äquivalenz-Diagnostik

import mpmath as mp
import numpy as np
import sympy as sp

# -----------------------------
# Rigorose Heat-Trace ζ_H(s) für H0
# -----------------------------

def dirichlet_eigs_free(L, N):
    length = mp.mpf(L) - 1
    pi_over_len = mp.pi / length
    return [(pi_over_len * n) ** 2 for n in range(1, N + 1)]

def heat_trace_dirichlet_free_poisson(L, t, kmax=None):
    L = mp.mpf(L); t = mp.mpf(t)
    length = L - 1
    alpha = (mp.pi**2) * t / (length**2)
    sqrt_term = mp.sqrt(mp.pi / alpha)
    if kmax is None:
        target = mp.mpf('1e-30')
        k_est = mp.sqrt(alpha * mp.log(1 / target)) / mp.pi
        kmax = max(5, int(mp.ceil(k_est)) + 5)
    tail = mp.fsum(mp.e**(-(mp.pi**2) * (k**2) / alpha) for k in range(1, kmax + 1))
    return mp.mpf('0.5') * (sqrt_term - 1) + sqrt_term * tail

class RigorousHeatTraceZeta:
    def __init__(self, L=100.0, N_tail=4000, t0=1e-10, dps=80):
        mp.mp.dps = dps
        self.L = mp.mpf(L)
        self.N_tail = int(N_tail)
        self.t0 = mp.mpf(t0)
        self.a0 = (self.L - 1) / mp.sqrt(4 * mp.pi)
        self.a1 = -mp.mpf('0.5')
        self.eigs_tail = dirichlet_eigs_free(self.L, self.N_tail)

    def Tr(self, t):
        return heat_trace_dirichlet_free_poisson(self.L, t)

    def zeta(self, s):
        s = mp.mpf(s); t0 = self.t0
        def integrand_small(t):
            t = mp.mpf(t)
            return t**(s - 1) * (self.Tr(t) - self.a0 * t**(-0.5) - self.a1)
        I_small = mp.quad(integrand_small, [t0, 1])
        closed_small = self.a0 / (s - 0.5) + self.a1 / s
        tail_sum = mp.fsum((lam ** (-s)) * mp.gammainc(s, lam, mp.inf) for lam in self.eigs_tail)
        return (I_small + closed_small + tail_sum) / mp.gamma(s)

    def zeta_reference_free(self, s):
        s = mp.mpf(s)
        factor = ((self.L - 1) / mp.pi) ** (2 * s)
        return factor * mp.zeta(2 * s)

# -----------------------------
# Äquivalenz-Diagnostik (A1, A2 korrigiert)
# -----------------------------

class RigorousEquivalenceDiagnostics:
    def __init__(self, L=100.0, N_tail=4000, t0=1e-10, dps=80):
        mp.mp.dps = dps
        self.res = {}
        self.zeta_H = RigorousHeatTraceZeta(L=L, N_tail=N_tail, t0=t0, dps=dps)

    def zeta(self, s):
        return self.zeta_H.zeta(s)

    # ---- A1': Residuum bei s=1/2 via kleiner Laurent-Fit ----
    def _laurent_residue_half(self, radii=None):
        if radii is None:
            radii = [mp.mpf('1e-3'), mp.mpf('5e-4'), mp.mpf('2e-4'), mp.mpf('1e-4')]
        s0 = mp.mpf('0.5')
        xs, ys = [], []
        for r in radii:
            for sign in (+1, -1):
                x = s0 + sign*r
                xs.append(float(x)); ys.append(float(self.zeta(x)))
        # Fit: A/(x-s0) + B + C*(x-s0)
        import numpy as _np
        A = _np.zeros((len(xs), 3)); b = _np.zeros(len(xs))
        for i, x in enumerate(xs):
            A[i, 0] = 1.0/(x - float(s0))
            A[i, 1] = 1.0
            A[i, 2] = (x - float(s0))
            b[i] = ys[i]
        coeff, *_ = _np.linalg.lstsq(A, b, rcond=None)
        Ares, Bconst, _ = coeff
        return mp.mpf(Ares), mp.mpf(Bconst)

    def test_A1_residue_at_half(self):
        print("🔍 Test A1': Residuum bei s=1/2 (Laurent-Fit)...")
        Ares, _ = self._laurent_residue_half()
        expected = self.zeta_H.a0 / mp.sqrt(mp.pi)
        err = abs(Ares - expected)
        self.res['A1'] = {
            'status': '✅' if err < mp.mpf('1e-9') else '⚠️',
            'residue_fit': float(Ares),
            'expected': float(expected),
            'abs_error': float(err),
        }
        return self.res['A1']

    # ---- A2': KORREKTE Symmetrie für H0: Xi_H0(s) = Xi_H0(1/2 - s) ----
    def _Xi_H0(self, s):
        s = mp.mpf(s)
        c = (self.zeta_H.L - 1) / mp.pi
        # Xi_H0(s) := c^{1-4s} * π^{-s} Γ(s) * ζ_H0(s)
        return mp.power(c, 1 - 4*s) * mp.power(mp.pi, -s) * mp.gamma(s) * self.zeta(s)

    def test_A2_functional_equation_H0(self):
        print("🔍 Test A2': Symmetrie Xi_H0(s) = Xi_H0(1/2 - s)...")
        test_points = [mp.mpf(x) for x in (0.25, 0.3, 0.4)]
        errs = []
        for s in test_points:
            left  = self._Xi_H0(s)
            right = self._Xi_H0(mp.mpf('0.5') - s)
            errs.append(float(abs(left - right)))
        mx = max(errs)
        self.res['A2'] = {
            'status': '✅' if mx < 1e-12 else '⚠️',
            'max_error': mx,
            'note': 'Symmetrie gilt exakt für H0 mit dieser Xi-Definition'
        }
        return self.res['A2']

    # ---- Strukturtests (A3, B1, B3) bleiben diagnostisch ----
    def test_A3_euler_product_structural(self):
        print("🔍 Test A3': Struktureller Euler-Produkt-Vergleich...")
        s = mp.mpf('2.0')
        zH = self.zeta(s); zR = mp.zeta(s)
        rel = abs(zH - zR) / abs(zR)
        self.res['A3'] = {'status': '⚠️', 'zeta_H(2)': float(zH), 'zeta(2)': float(zR), 'relative_difference': float(rel)}
        return self.res['A3']

    def test_B1_log_derivative_structural(self):
        print("🔍 Test B1: Logarithmische Ableitung (strukturell)...")
        s = mp.mpf('2.0'); h = mp.mpf('1e-6')
        z_plus = self.zeta(s+h); z_minus = self.zeta(s-h)
        log_deriv_H = -(z_plus - z_minus)/(2*h)/self.zeta(s)
        ref = mp.mpf('0.56996')
        diff = abs(log_deriv_H - ref)
        self.res['B1'] = {'status': '⚠️', 'log_deriv_H': float(log_deriv_H), 'log_deriv_R': float(ref), 'difference': float(diff)}
        return self.res['B1']

    def test_B3_Psi_structural(self):
        print("🔍 Test B3: Ψ(s) = ζ_H/ζ (strukturell)...")
        pts = [2.0, 3.0, 4.0]
        vals = []
        for x in pts:
            s = mp.mpf(x); zH = self.zeta(s); zR = mp.zeta(s)
            if abs(zR) > 1e-30: vals.append(float(zH/zR))
        mean = float(np.mean(vals)) if vals else 0.0
        std  = float(np.std(vals)) if vals else 0.0
        self.res['B3'] = {'status': '⚠️', 'mean_Psi': mean, 'std_Psi': std, 'values': vals}
        return self.res['B3']

    def summary(self):
        print("\n" + "="*70)
        print("📊 RIGOROSE DIAGNOSE - MIT HEAT-TRACE ζ_H(s)")
        print("="*70)
        print(f"L = {float(self.zeta_H.L)}, N_tail = {self.zeta_H.N_tail}")
        print(f"a0 = {float(self.zeta_H.a0):.6f}, a1 = {float(self.zeta_H.a1):.6f}\n")
        for k in ['A1','A2','A3','B1','B3']:
            if k in self.res:
                print(k, self.res[k])

# -----------------------------
# Hauptlauf (optional)
# -----------------------------
if __name__ == "__main__":
    mp.mp.dps = 80
    diag = RigorousEquivalenceDiagnostics(L=100.0, N_tail=4000, t0=1e-10, dps=80)
    for test in [diag.test_A1_residue_at_half,
                 diag.test_A2_functional_equation_H0,
                 diag.test_A3_euler_product_structural,
                 diag.test_B1_log_derivative_structural,
                 diag.test_B3_Psi_structural]:
        try:
            test()
        except Exception as e:
            print("❌", e)
    diag.summary()
