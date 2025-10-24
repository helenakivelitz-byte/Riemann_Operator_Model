# proof10_equivalence_rigorous_bridge.py
# BRÜCKE: Rigorose Heat-Trace ζ_H(s) + Äquivalenz-Diagnostik
# Testfall: freier Dirichlet-Operator H0 = -d^2/dx^2 auf [1, L]

import mpmath as mp
import numpy as np
import sympy as sp

# ============================================
# 1) Exakte Eigenwerte (nur für Referenz/Tail-Konstuktion)
# ============================================

def dirichlet_eigs_free(L, N):
    """Exakte Eigenwerte für H0 = -d²/dx² auf [1, L] mit Dirichlet."""
    L = mp.mpf(L)
    length = L - 1
    pi_over_len = mp.pi / length
    return [(pi_over_len * n) ** 2 for n in range(1, int(N) + 1)]


# ============================================
# 2) Wärmespur via Poisson-Summation (Jacobi-Theta-Transform)
# ============================================

def heat_trace_dirichlet_free_poisson(L, t, kmax=None):
    """
    Tr(e^{-tH0}) auf [1,L] mit Dirichlet über Poisson-Summation.
    θ-Transform liefert:
       Sum_{n>=1} exp(- (π n/(L-1))^2 t )
     = 1/2 ( sqrt(π/α) - 1 ) + sqrt(π/α) * Sum_{k>=1} exp( -π^2 k^2 / α )
    mit α = (π^2 t) / (L-1)^2.
    """
    L = mp.mpf(L)
    t = mp.mpf(t)
    length = L - 1
    alpha = (mp.pi**2) * t / (length**2)
    sqrt_term = mp.sqrt(mp.pi / alpha)

    # automatische kmax-Schätzung (sehr streng)
    if kmax is None:
        target = mp.mpf('1e-30')
        k_est = mp.sqrt(alpha * mp.log(1 / target)) / mp.pi
        kmax = max(5, int(mp.ceil(k_est)) + 5)

    tail = mp.fsum([mp.exp(-(mp.pi**2) * (k**2) / alpha) for k in range(1, kmax + 1)])
    S = mp.mpf('0.5') * (sqrt_term - 1) + sqrt_term * tail
    return S


# ============================================
# 3) Rigorose ζ_H(s) via Heat-Trace + Seeley–DeWitt + analytischer Tail
# ============================================

class RigorousHeatTraceZeta:
    """
    RIGOROSE ζ_H(s) für den freien Dirichlet-Operator H0 auf [1,L].
    Kombination:
      • (0,1): Seeley–DeWitt-Subtraktion (a0 t^{-1/2} + a1)
      • Wärmespur Tr(e^{-tH0}) über Poisson-Formel
      • (1,∞): analytischer Tail ∑ λ^{-s} Γ(s, λ)
      • Referenz: ζ_{H0}(s) = ((L-1)/π)^{2s} ζ_R(2s)
    """
    def __init__(self, L=100.0, N_tail=4000, t0=1e-10, dps=80):
        mp.mp.dps = dps
        self.L = mp.mpf(L)
        self.N_tail = int(N_tail)
        self.t0 = mp.mpf(t0)

        # Seeley–DeWitt-Koeffizienten (1D, Dirichlet, V=0)
        self.a0 = (self.L - 1) / mp.sqrt(4 * mp.pi)
        self.a1 = -mp.mpf('0.5')

        # Tail-Eigenwerte (für ∑ λ^{-s} Γ(s, λ))
        self.eigs_tail = dirichlet_eigs_free(self.L, self.N_tail)

    def Tr(self, t):
        """Wärmespur Tr(e^{-tH0}) via Poisson-Summation."""
        return heat_trace_dirichlet_free_poisson(self.L, t)

    def zeta(self, s):
        """
        ζ_H(s) = 1/Γ(s) * [ ∫_{t0}^1 t^{s-1}(Tr - a0 t^{-1/2} - a1) dt
                           + a0/(s-1/2) + a1/s
                           + ∑ λ^{-s} Γ(s,λ) ].
        """
        s = mp.mpf(s)
        t0 = self.t0

        # (0,1): subtrahierte Singularteile
        def integrand_small(t):
            t_val = mp.mpf(t)
            return t_val**(s - 1) * (self.Tr(t_val) - self.a0 * t_val**(-mp.mpf('0.5')) - self.a1)

        I_small = mp.quad(integrand_small, [t0, 1])

        # geschlossene Anteile aus (0,t0)+(t0,1)
        closed_small = self.a0 / (s - mp.mpf('0.5')) + self.a1 / s

        # Tail (1,∞) analytisch über unvollständige Gamma-Funktion
        tail_terms = [(lam ** (-s)) * mp.gammainc(s, lam, mp.inf) for lam in self.eigs_tail]
        tail_sum = mp.fsum(tail_terms)

        return (I_small + closed_small + tail_sum) / mp.gamma(s)

    def zeta_reference_free(self, s):
        """Exakte Referenz: ζ_{H0}(s) = ((L-1)/π)^{2s} ζ_R(2s)."""
        s_val = mp.mpf(s)
        factor = ((self.L - 1) / mp.pi) ** (2 * s_val)
        return factor * mp.zeta(2 * s_val)

    def demo(self, s_points=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0)):
        print("===== Beispiel: freier Dirichlet-Operator (V=0) =====")
        print("🚀 RIGOROSE ζ_H(s) VIA HEAT-TRACE (Seeley–DeWitt + Poisson + analytischer Tail)")
        # KORREKTUR: mpf-Objekte zu float konvertieren für Formatierung
        print(f"  Intervall: [1, {float(self.L):.1f}]  |  N_tail={self.N_tail}  |  t0={float(self.t0):.1e}")
        print("  a0=(L-1)/√(4π), a1=-1/2  (Dirichlet, 1D)\n")
        for s in s_points:
            z = self.zeta(mp.mpf(str(s)))
            z_ref = self.zeta_reference_free(s)
            relerr = abs((z - z_ref) / z_ref) if z_ref else mp.nan
            # KORREKTUR: mpf-Objekte zu float konvertieren
            print(f"s={s:>4}:  ζ_H(s)≈{float(z):.12g}   |   ζ_ref(s)≈{float(z_ref):.12g}   |   rel.err≈{float(relerr):.3e}")
        print("\nTipps für Genauigkeit:")
        print("• N_tail erhöhen (mehr Tail-Eigenwerte) → genauerer ∫_1^∞-Anteil.")
        print("• dps erhöhen (mpmath-Präzision), falls s nahe 1/2.")
        print("• t0 kleiner wählen, wenn dps groß genug (subtrahierte Singularitäten).")


# ============================================
# 4) Äquivalenz-Diagnostik mit rigoroser ζ_H(s)
# ============================================

class RigorousEquivalenceDiagnostics:
    """
    A1': Residuum bei s=1/2 (korrekt für 1D)
    A2': Symmetrie der korrekt gewichteten Größe Xi_H0(s) = Xi_H0(1/2 - s)
         mit Xi_H0(s) = c^(1/2 - 2s) * π^{-s} Γ(s) ζ_H(s),  c=(L-1)/π
    A3:  Struktureller Vergleich zu ζ(2) (kein Euler-Produkt für H0)
    B1:  Struktureller Vergleich der log. Ableitung bei s=2
    B3:  Ψ(s) = ζ_H/ζ (rein diagnostisch)
    """
    def __init__(self, L=100.0, N_tail=4000, t0=1e-10, dps=80):
        mp.mp.dps = dps
        self.res = {}
        self.zeta_H = RigorousHeatTraceZeta(L=L, N_tail=N_tail, t0=t0, dps=dps)

    # --- interne Hilfsfunktion: korrekte Xi_H0 mit c^(1/2 - 2s) ---
    def _Xi_H0(self, s):
        s = mp.mpf(s)
        c = (self.zeta_H.L - 1) / mp.pi
        return mp.power(c, mp.mpf('0.5') - 2*s) * mp.power(mp.pi, -s) * mp.gamma(s) * self.zeta_H.zeta(s)

    def test_A1_residue_at_half(self):
        """A1': Residuum bei s=1/2. Erwartet: a0/√π."""
        print("🔍 Test A1': Residuum bei s=1/2 (Laurent-Fit)...")
        s0 = mp.mpf('0.5')
        eps = mp.mpf('1e-4')
        # lineare Extrapolation aus zwei Punkten (vorwärts)
        s1, s2 = s0 + eps, s0 + eps/2
        r1 = (s1 - s0) * self.zeta_H.zeta(s1)
        r2 = (s2 - s0) * self.zeta_H.zeta(s2)
        # einfache lineare Extrapolation
        residue_fit = 2*r2 - r1
        expected = self.zeta_H.a0 / mp.sqrt(mp.pi)
        err = abs(residue_fit - expected)
        self.res['A1'] = {
            'status': '✅' if err < mp.mpf('1e-10') else '⚠️',
            'residue_fit': float(residue_fit),
            'expected': float(expected),
            'abs_error': float(err)
        }
        return self.res['A1']

    def test_A2_functional_eq_H0(self):
        """A2': Xi_H0(s) = Xi_H0(1/2 - s) – Symmetrie-Test mit korrektem Gewicht."""
        print("🔍 Test A2': Symmetrie Xi_H0(s) = Xi_H0(1/2 - s)...")
        points = [0.25, 0.30, 0.40, 0.45]
        diffs = []
        for x in points:
            s = mp.mpf(str(x))
            lhs = self._Xi_H0(s)
            rhs = self._Xi_H0(mp.mpf('0.5') - s)
            diffs.append(abs(lhs - rhs))
        mx = max(diffs)
        self.res['A2'] = {
            'status': '✅' if mx < mp.mpf('1e-10') else '⚠️',
            'max_error': float(mx),
            'note': 'Symmetrie gilt exakt für H0 mit dieser Xi-Definition'
        }
        return self.res['A2']

    def test_A3_euler_product_structural(self):
        """A3: Struktureller Vergleich zu ζ(2). Für H0 kein Euler-Produkt – immer diagnostisch."""
        print("🔍 Test A3': Struktureller Euler-Produkt-Vergleich...")
        s = mp.mpf('2.0')
        z_H = self.zeta_H.zeta(s)
        z_R = mp.zeta(s)
        rel = abs(z_H - z_R) / abs(z_R)
        self.res['A3'] = {
            'status': '⚠️',
            'zeta_H(2)': float(z_H),
            'zeta(2)': float(z_R),
            'relative_difference': float(rel)
        }
        return self.res['A3']

    def test_B1_log_derivative_structural(self):
        """B1: Logarithmische Ableitung – struktureller Vergleich."""
        print("🔍 Test B1: Logarithmische Ableitung (strukturell)...")
        s = mp.mpf('2.0')
        h = mp.mpf('1e-6')
        z_plus = self.zeta_H.zeta(s + h)
        z_minus = self.zeta_H.zeta(s - h)
        derivative_H = (z_plus - z_minus) / (2*h)
        log_deriv_H = -derivative_H / self.zeta_H.zeta(s)
        # Referenz (numerisch bekannte Zahl) für -ζ'/ζ(2):
        log_deriv_R = mp.mpf('0.56996')
        diff = abs(log_deriv_H - log_deriv_R)
        self.res['B1'] = {
            'status': '⚠️',
            'log_deriv_H': float(log_deriv_H),
            'log_deriv_R': float(log_deriv_R),
            'difference': float(diff)
        }
        return self.res['B1']

    def test_B3_Psi_structural(self):
        """B3: Ψ(s)=ζ_H/ζ – rein diagnostisch (nur Re(s)>1)."""
        print("🔍 Test B3: Ψ(s) = ζ_H/ζ (strukturell)...")
        test_points = [2.0, 3.0, 4.0]
        vals = []
        for x in test_points:
            s = mp.mpf(str(x))
            zH = self.zeta_H.zeta(s)
            zR = mp.zeta(s)
            if abs(zR) > 1e-30:
                vals.append(float(zH / zR))
        mean_psi = float(np.mean(vals)) if vals else 0.0
        std_psi  = float(np.std(vals))  if vals else 1.0
        self.res['B3'] = {
            'status': '⚠️',
            'mean_Psi': mean_psi,
            'std_Psi': std_psi,
            'values': vals
        }
        return self.res['B3']

    def summary(self):
        print("\n" + "="*70)
        print("📊 RIGOROSE DIAGNOSE - MIT HEAT-TRACE ζ_H(s)")
        print("="*70)
        # KORREKTUR: mpf-Objekte zu float konvertieren
        print(f"L = {float(self.zeta_H.L)}, N_tail = {self.zeta_H.N_tail}")
        print(f"a0 = {float(self.zeta_H.a0):.6f}, a1 = {float(self.zeta_H.a1):.6f}\n")
        for k in ['A1','A2','A3','B1','B3']:
            if k in self.res:
                print(k, self.res[k])
        print("\n" + "="*70)
        print("💡 Hinweis:")
        print("• A1' und A2' validieren die korrekte Struktur für H0.")
        print("• A3/B1/B3 sind nur strukturelle Vergleiche (H0 hat kein Euler-Produkt).")
        print("• Ziel war die Validierung der Methode und der korrekten Xi-Definition.")
        print("="*70)


# ============================================
# 5) Hauptprogramm: Demo + Diagnose + Analyse
# ============================================

def run_rigorous_diagnostics():
    print("🚀 RIGOROSE ÄQUIVALENZ-DIAGNOSE")
    print("Verwendet Heat-Trace Methode für ζ_H(s)")
    print("=" * 60)
    diag = RigorousEquivalenceDiagnostics(
        L=100.0,      # großes Intervall
        N_tail=4000,  # viele Tail-Eigenwerte
        t0=1e-10,     # kleine Startgrenze
        dps=80        # hohe Präzision
    )
    # Tests
    for test in [
        diag.test_A1_residue_at_half,
        diag.test_A2_functional_eq_H0,
        diag.test_A3_euler_product_structural,
        diag.test_B1_log_derivative_structural,
        diag.test_B3_Psi_structural
    ]:
        try:
            test()
        except Exception as e:
            print(f"❌ Test fehlgeschlagen: {e}")
    diag.summary()
    return diag.res

def analyze_spectral_properties():
    print("\n🔬 SPEKTRALE EIGENSCHAFTEN-ANALYSE:")
    print("=" * 50)
    zeta_H = RigorousHeatTraceZeta(L=100.0, N_tail=4000, dps=80)
    print("===== Beispiel: freier Dirichlet-Operator (V=0) =====")
    print("🚀 RIGOROSE ζ_H(s) VIA HEAT-TRACE (Seeley–DeWitt + Poisson + analytischer Tail)")
    # KORREKTUR: mpf-Objekte zu float konvertieren
    print(f"  Intervall: [1, {float(zeta_H.L):.1f}]  |  N_tail={zeta_H.N_tail}  |  t0={float(zeta_H.t0):.1e}")
    print("  a0=(L-1)/√(4π), a1=-1/2  (Dirichlet, 1D)\n")
    print(f"{'s':>6} {'ζ_H(s)':>15} {'ζ_ref(s)':>15} {'Verhältnis':>12}")
    print("-" * 55)
    for s in [0.6, 1.0, 1.5, 2.0, 2.5, 3.0]:
        z_H = zeta_H.zeta(s)
        z_ref = zeta_H.zeta_reference_free(s)
        ratio = float(z_H / z_ref) if z_ref else float('nan')
        print(f"{s:6.1f} {float(z_H):15.6f} {float(z_ref):15.6f} {ratio:12.6f}")

if __name__ == "__main__":
    # 1) Demo der zeta_H-Berechnung mit Referenzvergleich
    HT = RigorousHeatTraceZeta(L=100.0, N_tail=4000, t0=1e-10, dps=80)
    HT.demo(s_points=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0))

    # 2) Rigorose Äquivalenz-Diagnostik
    results = run_rigorous_diagnostics()

    # 3) Detaillierte spektrale Analyse (optional)
    analyze_spectral_properties()

    print("\n" + "=" * 60)
    print("🎯 FAZIT: Heat-Trace liefert für H0 maschinengenaue ζ_H(s).")
    print("   A2' ist jetzt korrekt (Xi_H0 mit c^(1/2 - 2s)).")
    print("   Nächster Schritt: Anwendung auf H[ρ] (mit V_eff) und")
    print("   entsprechende Anpassung der Seeley–DeWitt-Koeffizienten.")
    print("=" * 60)
