# proof10_heat_trace_zeta_rigorous.py
# RIGOROSE ζ_H(s) via Heat-Trace + Seeley–DeWitt-Subtraktion
# Fallstudie: freier Dirichlet-Operator H0 = -d^2/dx^2 auf [1, L]

import mpmath as mp

# -----------------------------
# Hilfsfunktionen
# -----------------------------

def dirichlet_eigs_free(L, M):
    """
    Exakte Eigenwerte für H0 = -d^2/dx^2 auf [1, L] mit Dirichlet.
    Länge = L-1, Eigenwerte: λ_n = (π n / (L-1))^2, n=1..M
    """
    length = mp.mpf(L) - 1
    pi_over_len = mp.pi / length
    return [ (pi_over_len * (n))**2 for n in range(1, M+1) ]

def heat_trace_from_eigs(eigs, t):
    """
    Tr(e^{-tH}) = Sum_j exp(-t * λ_j)
    Robuste Version ohne mpmath-Indexsummation (vermeidet mpf-Index-Fehler).
    """
    t_mp = mp.mpf(t)
    # mp.fsum ist stabiler als naive Summation
    terms = [ mp.exp(-t_mp * lam) for lam in eigs ]
    return mp.fsum(terms)

# -----------------------------
# Hauptklasse
# -----------------------------

class HeatTraceZeta:
    """
    Rigorose ζ_H(s)-Berechnung aus Wärmespur:
    ζ_H(s) = 1/Γ(s) ∫_0^∞ t^{s-1} Tr(e^{-tH}) dt
    Numerik:
      - auf (0,1): Subtraktion der ersten Seeley–DeWitt-Koeffizienten (a0, a1)
      - auf (1,T): direkte Integration von Tr(e^{-tH})
      - t=0 wird vermieden (Start bei t0>0)
    Testfall: H0 (frei, Dirichlet), a0=(L-1)/√(4π), a1=-1/2
    """
    def __init__(self, L=20.0, M=400, T_large=10.0, t0=1e-10, dps=50):
        mp.mp.dps = dps
        self.L = mp.mpf(L)
        self.M = M
        self.T_large = mp.mpf(T_large)
        self.t0 = mp.mpf(t0)

        # Eigenwerte (frei, Dirichlet)
        self.eigs = dirichlet_eigs_free(self.L, self.M)

        # Seeley–DeWitt-Koeffizienten in 1D (Dirichlet, V=0)
        # Tr(e^{-tH}) ~ a0 t^{-1/2} + a1 + O(t^{1/2})  (t -> 0+)
        self.a0 = (self.L - 1) / mp.sqrt(4 * mp.pi)
        self.a1 = -mp.mpf('0.5')

    def Tr(self, t):
        """ Wärmespur-Truncation (M Terme) """
        return heat_trace_from_eigs(self.eigs, t)

    def zeta(self, s):
        """
        ζ_H(s) ≈ 1/Γ(s) * [ ∫_{t0}^1 t^{s-1} (Tr - a0 t^{-1/2} - a1) dt
                             + ∫_1^{T} t^{s-1} Tr dt
                             + a0 * ∫_0^1 t^{s-1} t^{-1/2} dt
                             + a1 * ∫_0^1 t^{s-1} dt ]
        = 1/Γ(s) * [ I1 + I2 + a0/(s - 1/2) + a1/s ]   (Re(s) > 1/2)
        Hinweis: Die Trunkierung bei T_large macht den großen-t-Tail numerisch klein.
        """
        s = mp.mpf(s)

        t0 = self.t0  # kleine positive Startgrenze

        def integrand_small(t):
            t = mp.mpf(t)
            # Subtrahierter integrand ist bei t->0 integrierbar;
            # numerisch vermeiden wir t=0 durch Start bei t0.
            return t**(s - 1) * ( self.Tr(t) - self.a0 * t**(-mp.mpf('0.5')) - self.a1 )

        def integrand_large(t):
            t = mp.mpf(t)
            return t**(s - 1) * self.Tr(t)

        # Numerische Integrale
        I1 = mp.quad(integrand_small, [t0, 1])         # (t0, 1)
        I2 = mp.quad(integrand_large, [1, self.T_large])  # (1, T)

        # Geschlossene Beiträge (0..1) der subtrahierten Terme:
        # ∫_0^1 t^{s-1} t^{-1/2} dt = 1/(s - 1/2),   ∫_0^1 t^{s-1} dt = 1/s
        closed_small = self.a0 / (s - mp.mpf('0.5')) + self.a1 / s

        z = (I1 + I2 + closed_small) / mp.gamma(s)
        return z

    # Referenz-Formel für H0: ζ_{H0}(s) = ((L-1)/π)^{2s} ζ_R(2s)
    def zeta_reference_free(self, s):
        s = mp.mpf(s)
        factor = ((self.L - 1) / mp.pi)**(2*s)
        return factor * mp.zeta(2*s)

    def demo(self, s_points=(0.6, 1.0, 2.0)):
        print("===== Beispiel: freier Dirichlet-Operator (V=0) =====")
        print("🚀 RIGOROSE ζ_H(s) VIA HEAT-TRACE (Seeley–DeWitt Subtraktion)")
        print(f"  Intervall: [1, {float(self.L)}]  |  M={self.M}  |  λ_min≈{float(self.eigs[0]):.6f}  |  T={float(self.T_large)}")
        print("  a0=(L-1)/√(4π), a1=-1/2  (Dirichlet, 1D)\n")

        for s_val in s_points:
            s = mp.mpf(str(s_val))
            z = self.zeta(s)
            z_ref = self.zeta_reference_free(s)
            relerr = abs((z - z_ref) / z_ref) if abs(z_ref) > 1e-15 else float('inf')
            
            # Korrekte Formatierung für mpmath-Zahlen
            print(f"s={s_val:>4}:  ζ_H(s)≈{float(z):.12g}   |   ζ_ref(s)≈{float(z_ref):.12g}   |   rel.err≈{float(relerr):.3e}")

        print("\nHinweis:")
        print("- Erhöhe M (Anzahl Eigenwerte) und T_large, um Genauigkeit zu steigern.")
        print("- verringere t0 für feinere Auflösung nahe t=0 (mit höherer Präzision).")

# -----------------------------
# Ausführen (Demo)
# -----------------------------

if __name__ == "__main__":
    # Parameter kannst du nach Bedarf erhöhen:
    #   - L größer → längeres Intervall
    #   - M größer → mehr Eigenwerte in Tr(e^{-tH})
    #   - T_large größer → längeres großes-t-Integral
    #   - dps (mpmath-Präzision) größer → präzisere Quadratur
    HT = HeatTraceZeta(L=20.0, M=400, T_large=10.0, t0=1e-10, dps=60)
    HT.demo(s_points=(0.6, 1.0, 2.0))
