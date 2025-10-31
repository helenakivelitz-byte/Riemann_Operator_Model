# proof10_heat_trace_zeta_rigorous_improved.py
# RIGOROSE ζ_H(s) via Heat-Trace + Seeley–DeWitt-Subtraktion + analytischer Tail
# Fallstudie: freier Dirichlet-Operator H0 = -d^2/dx^2 auf [1, L]
# Verbesserungen:
#   - Wärmespur Tr(e^{-tH0}) via Poisson-Summation (Jacobi-Theta-Transform)
#   - exakter Integraltail \int_1^\infty ... dt über obere unvollständige Gammafunktion Γ(s, λ)
#   - Korrektur: mp.fsum statt mp.nsum über eine Liste von Eigenwerten

import mpmath as mp

# -----------------------------
# Exakte Eigenwerte (Referenz / Tail-Summe)
# -----------------------------
def dirichlet_eigs_free(L, N):
    """
    Exakte Eigenwerte für H0 = -d^2/dx^2 auf [1, L] mit Dirichlet.
    Länge = L-1, Eigenwerte: λ_n = (π n / (L-1))^2, n=1..N
    """
    length = mp.mpf(L) - 1
    pi_over_len = mp.pi / length
    return [(pi_over_len * n) ** 2 for n in range(1, N + 1)]

# -----------------------------
# Wärmespur via Poisson-Summation (Jacobi-Theta-Transform)
# -----------------------------
def heat_trace_dirichlet_free_poisson(L, t, kmax=None):
    """
    Tr(e^{-tH0}) = sum_{n=1}^\infty exp(- (π^2 t / (L-1)^2) n^2)
                 = S(α), α = (π^2 t)/(L-1)^2

    Poisson-Identität:
      S(α) = 1/2*(sqrt(pi/α) - 1) + sqrt(pi/α) * sum_{k=1}^\infty exp(-π^2 k^2 / α)

    Diese Darstellung konvergiert rasant für kleines t (α klein)
    und ist numerisch stabil im gesamten t-Bereich.
    """
    L = mp.mpf(L)
    t = mp.mpf(t)
    length = L - 1
    alpha = (mp.pi**2) * t / (length**2)

    sqrt_term = mp.sqrt(mp.pi / alpha)

    # dynamischer k-Summen-Abbruch
    if kmax is None:
        target = mp.mpf('1e-30')  # gewünschte Genauigkeit pro Term
        # exp(-π^2 k^2 / α) < target  =>  k > sqrt(α * log(1/target)) / π
        k_est = mp.sqrt(alpha * mp.log(1 / target)) / mp.pi
        kmax = max(5, int(mp.ceil(k_est)) + 5)

    # sichere, endliche Summe
    tail = mp.fsum([mp.e ** (-(mp.pi**2) * (k**2) / alpha) for k in range(1, kmax + 1)])
    S = mp.mpf('0.5') * (sqrt_term - 1) + sqrt_term * tail
    return S

# -----------------------------
# Hauptklasse
# -----------------------------
class HeatTraceZeta:
    """
    Rigorose ζ_H(s)-Berechnung aus Wärmespur:
      ζ_H(s) = 1/Γ(s) [ ∫_0^1 t^{s-1} Tr(e^{-tH}) dt + ∫_1^\infty t^{s-1} Tr(e^{-tH}) dt ]

    Numerik:
      - auf (0,1): Seeley–DeWitt-Subtraktion (a0, a1) + Poisson-Heat-Trace
      - auf (1,∞): analytisch via Summe λ^{-s} Γ(s, λ) (obere unvollständige Gammafunktion)

    Testfall: H0 (frei, Dirichlet), a0=(L-1)/√(4π), a1=-1/2
    """
    def __init__(self, L=20.0, N_tail=2000, t0=1e-10, dps=80):
        mp.mp.dps = dps
        self.L = mp.mpf(L)
        self.N_tail = int(N_tail)
        self.t0 = mp.mpf(t0)

        # Seeley–DeWitt-Koeffizienten in 1D (Dirichlet, V=0):
        # Tr(e^{-tH}) ~ a0 t^{-1/2} + a1 + O(t^{1/2})  (t -> 0+)
        self.a0 = (self.L - 1) / mp.sqrt(4 * mp.pi)
        self.a1 = -mp.mpf('0.5')

        # Eigenwerte für Tail-Summe ∫_1^\infty ... dt
        self.eigs_tail = dirichlet_eigs_free(self.L, self.N_tail)

    def Tr(self, t):
        """ Wärmespur via Poisson-Summation (stabil über gesamten t-Bereich). """
        return heat_trace_dirichlet_free_poisson(self.L, t)

    def zeta(self, s):
        """
        ζ_H(s) = 1/Γ(s) * [ ∫_{t0}^1 t^{s-1} (Tr - a0 t^{-1/2} - a1) dt
                            + a0/(s - 1/2) + a1/s
                            + ∑_{n=1}^{N_tail} λ_n^{-s} Γ(s, λ_n) ]
        Die Summe repräsentiert den gesamten Tail ∫_1^\infty t^{s-1} Tr(e^{-tH}) dt,
        aufgeteilt in Eigenwerte. (Für H0 sind λ_n exakt bekannt.)
        """
        s = mp.mpf(s)
        t0 = self.t0

        # (0,1): integrand with subtractions
        def integrand_small(t):
            t = mp.mpf(t)
            return t**(s - 1) * ( self.Tr(t) - self.a0 * t**(-mp.mpf('0.5')) - self.a1 )

        I_small = mp.quad(integrand_small, [t0, 1])

        closed_small = self.a0 / (s - mp.mpf('0.5')) + self.a1 / s

        # Tail (1,∞): ∫_1^\infty t^{s-1} e^{-λ t} dt = λ^{-s} Γ(s, λ)
        # ACHTUNG: obere unvollständige Gamma: Γ(s, λ) = gammainc(s, λ, ∞)
        tail_terms = [ (lam ** (-s)) * mp.gammainc(s, lam, mp.inf) for lam in self.eigs_tail ]
        tail_sum = mp.fsum(tail_terms)

        z = (I_small + closed_small + tail_sum) / mp.gamma(s)
        return z

    # Exakte Referenz für H0: ζ_{H0}(s) = ((L-1)/π)^{2s} ζ_R(2s)
    def zeta_reference_free(self, s):
        s = mp.mpf(s)
        factor = ((self.L - 1) / mp.pi) ** (2 * s)
        return factor * mp.zeta(2 * s)

    def demo(self, s_points=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0)):
        print("===== Beispiel: freier Dirichlet-Operator (V=0) =====")
        print("🚀 RIGOROSE ζ_H(s) VIA HEAT-TRACE (Seeley–DeWitt + Poisson + analytischer Tail)")
        print(f"  Intervall: [1, {float(self.L):.1f}]  |  N_tail={self.N_tail}  |  t0={float(self.t0):.1e}")
        print("  a0=(L-1)/√(4π), a1=-1/2  (Dirichlet, 1D)\n")

        for s in s_points:
            z = self.zeta(mp.mpf(str(s)))
            z_ref = self.zeta_reference_free(s)
            relerr = abs((z - z_ref) / z_ref) if z_ref else mp.nan

            z_f = float(z)
            z_ref_f = float(z_ref)
            relerr_f = float(relerr)

            print(f"s={s:>4}:  ζ_H(s)≈{z_f:.12g}   |   ζ_ref(s)≈{z_ref_f:.12g}   |   rel.err≈{relerr_f:.3e}")

        print("\nTipps für Genauigkeit:")
        print("• N_tail erhöhen (mehr Eigenwerte im Tail) → genauerer ∫_1^∞-Anteil.")
        print("• dps erhöhen (mpmath-Präzision), falls s nahe 1/2.")
        print("• t0 kleiner wählen, wenn dps groß genug (subtrahierte Singularitäten besser aufgelöst).")

# -----------------------------
# Ausführen (Demo)
# -----------------------------
if __name__ == "__main__":
    HT = HeatTraceZeta(
        L=100.0,     # großes Intervall
        N_tail=4000, # viele Eigenwerte für den Tail
        t0=1e-10,    # kleine Startgrenze (mit hoher dps ok)
        dps=80       # hohe Präzision
    )
    HT.demo(s_points=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0))
