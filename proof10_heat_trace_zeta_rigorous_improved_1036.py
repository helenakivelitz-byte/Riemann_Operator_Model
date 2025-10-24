# proof10_heat_trace_zeta_rigorous_improved.py
# VERBESSERTE RIGOROSE ζ_H(s) via Heat-Trace + Seeley–DeWitt-Subtraktion

import mpmath as mp
import time

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
    return [ (pi_over_len * n)**2 for n in range(1, M+1) ]

def heat_trace_from_eigs(eigs, t):
    """
    Tr(e^{-tH}) = Sum_j exp(-t * λ_j)
    Optimierte Version für große M.
    """
    t_mp = mp.mpf(t)
    total = mp.mpf(0)
    # Für sehr kleine t, breche früh ab wenn Terme vernachlässigbar werden
    for i, lam in enumerate(eigs):
        term = mp.exp(-t_mp * lam)
        total += term
        # Frühzeitiges Abbrechen für kleine Terme
        if term < 1e-20 and i > 50:
            break
    return total

# -----------------------------
# Hauptklasse
# -----------------------------

class ImprovedHeatTraceZeta:
    """
    VERBESSERTE rigorose ζ_H(s)-Berechnung
    """
    def __init__(self, L=100.0, M=2000, T_large=50.0, t0=1e-12, dps=80):
        mp.mp.dps = dps
        self.L = mp.mpf(L)
        self.M = M
        self.T_large = mp.mpf(T_large)
        self.t0 = mp.mpf(t0)

        print("🔧 Berechne Eigenwerte...")
        start_time = time.time()
        self.eigs = dirichlet_eigs_free(self.L, self.M)
        eigen_time = time.time() - start_time
        print(f"✅ {self.M} Eigenwerte berechnet in {eigen_time:.2f}s")

        # Seeley–DeWitt-Koeffizienten
        self.a0 = (self.L - 1) / mp.sqrt(4 * mp.pi)
        self.a1 = -mp.mpf('0.5')
        
        print(f"📊 Parameter: L={float(self.L)}, M={M}, T_large={float(self.T_large)}")
        print(f"   a0={float(self.a0):.6f}, a1={float(self.a1):.6f}")
        print(f"   λ_min={float(self.eigs[0]):.8f}, λ_max={float(self.eigs[-1]):.2f}")

    def Tr(self, t):
        """Optimierte Wärmespur-Berechnung"""
        return heat_trace_from_eigs(self.eigs, t)

    def zeta_improved(self, s):
        """
        VERBESSERTE ζ_H(s) Berechnung mit besserer Numerik
        """
        s = mp.mpf(s)
        
        print(f"🔍 Berechne ζ_H({float(s)})...")
        start_time = time.time()

        # 1. Integral (t0, 1) mit Subtraktion
        def integrand_small(t):
            t_val = mp.mpf(t)
            Tr_val = self.Tr(t_val)
            asympt = self.a0 * t_val**(-0.5) + self.a1
            return t_val**(s-1) * (Tr_val - asympt)

        # 2. Integral (1, T_large)
        def integrand_large(t):
            t_val = mp.mpf(t)
            return t_val**(s-1) * self.Tr(t_val)

        try:
            I1 = mp.quad(integrand_small, [self.t0, 1], error=True)
            I2 = mp.quad(integrand_large, [1, self.T_large], error=True)
            
            # Geschlossene Terme
            closed_small = self.a0 / (s - 0.5) + self.a1 / s
            
            z = (I1 + I2 + closed_small) / mp.gamma(s)
            
            calc_time = time.time() - start_time
            print(f"✅ ζ_H({float(s)}) = {float(z):.8f} (Berechnungszeit: {calc_time:.2f}s)")
            
            return z
            
        except Exception as e:
            print(f"❌ Fehler bei ζ_H({float(s)}): {e}")
            return mp.mpf('0')

    # Referenz-Formel für H0
    def zeta_reference_free(self, s):
        s_val = mp.mpf(s)
        factor = ((self.L - 1) / mp.pi)**(2*s_val)
        return factor * mp.zeta(2*s_val)

    def analyze_convergence(self, s_test=2.0):
        """
        Analysiere Konvergenz bezüglich M und T_large
        """
        print("\n🔬 KONVERGENZANALYSE:")
        print("=" * 50)
        
        s = mp.mpf(s_test)
        z_ref = self.zeta_reference_free(s)
        
        # Test verschiedene M Werte
        M_values = [500, 1000, 2000, 3000]
        print("📈 Konvergenz bezüglich M:")
        for M_test in M_values:
            test_eigs = dirichlet_eigs_free(self.L, M_test)
            original_eigs = self.eigs
            self.eigs = test_eigs
            
            z_test = self.zeta_improved(s)
            error = abs((z_test - z_ref) / z_ref) if abs(z_ref) > 1e-15 else float('inf')
            print(f"   M={M_test:4}: ζ_H = {float(z_test):.8f}, rel. Fehler = {float(error):.2e}")
            
            self.eigs = original_eigs

    def demo_comprehensive(self, s_points=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0)):
        """
        Umfassende Demo mit Genauigkeitsanalyse
        """
        print("\n" + "="*70)
        print("🎯 RIGOROSE ζ_H(s) - VERBESSERTE VERSION")
        print("="*70)
        print(f"Intervall: [1, {float(self.L)}] | Eigenwerte: {self.M} | T_large: {float(self.T_large)}")
        print(f"Seeley–DeWitt: a₀ = {float(self.a0):.6f}, a₁ = {float(self.a1):.6f}")
        print()

        results = []
        
        for s_val in s_points:
            s = mp.mpf(str(s_val))
            
            # Unsere Berechnung
            z_H = self.zeta_improved(s)
            
            # Referenz (exakt für H₀)
            z_ref = self.zeta_reference_free(s)
            
            # Fehleranalyse
            abs_error = abs(z_H - z_ref)
            rel_error = abs_error / abs(z_ref) if abs(z_ref) > 1e-15 else float('inf')
            
            results.append({
                's': s_val,
                'zeta_H': float(z_H),
                'zeta_ref': float(z_ref),
                'abs_error': float(abs_error),
                'rel_error': float(rel_error)
            })
            
            print(f"s={s_val:>4}: ζ_H = {float(z_H):12.6f} | ζ_ref = {float(z_ref):12.6f} | rel. Fehler = {float(rel_error):8.2e}")

        # Zusammenfassung
        print("\n" + "="*70)
        print("📊 ZUSAMMENFASSUNG:")
        avg_rel_error = np.mean([r['rel_error'] for r in results])
        max_rel_error = np.max([r['rel_error'] for r in results])
        
        print(f"Durchschnittlicher relativer Fehler: {avg_rel_error:.2e}")
        print(f"Maximaler relativer Fehler: {max_rel_error:.2e}")
        
        if avg_rel_error < 0.01:
            print("✅ EXZELLENTE GENAUIGKEIT")
        elif avg_rel_error < 0.1:
            print("⚠️  GUTE GENAUIGKEIT")
        else:
            print("❌ GENAUIGKEIT BENÖTIGT VERBESSERUNG")
            
        print("="*70)

    def test_pole_structure(self):
        """
        Teste die Polstruktur bei s = 1/2, -1/2, -3/2, ...
        """
        print("\n🔍 TEST DER POLSTRUKTUR:")
        print("=" * 50)
        
        pole_points = [0.5, -0.5, -1.5]
        
        for s_val in pole_points:
            s = mp.mpf(str(s_val))
            print(f"\nTest bei s = {s_val}:")
            
            try:
                z = self.zeta_improved(s + 0.001)  # Vermeide exakten Pol
                print(f"   ζ_H({float(s+0.001):.3f}) = {float(z):.6f}")
                
                # Teste Residuum
                if s_val == 0.5:
                    residue_expected = self.a0 / mp.sqrt(mp.pi)
                    print(f"   Erwartetes Residuum bei s=1/2: {float(residue_expected):.6f}")
                    
            except Exception as e:
                print(f"   ❌ Fehler: {e}")

# -----------------------------
# Ausführen (Verbesserte Demo)
# -----------------------------

if __name__ == "__main__":
    print("🚀 STARTE VERBESSERTE ζ_H(s) BERECHNUNG")
    print("Bitte warten... Dies kann einige Minuten dauern.\n")
    
    # VERBESSERTE PARAMETER:
    # - Größeres L für bessere L→∞ Approximation
    # - Mehr Eigenwerte M für genauere Heat-Trace
    # - Höhere Präzision dps
    # - Größeres T_large für bessere Integration
    
    HT = ImprovedHeatTraceZeta(
        L=100.0,      # Größeres Intervall
        M=2000,       # Mehr Eigenwerte  
        T_large=50.0, # Längere Integration
        t0=1e-12,     # Feinere Auflösung nahe 0
        dps=80        # Höhere Präzision
    )
    
    # Umfassende Demo
    HT.demo_comprehensive(s_points=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0))
    
    # Optional: Konvergenzanalyse
    # HT.analyze_convergence(s_test=2.0)
    
    # Optional: Polstruktur-Test
    HT.test_pole_structure()
    
    print("\n💡 TIPPS FÜR NOCH BESSERE GENAUIGKEIT:")
    print("   - Erhöhe M auf 5000+ für sehr kleine t")
    print("   - Erhöhe T_large auf 100+ für große s")
    print("   - Verwende dps=100+ für höchste Präzision")
    print("   - Für L→∞: verwende L=1000+")
