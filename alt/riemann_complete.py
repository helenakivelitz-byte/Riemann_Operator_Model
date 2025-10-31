# riemann_complete.py - KORRIGIERTE VERSION
"""
VOLLSTÄNDIGES RIEMANN OPERATOR MODELL
Alles in einer Datei - Einfache Ausführung
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigs
from scipy import stats
import mpmath as mp

print("🚀 Starte komplettes Riemann Operator Modell...")

class RiemannOperator:
    def __init__(self, L=12.0, n_points=200):
        self.L = L
        self.n_points = n_points
        self.x = np.linspace(1.0, L, n_points)
        self.h = (L - 1.0) / (n_points - 1)
        
        # Parameter aus dem Paper
        self.epsilon = 0.1      # Regularisierung
        self.gamma = 0.35       # Gravitative Kopplung
        self.sigma = 0.5        # Wechselwirkungsreichweite
        self.eta = 0.25         # Austauschkoeffizient
        self.rho_min = 1e-6     # Dichte-Minimum
        
    def build_laplacian(self):
        n = self.n_points - 2
        main_diag = 2.0 * np.ones(n) / self.h**2
        off_diag = -1.0 * np.ones(n-1) / self.h**2
        return sparse.diags([main_diag, off_diag, off_diag], [0, 1, -1])
    
    def V_prim(self, x):
        """Primzahl-Potential aus Paper"""
        primes = [2, 3, 5, 7, 11]  # p ≤ 12
        V = 0.0
        for p in primes:
            V += (1/p) * np.exp(-(x - np.log(p))**2 / (2 * self.epsilon**2))
        return V
    
    def V_grav(self, rho):
        """Gravitatives Potential"""
        x_grid, y_grid = np.meshgrid(self.x, self.x, indexing='ij')
        r = np.abs(x_grid - y_grid)
        K = 1.0 / np.sqrt(1 + r**2) * np.exp(-r**2 / (2 * self.sigma**2))
        
        weights = np.ones_like(self.x) * self.h
        weights[0] = weights[-1] = self.h / 2.0
        
        return self.gamma * np.dot(K, rho * weights)
    
    def V_exch(self, rho):
        """Austausch-Potential"""
        rho_safe = np.maximum(rho, self.rho_min)
        return -self.eta * rho_safe**(1/3)
    
    def build_hamiltonian(self, rho):
        H0 = self.build_laplacian()
        
        # Potentials auf inneren Punkten
        V_prim_inner = self.V_prim(self.x[1:-1])
        V_grav_inner = self.V_grav(rho)[1:-1] 
        V_exch_inner = self.V_exch(rho)[1:-1]
        
        V_total = V_prim_inner + V_grav_inner + V_exch_inner
        V_diag = sparse.diags(V_total, 0)
        
        return H0 + V_diag

class SCFSolver:
    def __init__(self, operator, max_iter=30, tolerance=1e-5, beta=0.3, N_states=8):
        self.operator = operator
        self.max_iter = max_iter
        self.tolerance = tolerance
        self.beta = beta
        self.N_states = N_states
        
    def boltzmann_weights(self, energies, T=0.75):
        energies_shifted = energies - np.min(energies)
        weights = np.exp(-energies_shifted / T)
        return weights / np.sum(weights)
    
    def compute_density(self, eigenvectors, weights):
        """KORRIGIERT: Dichte aus gewichteten Eigenfunktionen berechnen"""
        # eigenvectors shape: (n_points, N_states)
        # weights shape: (N_states,)
        density = np.zeros_like(eigenvectors[:, 0])
        for i in range(len(weights)):
            density += weights[i] * np.abs(eigenvectors[:, i])**2
        return density
    
    def solve(self):
        n_points = self.operator.n_points
        rho = np.ones(n_points) / (self.operator.L - 1.0)
        
        print("🔄 Starte SCF-Iteration...")
        
        for iteration in range(self.max_iter):
            # Hamiltonian für aktuelle Dichte
            H = self.operator.build_hamiltonian(rho)
            
            try:
                # Eigenwerte berechnen
                energies, vectors = eigs(H, k=self.N_states, which='SR')
                energies = np.real(energies)
                vectors = np.real(vectors)
                
                # Vollständige Eigenfunktionen (mit Rand=0)
                vectors_full = np.zeros((n_points, self.N_states))
                vectors_full[1:-1, :] = vectors
                
            except Exception as e:
                print(f"⚠️ Iteration {iteration}: {e}")
                continue
            
            # Neue Dichte berechnen - KORRIGIERT
            weights = self.boltzmann_weights(energies)
            rho_new = self.compute_density(vectors_full, weights)
            
            # Dämpfung und Renormierung
            rho_mixed = self.beta * rho + (1 - self.beta) * rho_new
            rho_mixed = rho_mixed / np.trapz(rho_mixed, self.operator.x)
            
            # Konvergenz checken
            diff = np.linalg.norm(rho_mixed - rho)
            if iteration % 5 == 0:
                print(f"   Iteration {iteration+1}: ||Δρ|| = {diff:.2e}")
            
            if diff < self.tolerance:
                print(f"✅ SCF konvergiert nach {iteration+1} Iterationen!")
                break
                
            rho = rho_mixed
        
        return rho, energies, vectors_full

class SpectralAnalyzer:
    def __init__(self, eigenvalues):
        self.eigenvalues = np.sort(eigenvalues)
        
    def level_spacing_statistics(self):
        spacings = np.diff(self.eigenvalues)
        if len(spacings) < 3:
            return {'error': 'Zu wenige Eigenwerte'}
            
        normalized_spacings = spacings / np.mean(spacings)
        
        # GOE vs Poisson
        def goe_cdf(s):
            return 1 - np.exp(-np.pi * s**2 / 4)
            
        ks_goe = stats.kstest(normalized_spacings, goe_cdf).statistic
        ks_poisson = stats.kstest(normalized_spacings, 'expon').statistic
        
        # Gap-Korrelation
        gap_correlation = np.corrcoef(spacings[:-1], spacings[1:])[0,1] if len(spacings) > 2 else 0
        
        return {
            'mean_spacing': float(np.mean(spacings)),
            'ks_goe': float(ks_goe),
            'ks_poisson': float(ks_poisson),
            'gap_correlation': float(gap_correlation)
        }
    
    def prime_correlation(self):
        primes = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29])
        
        if len(self.eigenvalues) >= len(primes):
            ev_subset = self.eigenvalues[:len(primes)]
            correlation = np.corrcoef(ev_subset, primes)[0,1]
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(ev_subset, primes)
            
            return {
                'correlation_r2': float(r_value**2),
                'linear_fit_slope': float(slope),
                'n_primes_compared': len(primes)
            }
        else:
            return {'error': 'Nicht genug Eigenwerte'}

def main():
    print("=" * 60)
    print("RIEMANN OPERATOR MODELL - KOMPLETTE VERSION")
    print("=" * 60)
    
    # 1. Operator erstellen
    print("1. Erstelle Riemann Operator...")
    operator = RiemannOperator(L=12.0, n_points=200)
    
    # 2. SCF-Lösung berechnen
    print("2. Berechne self-consistent density...")
    solver = SCFSolver(operator, max_iter=30, N_states=8)
    rho_star, eigenvalues, eigenvectors = solver.solve()
    
    # 3. Spektrale Analyse
    print("3. Analysiere Eigenwert-Spektrum...")
    analyzer = SpectralAnalyzer(eigenvalues)
    
    level_stats = analyzer.level_spacing_statistics()
    prime_stats = analyzer.prime_correlation()
    
    # 4. Ergebnisse anzeigen
    print("\n" + "=" * 60)
    print("🎯 ERGEBNISSE:")
    print("=" * 60)
    print(f"• Eigenwerte berechnet: {len(eigenvalues)}")
    print(f"• Eigenwert-Bereich: {np.min(eigenvalues):.3f} - {np.max(eigenvalues):.3f}")
    
    if 'error' not in level_stats:
        print(f"• Level-Statistik: KS(GOE) = {level_stats['ks_goe']:.3f}")
        print(f"• Gap-Korrelation: {level_stats['gap_correlation']:.3f}")
    
    if 'error' not in prime_stats:
        print(f"• Primzahl-Korrelation: R² = {prime_stats['correlation_r2']:.4f}")
    
    print(f"• Dichte-Integral: {np.trapz(rho_star, operator.x):.6f} (sollte ≈1.0)")
    
    print("\n✅ MODELL ERFOLGREICH AUSGEFÜHRT!")
    
    # Eigenwerte speichern für spätere Analyse
    np.savetxt('eigenvalues.txt', eigenvalues)
    print("💾 Eigenwerte gespeichert in 'eigenvalues.txt'")

if __name__ == "__main__":
    main()
