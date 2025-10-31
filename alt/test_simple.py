# test_simple.py
print("🚀 Starte einfachen Test des Riemann Operators...")
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigs

class SimpleRiemannOperator:
    def __init__(self, L=12.0, n_points=100):
        self.L = L
        self.n_points = n_points
        self.x = np.linspace(1.0, L, n_points)
        self.h = (L - 1.0) / (n_points - 1)
        
    def build_laplacian(self):
        n = self.n_points - 2
        main_diag = 2.0 * np.ones(n) / self.h**2
        off_diag = -1.0 * np.ones(n-1) / self.h**2
        return sparse.diags([main_diag, off_diag, off_diag], [0, 1, -1])
    
    def build_hamiltonian(self, rho):
        H0 = self.build_laplacian()
        # Einfaches Test-Potential
        V_test = 0.1 * np.sin(self.x[1:-1])  # Beispiel-Potential
        V_diag = sparse.diags(V_test, 0)
        return H0 + V_diag

# Test durchführen
print("1. Erstelle Operator...")
operator = SimpleRiemannOperator(L=12.0, n_points=100)

print("2. Starte mit uniformer Dichte...")
rho = np.ones(100) / 11.0  # ∫ρ = 1

print("3. Berechne Hamiltonian...")
H = operator.build_hamiltonian(rho)

print("4. Berechne erste 5 Eigenwerte...")
try:
    eigenvalues, eigenvectors = eigs(H, k=5, which='SR')
    eigenvalues = np.real(eigenvalues)
    print(f"✅ Erfolg! Eigenwerte: {eigenvalues}")
    
    # Einfache Analyse
    print(f"   - Anzahl Eigenwerte: {len(eigenvalues)}")
    print(f"   - Kleinster Eigenwert: {np.min(eigenvalues):.4f}")
    print(f"   - Größter Eigenwert: {np.max(eigenvalues):.4f}")
    
except Exception as e:
    print(f"❌ Fehler: {e}")

print("\n🎉 Test abgeschlossen!")
