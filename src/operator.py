# src/operator.py
import numpy as np
import math

def _primes_upto(N):
    N = int(N)
    if N < 2:
        return []
    sieve = np.ones(N+1, dtype=bool)
    sieve[0:2] = False
    for p in range(2, int(N**0.5)+1):
        if sieve[p]:
            sieve[p*p:N+1:p] = False
    return np.nonzero(sieve)[0].tolist()

class RiemannOperator:
    """
    Builds H[rho] on [1, L] with Dirichlet BC.
    Grid: x[0]=1, x[-1]=L. Interior size M = n_points-2.
    """
    def __init__(self, L=12.0, n_points=1201,
                 eps=None, epsilon=None,          # <-- accept both names
                 gamma=0.35, eta=0.25,
                 sigma=0.5, rho_min=1e-6,
                 P_max=None,
                 **kwargs):                        # <-- ignore extra keys safely
        # handle eps / epsilon alias with a sensible default
        if eps is None and epsilon is not None:
            eps = epsilon
        if eps is None:
            eps = 0.1

        self.L = float(L)
        self.n = int(n_points)
        assert self.n >= 3, "n_points must be >= 3"
        self.x = np.linspace(1.0, self.L, self.n)
        self.h = (self.L - 1.0) / (self.n - 1)
        self.M = self.n - 2  # interior points for Dirichlet

        self.eps = float(eps)
        self.gamma = float(gamma)
        self.eta = float(eta)
        self.sigma = float(sigma)
        self.rho_min = float(rho_min)
        self.P_max = int(P_max) if P_max is not None else int(self.L)

        # precompute kinetic (Dirichlet: interior nodes only)
        self._K = self._build_kinetic_matrix()

    def _build_kinetic_matrix(self):
        """Tri-diagonal matrix for -d^2/dx^2 with Dirichlet BC on interior nodes."""
        M = self.M
        h2 = self.h * self.h
        main = np.full(M, 2.0 / h2, dtype=np.float64)
        off  = np.full(M-1, -1.0 / h2, dtype=np.float64)
        K = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
        return K  # discretization of -d^2/dx^2

    # --------- potentials (unchanged) ----------
    def V_prim(self):
        xp = self.x
        V = np.zeros_like(xp, dtype=np.float64)
        primes = _primes_upto(self.P_max)
        if not primes:
            return V
        inv_two_eps2 = 1.0 / (2.0 * self.eps * self.eps)
        for p in primes:
            lp = math.log(p)
            V += (1.0 / p) * np.exp(-(xp - lp) * (xp - lp) * inv_two_eps2)
        return V

    def V_grav(self, rho):
        xp = self.x
        rho = np.asarray(rho, dtype=np.float64)
        dx = self.h
        V = np.zeros_like(xp, dtype=np.float64)
        for i, xi in enumerate(xp):
            d = xi - xp
            Kxy = (1.0 / np.sqrt(1.0 + d*d)) * np.exp(-0.5 * (d*d) / (self.sigma*self.sigma))
            V[i] = self.gamma * np.sum(Kxy * rho) * dx
        return V

    def V_exch(self, rho):
        rho = np.asarray(rho, dtype=np.float64)
        rho_safe = np.maximum(rho, self.rho_min)
        return -self.eta * np.power(rho_safe, 1.0/3.0)

    def build_hamiltonian(self, rho):
        rho = np.asarray(rho, dtype=np.float64)
        if rho.shape[0] != self.n:
            if rho.shape[0] == self.M:
                rho_full = np.zeros(self.n, dtype=np.float64)
                rho_full[1:-1] = rho
                rho = rho_full
            else:
                raise ValueError(f"rho has shape {rho.shape}, expected ({self.n},) or interior ({self.M},)")
        Vsum = (self.V_prim()[1:-1] + self.V_grav(rho)[1:-1] + self.V_exch(rho)[1:-1]).astype(np.float64)
        H = self._K.copy().astype(np.float64) + np.diag(Vsum)
        H = np.ascontiguousarray(H, dtype=np.float64)
        assert H.ndim == 2 and H.shape[0] == H.shape[1] == self.M, f"H bad shape: {H.shape}"
        return H

# Optional back-compat alias if some code imports OperatorBuilder from src.operator
OperatorBuilder = RiemannOperator
