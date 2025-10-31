# src/scf_solver.py
import numpy as np

class SCFSolver:
    """
    Self-Consistent Field (SCF) solver for the Riemann Operator Model.
    Iterates rho_{n+1} = F(rho_n) until convergence.
    """
    def __init__(self, operator, max_iter=100, tol=1e-5, mixing=0.5,
                 N_states=32, beta=None, **kwargs):
        """
        Parameters
        ----------
        operator : RiemannOperator
        max_iter : int
        tol      : float
        mixing   : float
        N_states : int     number of lowest states used to build rho
        beta     : float|None  optional inverse temperature (not used yet)
        **kwargs : ignore any extra config keys to stay forward-compatible
        """
        self.operator = operator
        self.max_iter = max_iter
        self.tol = tol
        self.mixing = mixing
        self.N_states = int(N_states)
        self.beta = beta  # currently unused; reserved for Fermi weights
        self.history = []

    def solve(self, rho_init=None):
        print("🔄 Starting SCF...")
        if rho_init is None:
            x = self.operator.x
            rho = np.ones_like(x)
            rho /= np.trapz(rho, x)
        else:
            rho = rho_init.copy()

        for k in range(1, self.max_iter + 1):
            H = self.operator.build_hamiltonian(rho)
            evals, evecs = np.linalg.eigh(H)

            # evecs are on interior nodes (length M = n-2)
            N = min(self.N_states, len(evals))
            psi = evecs[:, :N]                             # shape (M, N)

            # density on interior nodes
            rho_new_int = np.sum(np.abs(psi)**2, axis=1)  # shape (M,)

            # pad to full grid (Dirichlet zeros at boundaries) for consistent normalization
            rho_new = np.zeros_like(self.operator.x, dtype=float)  # length n
            rho_new[1:-1] = rho_new_int

            # normalize on full grid
            area = np.trapz(rho_new, self.operator.x)
            if area <= 0:
                raise ValueError("Non-positive density area during SCF normalization.")
            rho_new /= area

            # keep the mixing as you had it:
            rho_mixed = (1 - self.mixing) * rho + self.mixing * rho_new
            d_norm = np.linalg.norm(rho_mixed - rho)

            self.history.append({
                "iter": k,
                "norm": float(d_norm),
                "mix": float(self.mixing),
                "tol": float(self.tol),
                "N_states": int(N),
                "beta": None if self.beta is None else float(self.beta),
            })

            if k % 5 == 1 or d_norm < self.tol:
                print(f"   Iter {k:3d}: ||Δρ||₂ = {d_norm:8.3e}")

            if d_norm < self.tol:
                print(f"✅ SCF converged in {k} iterations.")
                break

            rho = rho_mixed
        else:
            print("⚠️ SCF did not converge within max_iter limit.")

        return rho, evals, evecs, self.history
