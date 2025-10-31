# src/spectral_analyzer.py
import numpy as np
from scipy import stats

class SpectralAnalyzer:
    def __init__(self, x, eigenvalues):
        self.x = np.asarray(x, dtype=float)
        self.E = np.sort(np.asarray(eigenvalues, dtype=float))

    def weyl_sanity(self, Emax=None):
        if Emax is None:
            Emax = np.percentile(self.E, 90.0)
        count = np.sum(self.E <= Emax)
        L = self.x[-1] - self.x[0]
        weyl = (L / np.pi) * np.sqrt(max(Emax, 0))
        return {"Emax": float(Emax), "count": int(count), "weyl_leading": float(weyl)}

    def spacing_statistics(self):
        d = np.diff(self.E)
        d = d[d > 0]
        if d.size < 5:
            return {"n": int(d.size)}
        d_mean = np.mean(d)
        s = d / d_mean

        # Wigner surmise (GOE)
        def cdf_goe(x):
            return 1.0 - np.exp(-np.pi * x**2 / 4.0)
        # Poisson
        def cdf_pois(x):
            return 1.0 - np.exp(-x)

        stat_goe, _ = stats.kstest(s, cdf_goe)
        stat_poi, _ = stats.kstest(s, cdf_pois)

        return {"n": int(s.size), "KS_GOE": float(stat_goe), "KS_Poisson": float(stat_poi)}
