# Proof 13 – L→∞ Limitcheck (CSV)
- CSV-Muster: `zetaH_L{L}.csv`
- L-Werte: `[8.0, 10.0, 12.0, 14.0]`
- s_ref: `3.0`
- **Fit:** `sup ≈ a*(1/L) + b`  mit  `a=32.2205`, `b=-1.45651`, `R²=0.904007`

## Kennzahlen
|   L |      invL |   sup_dev |   mean_dev |   median_dev |   s_ref_eff |
|----:|----------:|----------:|-----------:|-------------:|------------:|
|   8 | 0.125     |   2.75368 |   1.2324   |      1.00094 |           3 |
|  10 | 0.1       |   1.45734 |   1.02688  |      1.00585 |           3 |
|  12 | 0.0833333 |   1.14642 |   0.963392 |      1.00205 |           3 |
|  14 | 0.0714286 |   1.05265 |   0.94457  |      1.0007  |           3 |

## Bewertung
- **Konvergenz bestätigt**, wenn `b ≈ 0` (klein) und `R²` nahe 1, sowie die Kurven `Ψ_rel(s)` an 1 kleben.