# src/paper_results_summary.py
import json
import os
import csv
import time

def run():
    """Collects all numerical results (raw + rescaled) into one summary file."""
    print("🧮 Building paper summary bundle...")
    base = "data"
    raw_path = os.path.join(base, "results.json")
    res_path = os.path.join(base, "results_rescaled.json")

    if not os.path.exists(raw_path) or not os.path.exists(res_path):
        print("❌ Missing results. Run main pipeline and rescaling first.")
        return

    with open(raw_path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    with open(res_path, "r", encoding="utf-8") as f:
        res = json.load(f)

    # Extract some top-level values
    zeta_raw = raw.get("zeta", {})
    zeta_rescaled = res.get("zeta_rescaled", {})

    first_raw_key = next(iter(zeta_raw), None)
    trunc_K_raw = zeta_raw[first_raw_key].get("truncation_K", None) if first_raw_key else None

    summary = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "domain": raw.get("precheck", {}).get("domain", {}),
        "bounds": raw.get("precheck", {}).get("bounds", {}),
        "spectral": raw.get("spectral", {}),
        "K_raw": trunc_K_raw,
        "K_rescaled": res.get("K", None),
        "lambda_used": res.get("lambda_used", None),
        "lambda_theoretical": res.get("lambda_theoretical", None),
        "zeta_raw": zeta_raw,
        "zeta_rescaled": zeta_rescaled,
    }

    os.makedirs(base, exist_ok=True)
    json_path = os.path.join(base, "paper_summary.json")
    csv_path = os.path.join(base, "paper_summary.csv")

    # Write JSON
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"💾 JSON summary written to {json_path}")

    # Write quick CSV summary for LaTeX
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        dom = summary.get("domain", {})
        bnd = summary.get("bounds", {})
        spec = summary.get("spectral", {})
        writer.writerow(["L", dom.get("L")])
        writer.writerow(["n_points", dom.get("n_points")])
        writer.writerow(["||V_prim||∞", bnd.get("Vprim_max")])
        writer.writerow(["||V_grav||∞", bnd.get("Vgrav_max")])
        writer.writerow(["||V_exch||∞", bnd.get("Vexch_max")])
        writer.writerow(["K_raw", summary.get("K_raw")])
        writer.writerow(["K_rescaled", summary.get("K_rescaled")])
        writer.writerow(["lambda_used", summary.get("lambda_used")])
        writer.writerow(["lambda_theoretical", summary.get("lambda_theoretical")])
        if spec:
            weyl = spec.get("weyl", {})
            gaps = spec.get("gaps", {})
            writer.writerow(["Emax", weyl.get("Emax")])
            writer.writerow(["count", weyl.get("count")])
            writer.writerow(["weyl_leading", weyl.get("weyl_leading")])
            writer.writerow(["KS_GOE", gaps.get("KS_GOE")])
            writer.writerow(["KS_Poisson", gaps.get("KS_Poisson")])
    print(f"💾 CSV summary written to {csv_path}")
    print("✅ Paper summary bundle complete.")

if __name__ == "__main__":
    run()
