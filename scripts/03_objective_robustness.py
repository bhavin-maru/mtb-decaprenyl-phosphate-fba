"""Phase 3: stability of FCC across four cell wall objectives (Equation 7).

Input : data/models/iEK1011_m7H10_media.json
Output: results/phase3_objective_robustness.csv
Run   : python scripts/03_objective_robustness.py
"""
import warnings

import pandas as pd

from common import OBJECTIVES, RESULTS_DIR, compute_fcc, setup_model

warnings.filterwarnings("ignore")

SD_THRESHOLD = 0.05

rows = []
for name, coeffs in OBJECTIVES.items():
    print(f"\nObjective: {name}")
    model, consumers, oid = setup_model(coeffs, name)
    with model as m:
        m.objective = oid
        base = m.slim_optimize()
    print(f"  Baseline flux: {base:.6f}")
    if base < 1e-6:
        print("  Zero baseline; objective skipped")
        continue
    for rid in consumers:
        fcc, _ = compute_fcc(model, rid, oid)
        if fcc is not None:
            rows.append({"objective": name, "reaction": rid, "FCC": round(fcc, 4)})

df = pd.DataFrame(rows)
df.to_csv(RESULTS_DIR / "phase3_objective_robustness.csv", index=False)

if not df.empty:
    pivot = df.pivot_table(index="reaction", columns="objective",
                           values="FCC", aggfunc="first").round(3)
    print("\nFCC across objectives:")
    print(pivot.to_string())

    tested = stable = 0
    for rid in pivot.index:
        vals = pivot.loc[rid].dropna()
        if len(vals) >= 2:                 # active under at least two objectives
            tested += 1
            if vals.std() < SD_THRESHOLD:  # sample SD, n - 1
                stable += 1
    print(f"\nStable reactions (SD < {SD_THRESHOLD}): {stable}/{tested}")

print(f"\nSaved: {RESULTS_DIR / 'phase3_objective_robustness.csv'}")
