"""Phase 2: single-gene deletion of all genes under the CW_BIOMASS objective.

Input : data/models/iEK1011_m7H10_media.json
Output: results/phase2_deletions.csv
Run   : python scripts/02_gene_deletion.py   (about 5-15 minutes)
"""
import warnings

import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from common import (ESSENTIAL_THRESHOLD, OBJECTIVES, REDUCED_THRESHOLD,
                    RESULTS_DIR, setup_model)

warnings.filterwarnings("ignore")

model, _, obj_id = setup_model(OBJECTIVES["CW_BIOMASS"], "CW_BIOMASS")
with model as m:
    m.objective = obj_id
    baseline = m.slim_optimize()
print(f"Baseline CW_BIOMASS flux: {baseline:.6f}")
print(f"Deleting {len(model.genes)} genes...")

res = single_gene_deletion(model, model.genes, processes=1)
res["growth"] = res["growth"].fillna(0.0)
res["drop_pct"] = ((baseline - res["growth"]) / baseline * 100).round(2)


def classify(drop):
    if drop >= ESSENTIAL_THRESHOLD:
        return "essential"
    if drop >= REDUCED_THRESHOLD:
        return "reduced"
    return "non-essential"


res["status"] = res["drop_pct"].apply(classify)
res["gene"] = res["ids"].apply(lambda x: ",".join(list(x)) if x else "")
res = res[["gene", "growth", "drop_pct", "status"]].sort_values("drop_pct", ascending=False)
res.to_csv(RESULTS_DIR / "phase2_deletions.csv", index=False)

print("\nDeletion classes:")
print(res["status"].value_counts().to_string())
print(f"\nSaved: {RESULTS_DIR / 'phase2_deletions.csv'}")
