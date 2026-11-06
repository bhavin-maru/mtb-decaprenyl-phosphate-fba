"""Phase 1: flux control coefficients of the carrier-consuming reactions.

Input : data/models/iEK1011_m7H10_media.json
Output: results/phase1_fcc.csv
Run   : python scripts/01_flux_control.py   (from the repository root)
"""
import warnings

import pandas as pd

from common import OBJECTIVES, RESULTS_DIR, classify_fcc, compute_fcc, setup_model

warnings.filterwarnings("ignore")

model, carrier_consumers, obj_id = setup_model(OBJECTIVES["CW_BIOMASS"], "CW_BIOMASS")

with model as m:
    m.objective = "biomass"
    wt_growth = m.slim_optimize()
    m.objective = obj_id
    wt_cw = m.slim_optimize()
print(f"Wild-type growth: {wt_growth:.4f}; wild-type CW_BIOMASS: {wt_cw:.4f}")
print(f"Carrier-consuming reactions: {len(carrier_consumers)}")

rows = []
for rid in carrier_consumers:
    fcc, flux = compute_fcc(model, rid, obj_id)
    rxn = model.reactions.get_by_id(rid)
    rows.append({
        "reaction": rid,
        "name": rxn.name[:45],
        "gene": rxn.gene_reaction_rule or "(gap-fill)",
        "baseline_flux": round(flux, 6),
        "FCC": round(fcc, 4) if fcc is not None else None,
        "class": classify_fcc(fcc),
    })

df = pd.DataFrame(rows)
df["abs_fcc"] = df["FCC"].abs()
df = df.sort_values("abs_fcc", ascending=False, na_position="last")
df.to_csv(RESULTS_DIR / "phase1_fcc.csv", index=False)

print("\nFCC classes:")
print(df["class"].value_counts().to_string())
print(f"\nSaved: {RESULTS_DIR / 'phase1_fcc.csv'}")
