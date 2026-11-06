"""Phase 6: pairwise inhibition of the active carrier-consuming enzymes.

Inhibition is modelled by diverting a fraction of each product of the target
reaction to a waste metabolite. Combination effects are scored against the
Bliss independence (Equation 4) and highest single agent (Equation 5) models.

Input : data/models/iEK1011_m7H10_media.json
Output: results/pairwise_synergy.csv
Run   : python scripts/06_pairwise_synergy.py
"""
import itertools
import warnings

import pandas as pd
from cobra import Metabolite, Reaction

from common import CW_COEFFS, RESULTS_DIR, add_objective, load_model

warnings.filterwarnings("ignore")

ACTIVE_ENZYMES = ["DCPE", "DCPT", "DCPDPP", "UDCPDP", "PPM1",
                  "PAPPT1", "PAPPT2", "EMB", "AFTA", "ACGAMT"]
INHIB_FRACS = [0.3, 0.5, 0.7]
SYN_THRESHOLD = 10.0     # percentage points

model = load_model()
add_objective(model, "CW_BIOMASS", CW_COEFFS)

ACTIVE_ENZYMES = [r for r in ACTIVE_ENZYMES if r in model.reactions]
print(f"Enzymes tested: {len(ACTIVE_ENZYMES)}")

# One waste sink per enzyme
for rid in ACTIVE_ENZYMES:
    rxn = model.reactions.get_by_id(rid)
    if rxn.lower_bound < 0:
        rxn.lower_bound = 0
    waste = Metabolite(f"waste_{rid}")
    sink = Reaction(f"DM_waste_{rid}")
    sink.add_metabolites({waste: -1})
    sink.lower_bound, sink.upper_bound = 0, 1000
    model.add_reactions([sink])

with model as m:
    m.objective = "CW_BIOMASS"
    base_cw = m.optimize().fluxes["CW_BIOMASS"]
print(f"Baseline CW_BIOMASS flux: {base_cw:.6f}")


def perturb(m, rid, frac):
    """Divert `frac` of every product of reaction `rid` to its waste metabolite."""
    rxn = m.reactions.get_by_id(rid)
    waste = m.metabolites.get_by_id(f"waste_{rid}")
    changes, total = {}, 0.0
    for met, coeff in rxn.metabolites.items():
        if coeff > 0:
            changes[met] = -coeff * frac
            total += coeff * frac
    if total > 0:
        changes[waste] = total
        rxn.add_metabolites(changes, combine=True)


def run_fba(perts):
    """perts = [(reaction ID, inhibition fraction), ...]; returns flux or None."""
    with model as m:
        for rid, frac in perts:
            if frac > 0:
                perturb(m, rid, frac)
        m.objective = "CW_BIOMASS"
        sol = m.optimize()
        return sol.fluxes["CW_BIOMASS"] if sol.status == "optimal" else None


# Single-inhibitor effects
singles = {}
for rid in ACTIVE_ENZYMES:
    for d in INHIB_FRACS:
        flux = run_fba([(rid, d)])
        singles[(rid, d)] = None if flux is None else max(0.0, (base_cw - flux) / base_cw)
print(f"Single-inhibitor runs feasible: "
      f"{sum(v is not None for v in singles.values())}/{len(singles)}")

# Pairwise combinations
rows = []
for r1, r2 in itertools.combinations(ACTIVE_ENZYMES, 2):
    for d1 in INHIB_FRACS:
        for d2 in INHIB_FRACS:
            s1, s2 = singles[(r1, d1)], singles[(r2, d2)]
            if s1 is None or s2 is None:
                continue
            flux = run_fba([(r1, d1), (r2, d2)])
            if flux is None:
                rows.append({"Drug_1": r1, "Drug_2": r2, "Dose_1": d1, "Dose_2": d2,
                             "D1_pct": round(s1 * 100, 2), "D2_pct": round(s2 * 100, 2),
                             "Combo_pct": None, "Bliss_expected": None,
                             "HSA_expected": round(max(s1, s2) * 100, 2),
                             "Bliss_synergy": None, "HSA_synergy": None,
                             "Status": "Infeasible"})
                continue
            combo = max(0.0, (base_cw - flux) / base_cw)
            bliss = s1 + s2 - s1 * s2
            hsa = max(s1, s2)
            rows.append({"Drug_1": r1, "Drug_2": r2, "Dose_1": d1, "Dose_2": d2,
                         "D1_pct": round(s1 * 100, 2), "D2_pct": round(s2 * 100, 2),
                         "Combo_pct": round(combo * 100, 2),
                         "Bliss_expected": round(bliss * 100, 2),
                         "HSA_expected": round(hsa * 100, 2),
                         "Bliss_synergy": round((combo - bliss) * 100, 2),
                         "HSA_synergy": round((combo - hsa) * 100, 2),
                         "Status": "Optimal"})

df = pd.DataFrame(rows)
df.to_csv(RESULTS_DIR / "pairwise_synergy.csv", index=False)

ok = df[df["Status"] == "Optimal"]
print(f"\nPair-dose combinations: {len(df)} (feasible: {len(ok)})")
print(f"Bliss range: {ok['Bliss_synergy'].min():.1f}% to {ok['Bliss_synergy'].max():.1f}%")
print(f"Bliss >= {SYN_THRESHOLD}%: {(ok['Bliss_synergy'] >= SYN_THRESHOLD).sum()}")
print(f"HSA   >= {SYN_THRESHOLD}%: {(ok['HSA_synergy'] >= SYN_THRESHOLD).sum()}")
print("\nTop 10 by Bliss synergy:")
print(ok.nlargest(10, "Bliss_synergy")[
    ["Drug_1", "Drug_2", "Dose_1", "Dose_2", "Combo_pct", "Bliss_synergy"]
].to_string(index=False))
print(f"\nSaved: {RESULTS_DIR / 'pairwise_synergy.csv'}")
