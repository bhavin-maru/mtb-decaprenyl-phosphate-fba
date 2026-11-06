"""Phase 4: compare model-essential genes with Tn-seq essentiality (Equation 6).

Tn-seq data: DeJesus et al. 2017, mBio (supplementary file). Save the table
as data/tnseq_essentiality.csv. The first row of that file is a descriptive
text blob, and the real header is on the second row (skiprows=1).

Input : results/phase2_deletions.csv, data/tnseq_essentiality.csv
Output: results/phase4_tnseq_comparison.csv, results/phase4_tnseq_detail.csv
Run   : python scripts/04_tnseq_validation.py
"""
import sys

import pandas as pd

from common import RESULTS_DIR, TNSEQ_PATH, extract_rv

STRICT_CALLS = {"ES"}
LENIENT_CALLS = {"ES", "ESD", "GD"}

tn = pd.read_csv(TNSEQ_PATH, skiprows=1)
print(f"Loaded Tn-seq table: {len(tn)} rows")

gene_col = "ORF ID" if "ORF ID" in tn.columns else next(
    (c for c in tn.columns if "ORF" in str(c)), tn.columns[0])
call_col = "Final Call" if "Final Call" in tn.columns else next(
    (c for c in tn.columns if "call" in str(c).lower()), tn.columns[-1])
print(f"Gene column: {gene_col}; call column: {call_col}")
if call_col not in tn.columns or gene_col not in tn.columns:
    sys.exit("Header not recognised. Check skiprows and the file format.")

calls = tn[call_col].astype(str).str.strip()
print("\nFinal Call distribution:")
print(calls.value_counts().to_string())

tn_strict = extract_rv(tn[calls.isin(STRICT_CALLS)][gene_col])
tn_lenient = extract_rv(tn[calls.isin(LENIENT_CALLS)][gene_col])
print(f"\nTn-seq strict (ES): {len(tn_strict)} genes")
print(f"Tn-seq lenient (ES+ESD+GD): {len(tn_lenient)} genes")

deletions = pd.read_csv(RESULTS_DIR / "phase2_deletions.csv")
model_ess = extract_rv(deletions[deletions["status"] == "essential"]["gene"])
print(f"Model essential: {len(model_ess)} genes")


def compare(model_set, tn_set, label):
    overlap = model_set & tn_set
    precision = len(overlap) / len(model_set) if model_set else 0
    recall = len(overlap) / len(tn_set) if tn_set else 0
    print(f"\n{label}: TP={len(overlap)}  FP={len(model_set - tn_set)}  "
          f"FN={len(tn_set - model_set)}  precision={precision:.1%}  recall={recall:.1%}")
    return {
        "definition": label,
        "model_n": len(model_set),
        "tnseq_n": len(tn_set),
        "overlap_TP": len(overlap),
        "model_only_FP": len(model_set - tn_set),
        "tnseq_only_FN": len(tn_set - model_set),
        "precision": round(precision, 4),
        "recall": round(recall, 4),
    }


summary = pd.DataFrame([
    compare(model_ess, tn_strict, "Strict (ES)"),
    compare(model_ess, tn_lenient, "Lenient (ES+ESD+GD)"),
])
summary.to_csv(RESULTS_DIR / "phase4_tnseq_comparison.csv", index=False)

detail = pd.DataFrame([
    {"gene": g,
     "model_essential": g in model_ess,
     "tnseq_strict": g in tn_strict,
     "tnseq_lenient": g in tn_lenient}
    for g in sorted(model_ess | tn_lenient)
])
detail.to_csv(RESULTS_DIR / "phase4_tnseq_detail.csv", index=False)

print(f"\nSaved: {RESULTS_DIR / 'phase4_tnseq_comparison.csv'}")
print(f"Saved: {RESULTS_DIR / 'phase4_tnseq_detail.csv'}")
