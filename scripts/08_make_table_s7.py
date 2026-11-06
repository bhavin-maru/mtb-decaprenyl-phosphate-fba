"""Phase 8: build the human-homology table (Supplementary Table S7) for the
genes essential in both the model and strict Tn-seq.

Input : results/phase4_tnseq_detail.csv, results/human_homology_calls.tsv
Output: results/Table_S7_74genes.md
Run   : python scripts/08_make_table_s7.py
"""
import csv

from common import RESULTS_DIR

strict_overlap = set()
with open(RESULTS_DIR / "phase4_tnseq_detail.csv") as fh:
    for row in csv.DictReader(fh):
        if row["model_essential"] == "True" and row["tnseq_strict"] == "True":
            strict_overlap.add(row["gene"])

calls = {}
with open(RESULTS_DIR / "human_homology_calls.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        calls[row["query"]] = row

out_path = RESULTS_DIR / "Table_S7_74genes.md"
with open(out_path, "w") as out:
    out.write("| Gene | Human homolog | Best E-value | Query coverage (%) | Best % identity |\n")
    out.write("|---|---|---|---|---|\n")
    for gene in sorted(strict_overlap):
        res = calls.get(gene, {})
        if res.get("human_homolog") == "YES":
            out.write(f"| {gene} | YES | {res['best_evalue']} | "
                      f"{res['best_qcov']} | {res['best_pident']} |\n")
        else:
            out.write(f"| {gene} | NO | - | - | - |\n")

print(f"Genes in table: {len(strict_overlap)}")
print(f"Saved: {out_path}")
