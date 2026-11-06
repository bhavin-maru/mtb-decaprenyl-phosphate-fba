"""Phase 7: apply the human-homology filter to the BLASTp hits.

A gene is called as having a human homologue when at least one hit has an
E-value below 1e-5 and query coverage of at least 70%.

Input : blast/queries.faa            (from script 05)
        blast/blast_hits.tsv         (tab-separated, no header, six columns:
                                      query, subject, evalue, qcov, pident, length)
Output: results/human_homology_calls.tsv
Run   : python scripts/07_human_homology_calls.py
"""
import csv
from pathlib import Path

from common import BLAST_DIR, RESULTS_DIR

EVALUE_MAX = 1e-5
QCOV_MIN = 70.0

queries = [line[1:].split()[0]
           for line in (BLAST_DIR / "queries.faa").read_text().splitlines()
           if line.startswith(">")]

best = {}
with open(BLAST_DIR / "blast_hits.tsv") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if len(row) < 6:
            continue
        q, subject, evalue, qcov, pident, length = row[:6]
        if float(evalue) < EVALUE_MAX and float(qcov) >= QCOV_MIN:
            if q not in best or float(evalue) < best[q][1]:
                best[q] = (subject, float(evalue), float(qcov), float(pident))

out_path = RESULTS_DIR / "human_homology_calls.tsv"
with open(out_path, "w") as out:
    out.write("query\thuman_homolog\tbest_subject\tbest_evalue\tbest_qcov\tbest_pident\n")
    for q in queries:
        if q in best:
            s, e, qc, pi = best[q]
            out.write(f"{q}\tYES\t{s}\t{e:g}\t{qc:g}\t{pi:g}\n")
        else:
            out.write(f"{q}\tNO\t\t\t\t\n")

n_yes = sum(q in best for q in queries)
print(f"{len(queries)} queries: {n_yes} with a human homologue, "
      f"{len(queries) - n_yes} without")
print(f"Saved: {out_path}")
