"""Phase 5: fetch protein sequences of the genes essential in both the model
and strict Tn-seq (the 74-gene set) and write the BLASTp query file.

Sequences come from UniProt (M. tuberculosis H37Rv, taxon 83332). Each FASTA
header is replaced by the bare Rv ID so that script 07 can match hits to genes.

Input : results/phase4_tnseq_detail.csv
Output: blast/queries.faa, results/phase5_query_sequences.csv
Run   : python scripts/05_fetch_query_sequences.py   (needs internet)

Manual step afterwards: BLASTp blast/queries.faa against Homo sapiens
(taxid 9606), RefSeq Select proteins, default parameters, and save the hits
as blast/blast_hits.tsv (format described in script 07).
"""
import time

import pandas as pd
import requests

from common import BLAST_DIR, RESULTS_DIR

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"

detail = pd.read_csv(RESULTS_DIR / "phase4_tnseq_detail.csv")
genes = sorted(detail[(detail["model_essential"]) & (detail["tnseq_strict"])]["gene"])
print(f"Genes essential in the model and strict Tn-seq: {len(genes)}")

records, rows = [], []
for gene in genes:
    params = {"query": f"gene:{gene} AND organism_id:83332",
              "format": "fasta", "size": 1}
    try:
        r = requests.get(UNIPROT_URL, params=params, timeout=30)
        ok = r.status_code == 200 and r.text.startswith(">")
    except requests.RequestException as exc:
        print(f"  {gene}: request failed ({exc})")
        ok = False
    if ok:
        sequence = "".join(r.text.splitlines()[1:])
        records.append(f">{gene}\n{sequence}\n")
        print(f"  {gene}: fetched ({len(sequence)} aa)")
    else:
        print(f"  {gene}: NOT FOUND")
    rows.append({"gene": gene, "fetched": ok})
    time.sleep(0.2)

(BLAST_DIR / "queries.faa").write_text("".join(records))
pd.DataFrame(rows).to_csv(RESULTS_DIR / "phase5_query_sequences.csv", index=False)

n_ok = sum(r["fetched"] for r in rows)
print(f"\nWrote {n_ok}/{len(genes)} sequences to {BLAST_DIR / 'queries.faa'}")
if n_ok < len(genes):
    print("Genes not found must be added to queries.faa by hand before running BLAST.")
