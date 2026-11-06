# mtb-decaprenyl-phosphate-fba

Python pipeline for flux balance analysis (FBA) of decaprenyl phosphate-dependent cell wall biosynthesis in *Mycobacterium tuberculosis* H37Rv, using the genome-scale model iEK1011.

This repository contains the code for the manuscript:

> Maru B, Shah A. *Flux Balance Analysis of Decaprenyl Phosphate-Dependent Cell Wall Biosynthesis in Mycobacterium tuberculosis Identifies Constraints on Drug Synergy and Prioritizes Candidate Drug Targets.* (Manuscript in preparation.)

Archived release: [Zenodo DOI to be added]

## What the pipeline does

1. Adds a cell wall biomass objective (`CW_BIOMASS`) to iEK1011 and identifies the reactions that consume decaprenyl phosphate and related carrier metabolites.
2. Calculates flux control coefficients (FCC) for the active carrier-consuming reactions and tests their stability across four objectives.
3. Deletes all 1,011 genes in silico and classifies them as essential, reduced or non-essential.
4. Compares model-essential genes with Tn-seq essentiality calls (strict and lenient definitions).
5. Simulates pairwise inhibition of the active enzymes and scores synergy with the Bliss independence and highest single agent (HSA) models.
6. Filters candidate targets for human homologues (BLASTp) and prepares the supplementary table.

## Repository structure

```
mtb-decaprenyl-phosphate-fba/
├── README.md
├── LICENSE
├── requirements.txt
├── data/
│   ├── models/iEK1011_m7H10_media.json    # iEK1011 with modified Middlebrook 7H10 medium
│   └── tnseq_essentiality.csv             # Tn-seq calls, DeJesus et al. 2017 (see "Data")
├── scripts/
│   ├── common.py                          # shared settings and helper functions
│   ├── 01_flux_control.py
│   ├── 02_gene_deletion.py
│   ├── 03_objective_robustness.py
│   ├── 04_tnseq_validation.py
│   ├── 05_fetch_query_sequences.py
│   ├── 06_pairwise_synergy.py
│   ├── 07_human_homology_calls.py
│   └── 08_make_table_s7.py
├── blast/
│   ├── queries.faa                        # protein sequences used as BLASTp queries
│   └── blast_hits.tsv                     # BLASTp hits against human proteins
└── results/                               # output tables
```

## Requirements

- Python 3.12
- COBRApy 0.32.1
- pandas 2.3.3
- requests (for script 05)
- GLPK solver (installed with COBRApy)

Install the Python packages with:

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Data

- **Metabolic model.** The model iEK1011 (Kavvas et al. 2018, *BMC Systems Biology* 12:25) is available at <https://github.com/erolkavvas/iEK1011_materials>. The version used here is the one with modified Middlebrook 7H10 medium, saved as `data/models/iEK1011_m7H10_media.json`. The model is distributed under its own licence.
- **Tn-seq essentiality.** The calls come from DeJesus et al. 2017, *mBio* 8(1) (Comprehensive Essentiality Analysis of the *Mycobacterium tuberculosis* Genome via Saturating Transposon Mutagenesis). Download the supplementary table from the article and save it as `data/tnseq_essentiality.csv`. The script expects a descriptive text line in the first row and the column header (including `ORF ID` and `Final Call`) in the second row.

## Usage

Run every script from the repository root, because all paths are relative.

```bash
python scripts/01_flux_control.py
python scripts/02_gene_deletion.py          # about 5-15 minutes
python scripts/03_objective_robustness.py
python scripts/04_tnseq_validation.py
python scripts/05_fetch_query_sequences.py  # needs internet access
# Manual step: run BLASTp (see below), then:
python scripts/07_human_homology_calls.py
python scripts/08_make_table_s7.py

python scripts/06_pairwise_synergy.py       # independent; can be run at any time
```

| Script | Input | Output |
|---|---|---|
| `01_flux_control.py` | model | `results/phase1_fcc.csv` |
| `02_gene_deletion.py` | model | `results/phase2_deletions.csv` |
| `03_objective_robustness.py` | model | `results/phase3_objective_robustness.csv` |
| `04_tnseq_validation.py` | `phase2_deletions.csv`, Tn-seq table | `results/phase4_tnseq_comparison.csv`, `results/phase4_tnseq_detail.csv` |
| `05_fetch_query_sequences.py` | `phase4_tnseq_detail.csv` | `blast/queries.faa`, `results/phase5_query_sequences.csv` |
| `06_pairwise_synergy.py` | model | `results/pairwise_synergy.csv` |
| `07_human_homology_calls.py` | `queries.faa`, `blast_hits.tsv` | `results/human_homology_calls.tsv` |
| `08_make_table_s7.py` | `phase4_tnseq_detail.csv`, `human_homology_calls.tsv` | `results/Table_S7_74genes.md` |

### BLASTp step (manual)

Script 05 writes `blast/queries.faa`, with one sequence per gene and the Rv identifier as the header. Search these sequences with NCBI BLASTp against *Homo sapiens* (taxon 9606), RefSeq Select proteins, with default algorithm parameters. Save the hits as `blast/blast_hits.tsv`, a tab-separated file without a header and with six columns:

```
query   subject   evalue   query_coverage   percent_identity   alignment_length
```

A gene is classified as having a human homologue when at least one hit has an E-value below 1 × 10⁻⁵ and a query coverage of at least 70%.

## Main analysis settings

| Parameter | Value |
|---|---|
| Model | iEK1011, modified m7H10 medium (1,011 genes) |
| Objective | `CW_BIOMASS` (arabinogalactan, mycolic acid, mycolylated arabinan, peptidoglycan) |
| Flux control perturbation | upper bound lowered by 1% (δ = 0.01) |
| Activity threshold | flux above 10⁻⁵ mmol gDW⁻¹ h⁻¹ |
| Essentiality threshold | at least 90% drop in objective flux |
| Inhibition levels | 30%, 50%, 70% (45 pairs, 405 pair-dose combinations) |
| Synergy threshold | at least 10% excess over the expected effect |
| Tn-seq definitions | strict (ES); lenient (ES + ESD + GD) |
| Solver | GLPK |

## Summary of results reproduced by the pipeline

- 19 carrier-consuming reactions were identified; 10 were active, encoded by 14 genes.
- 100 genes were model-essential; 74 were essential in strict Tn-seq (precision 74.0%) and 84 in the lenient set (precision 84.0%).
- No pair was synergistic under Bliss independence (0 of 405; maximum +4.3%), whereas the HSA model called 147 of 405.
- 45 of the 74 genes had no qualifying human homologue.

## Reproducibility notes

- Results depend on the model file and on the versions of COBRApy and GLPK. Use the versions listed above.
- Script 05 retrieves sequences from UniProt at run time. The sequences used in the paper are stored in `blast/queries.faa`.
- Tn-seq data and the iEK1011 model are third-party resources and are subject to their own licences.

## Citation

If you use this code, please cite the manuscript above and the archived release:

```
Maru B, Shah A. mtb-decaprenyl-phosphate-fba (version 1.0.0). Zenodo. [DOI to be added]
```

## Licence

The code in this repository is released under the MIT Licence (see `LICENSE`).

## Contact

Bhavin Maru (GitHub: [@bhavin-maru](https://github.com/bhavin-maru))
Department of Pharmacy, Sumandeep Vidyapeeth Deemed to be University, Vadodara, Gujarat, India
