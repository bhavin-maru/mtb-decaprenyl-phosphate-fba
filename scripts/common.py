"""Shared configuration and helpers for the iEK1011 cell wall FBA pipeline.

Run every script from the repository root, for example:
    python scripts/01_flux_control.py
"""
import re
from pathlib import Path

import cobra
from cobra import Reaction

# ---------------------------------------------------------------- paths
MODEL_PATH = Path("data/models/iEK1011_m7H10_media.json")
TNSEQ_PATH = Path("data/tnseq_essentiality.csv")   # DeJesus et al. 2017, mBio
RESULTS_DIR = Path("results")
BLAST_DIR = Path("blast")
RESULTS_DIR.mkdir(exist_ok=True)
BLAST_DIR.mkdir(exist_ok=True)

# ----------------------------------------------------------- parameters
FCC_DELTA = 0.01              # 1% reduction of the upper bound
MIN_FLUX = 1e-5               # activity threshold (mmol gDW-1 h-1)
ESSENTIAL_THRESHOLD = 90.0    # % drop in objective flux
REDUCED_THRESHOLD = 10.0

# Decaprenyl phosphate and related carrier metabolites
CARRIER_METS = {
    "decd__tb_c", "decda__tb_c", "decdp__tb_c",
    "decdr__tb_c", "decdpr__tb_c", "decdman1p__tb_c",
    "udcpp_c", "udcpdp_c",
}

# CW_BIOMASS stoichiometry (Supplementary Table S1)
CW_COEFFS = {
    "arabinanagalfragund_c": -0.000731,
    "c78mycolatepp_c": -0.0295,
    "peptido__TB1_c": -0.01885,
    "peptido__TB2_c": -0.01885,
    "tmha1_c": -0.001106,
    "tmha2_c": -0.0011,
    "tmha3_c": -0.001104,
    "tmha4_c": -0.001097,
}

# The key is also used as the objective reaction ID
OBJECTIVES = {
    "CW_BIOMASS": CW_COEFFS,
    "MYCOLATE_ONLY": {"c78mycolatepp_c": -1.0},
    "ARABINOGALACTAN_ONLY": {"arabinanagalfragund_c": -1.0},
    "PEPTIDOGLYCAN_ONLY": {"peptido__TB1_c": -1.0, "peptido__TB2_c": -1.0},
}


# ------------------------------------------------------------- helpers
def load_model():
    """Load iEK1011 (m7H10 medium) and use the GLPK solver."""
    model = cobra.io.load_json_model(MODEL_PATH)
    model.solver = "glpk"
    return model


def add_objective(model, name, coeffs):
    """Add a cell wall objective reaction and make it the model objective."""
    rxn = Reaction(name)
    rxn.name = f"Cell wall objective ({name})"
    rxn.lower_bound, rxn.upper_bound = 0, 1000
    rxn.add_metabolites({model.metabolites.get_by_id(k): v
                         for k, v in coeffs.items()})
    model.add_reactions([rxn])
    model.objective = name
    return rxn


def setup_model(coeffs, name):
    """Load the model, add the objective, force carrier consumers forward.

    Returns (model, list of carrier-consuming reaction IDs, objective ID).
    """
    model = load_model()
    add_objective(model, name, coeffs)

    consumers = []
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        for met in rxn.metabolites:
            if met.id in CARRIER_METS and rxn.metabolites[met] < 0:
                consumers.append(rxn.id)
                break

    for rid in consumers:
        rxn = model.reactions.get_by_id(rid)
        if rxn.lower_bound < 0:
            rxn.lower_bound = 0
    return model, consumers, name


def compute_fcc(model, rxn_id, objective_id, delta=FCC_DELTA):
    """Flux control coefficient (Equation 2). Returns (FCC or None, baseline flux)."""
    with model as m:
        m.objective = objective_id
        sol = m.optimize()
        base = sol.fluxes[objective_id]
        flux = sol.fluxes.get(rxn_id, 0.0)
    if base < 1e-9 or flux < MIN_FLUX:
        return None, flux
    with model as m:
        m.reactions.get_by_id(rxn_id).upper_bound = flux * (1 - delta)
        m.objective = objective_id
        down = m.optimize().fluxes[objective_id]
    return (base - down) / base / delta, flux


def classify_fcc(f):
    if f is None:
        return "excluded"
    if abs(f) > 0.5:
        return "dominant"
    if abs(f) > 0.1:
        return "strong"
    if abs(f) > 0.01:
        return "moderate"
    return "negligible"


def extract_rv(series):
    """Return the set of Rv IDs found in a pandas Series."""
    return set(series.astype(str).str.strip()
               .str.extract(r"(Rv\d+[A-Za-z]?)", expand=False).dropna())
