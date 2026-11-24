"""Shared helpers for STAGE 2 (malignant compartment) SCENIC+ scripts."""
import os
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent                       # .../paml/malignant/scripts
REPO_ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))   # bundle root; run_all.sh exports AMLGRN_ROOT

STAGES = ["Diagnosis", "Relapse"]
TOKEN2TP = {"DX": "Diagnosis", "REM": "Remission", "REL": "Relapse"}

def load_config():
    cfg_name = os.environ.get("AMLGRN_CONFIG", "config.yaml")
    cfg = yaml.safe_load((HERE / cfg_name).read_text(encoding="utf-8"))
    cfg["P"] = {k: (REPO_ROOT / v) for k, v in cfg["paths"].items()}
    cfg["REPO_ROOT"] = REPO_ROOT
    return cfg

def compartment_label(cfg) -> str:
    "The value of the ATAC/RNA `Malignant` field to KEEP. Stage 2 = 'Malignant'\n    (Stage 1 kept 'Normal'). Read from config so it is a single source of native reference."
    return str(cfg.get("compartment_label", "Malignant"))

def rna_prefix(cfg) -> str:
    """RNA h5ad filename prefix: '{prefix}_{stage}.h5ad'. Stage 2 = 'malignant'."""
    return str(cfg.get("rna_prefix", "malignant"))

def config_stages(cfg):
    """Stages actually run, from config (defaults to the module STAGES)."""
    return list(cfg.get("stages", STAGES))

def stage_key(stage: str) -> str:
    """'Diagnosis' -> 'diagnosis' (matches the RNA h5ad filenames)."""
    return stage.lower()

def cell_id(barcode: str, sample_id: str, sep: str = "-") -> str:
    """Globally-unique cell name. Must match the cisTopic object's cell names
    after merge (asserted in the ATAC chain)."""
    return f"{barcode}{sep}{sample_id}"

def parse_sample_stem(name: str):
    """'GSM7498708_AML16_DX' -> dict(gsm, patient_token, tp_token, timepoint)."""
    parts = name.split("_")
    tok = parts[-1].upper()
    return {
        "gsm": parts[0],
        "patient_token": "_".join(parts[1:-1]),
        "tp_token": tok,
        "timepoint": TOKEN2TP.get(tok, tok),
    }

def log(msg: str):
    print(f"[stage2] {msg}", flush=True)
