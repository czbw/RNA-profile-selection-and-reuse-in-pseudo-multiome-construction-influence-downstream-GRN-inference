'Shared path, configuration and cell-identifier helpers for the PBMC bundle. AMLGRN_ROOT can override the project root.'
import os
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parent
REPO_ROOT = Path(os.environ.get("AMLGRN_ROOT", str(DEFAULT_ROOT)))

STAGES = ["Diagnosis", "Remission", "Relapse"]
HEALTHY_LABEL = "Normal"
TOKEN2TP = {"DX": "Diagnosis", "REM": "Remission", "REL": "Relapse"}

def load_config():
    # An environment override selects an alternative bundle configuration.
    cfg_name = os.environ.get("AMLGRN_CONFIG", "config.yaml")
    cfg = yaml.safe_load((HERE / cfg_name).read_text(encoding="utf-8"))
    cfg["P"] = {k: (REPO_ROOT / v) for k, v in cfg["paths"].items()}
    cfg["REPO_ROOT"] = REPO_ROOT
    return cfg

def stage_key(stage: str) -> str:
    """'Diagnosis' -> 'diagnosis' (matches the RNA h5ad filenames)."""
    return stage.lower()

def cell_id(barcode: str, sample_id: str, sep: str = "-") -> str:
    """Globally-unique cell name. Must match the cisTopic object's cell names
    after merge (asserted in step 05)."""
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
    print(f"[armA] {msg}", flush=True)
