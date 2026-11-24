"""Extract the per-cell-type metacell counts SCENIC+ reported for the Stage 2 AML pseudomultiome arm.

SCENIC+ logs "Automatically set `nr_metacells` to: <celltype>: <n>, ..." when it builds the
pseudo-multiome data. The count of metacells built is the sum of those per-cell-type values. Counts are read from the original SCENIC+ job logs supplied under data/logs/.

Writes analysis/out/aml_n_metacells.json.
"""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOGS = {"Diagnosis": ROOT / "data/logs/diagnosis_metacells.log",
        "Relapse": ROOT / "data/logs/relapse_metacells.log"}
PATTERN = re.compile(r"`nr_metacells` to: ([^\n]+)")
OUT = ROOT / "analysis/out/aml_n_metacells.json"

def main():
    result = {"k_cells_per_metacell_note": "per-cell-type counts are exactly as SCENIC+ logged them",
              "stages": {}}
    for stage, log in LOGS.items():
        text = log.read_text(encoding="utf-8", errors="replace")
        matches = PATTERN.findall(text)
        if len(matches) != 1:
            raise ValueError("expected one nr_metacells line in %s, found %d" % (log, len(matches)))
        pairs = re.findall(r"([A-Za-z0-9._\-]+): (\d+)", matches[0])
        counts = {name: int(value) for name, value in pairs}
        result["stages"][stage] = {
            "source_log": str(log.relative_to(ROOT)),
            "source_log_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            "n_cell_types": len(counts),
            "per_cell_type": counts,
            "n_metacells_total": sum(counts.values())}
    result["caveat"] = ("Diagnosis and Relapse counts come from the first run that built each stage's "
                        "pseudomultiome data (jobs 16435359 and 16882937). Later Diagnosis jobs "
                        "(16882978, 16949868) contain no ingestion step and reused it; that reuse is "
                        "inferred from the absence of the step in their logs, not from a recorded statement.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(json.dumps({k: v["n_metacells_total"] for k, v in result["stages"].items()}))

if __name__ == "__main__":
    main()
