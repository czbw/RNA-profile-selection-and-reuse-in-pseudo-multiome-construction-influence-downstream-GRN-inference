"""Add the arm_keyed/vocabularies fields to a legacy module definition so it can be scored under one arm.

Verifies that exactly those two keys were added and every other key is unchanged."""
import argparse
import json
from pathlib import Path

ADDED = ("arm_keyed", "vocabularies")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--arm", required=True, choices=["link", "repair"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    src = Path(a.src)
    raw = json.loads(src.read_text())
    if any(k in raw for k in ADDED):
        raise SystemExit(f"[shim] {src} already carries {ADDED}; it does not need this shim")

    shimmed = dict(raw)
    shimmed["arm_keyed"] = False
    shimmed["vocabularies"] = [a.arm]
    out = Path(a.out)
    out.write_text(json.dumps(shimmed, indent=1))

    # ---- verify: exactly two keys added, every shared key equal ----
    back = json.loads(out.read_text())
    added = set(back) - set(raw)
    removed = set(raw) - set(back)
    if added != set(ADDED) or removed:
        raise SystemExit(f"[shim] key-set mismatch: added={sorted(added)} removed={sorted(removed)}")
    bad = [k for k in raw if back[k] != raw[k]]
    if bad:
        raise SystemExit(f"[shim] {len(bad)} shared keys changed value: {bad[:5]}")
    print(f"[shim] {src.name} -> {out}")
    print(f"[shim] added exactly {sorted(added)}; {len(raw)} pre-existing keys all compare equal")
    print(f"[shim] module {raw['name']!r}: {raw['n_regulons']} regulons / {raw['n_TFs']} TFs, "
          f"background {len(raw['background_regulons'])}, arm '{a.arm}'")

if __name__ == "__main__":
    main()
