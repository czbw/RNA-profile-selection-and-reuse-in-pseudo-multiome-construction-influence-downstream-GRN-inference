"""Random-gene-set null pool (null B), re-drawn size-matched to the given module."""
import argparse
import json

import numpy as np
import pandas as pd

import scoring_common as C

N_POOL_PER_SLOT = 200
MIN_DETECT_FRAC = 0.10
AUC_THRESHOLD = 0.05

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--arm", required=True, choices=["link", "repair"])
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    log = lambda m: print(f"[nullB:{a.tag}] {m}", flush=True)

    M = C.resolve_module(a.module, a.arm)
    OUT = C.out_dir(a.arm)
    log(f"{M['name']} {M['n_regulons']} regulons | arm vocabulary '{M['arm']}'")

    E = C.Expr(log=log)
    E.selftest_auc(AUC_THRESHOLD, log=log)

    rng = np.random.default_rng(C.SEED)          # shared scoring seed
    detected = E.X.index[(E.X > 0).mean(axis=1) >= MIN_DETECT_FRAC].astype(str)
    pool_genes = [g for g in detected if g in E.gset]
    log(f"NULL B universe: {len(pool_genes):,} genes detected in >={MIN_DETECT_FRAC:.0%} of samples")

    meta = {"provenance": C.provenance(M), "seed": C.SEED,
            "n_pool_per_slot": N_POOL_PER_SLOT, "min_detect_frac": MIN_DETECT_FRAC,
            "nullB_universe_size": len(pool_genes), "auc_threshold": AUC_THRESHOLD,
            "by_vocabulary": {}}

    for vocab in C.DEFS:
        tmap = M["targets"][vocab]
        sizes = [len([x for x in tmap[r] if x in E.gset]) for r in M["regulons"]]
        log(f"{vocab}: {len(sizes)} slots (one per module member), sizes {min(sizes)}-{max(sizes)}, "
            f"{N_POOL_PER_SLOT} sets/slot = {len(sizes)*N_POOL_PER_SLOT:,} random gene sets")
        pool = np.empty((len(E.samples), len(sizes), N_POOL_PER_SLOT), dtype=np.float32)
        for si, k in enumerate(sizes):
            for j in range(N_POOL_PER_SLOT):
                pick = rng.choice(len(pool_genes), size=min(k, len(pool_genes)), replace=False)
                pool[:, si, j] = E.auc_fast([pool_genes[i] for i in pick], AUC_THRESHOLD)
            if (si + 1) % 20 == 0:
                log(f"  slot {si+1}/{len(sizes)} done")
        f = OUT / f"nullB_pool_{a.tag}_{vocab}.npz"
        np.savez_compressed(f, pool=pool, sizes=np.array(sizes),
                            samples=np.array(E.samples, dtype=object))
        meta["by_vocabulary"][vocab] = {"n_slots": len(sizes), "slot_sizes_min": int(min(sizes)),
                                        "slot_sizes_max": int(max(sizes)),
                                        "slot_sizes_median": float(np.median(sizes)),
                                        "n_random_gene_sets": len(sizes) * N_POOL_PER_SLOT,
                                        "file": f.name}
        log(f"  saved {f.name} {pool.shape}")

    (OUT / f"nullB_pool_{a.tag}_meta.json").write_text(json.dumps(meta, indent=2))
    log("done")

if __name__ == "__main__":
    main()
