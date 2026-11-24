"""Build three-axis candidate sets and verify historical membership.

The 55-member corrected set is exploratory: its adoption followed inspection of
external results. Definition timestamps do not establish prospective validation."""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[2] / 'modules')          # <repo>/modules
REPO = ROOT.parent                                   # <repo>
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)

BUILD_DATE = "2026-08-20"
STAGES = ("Diagnosis", "Relapse")

# per-arm wiring. Everything arm-specific lives here and nowhere else.
ARMS = {
    "link": {
        "funnel": REPO / "longitudinal/published" / "out" / "nomination_funnel.csv",
        "committed": REPO / "adult/published" / "out" / "module_definition.json",
        "eregulons": lambda s: (REPO / "paml/malignant" / "out"
                                / f"eRegulons_{s}_link.tsv"),
        "vocab_tag": "link",
        "expect_universe": 431,
        "expect_old": (26, 14),
    },
    "repair": {
        "funnel": REPO / "longitudinal/corrected" / "out" / "nomination_funnel.csv",
        "committed": REPO / "adult/corrected" / "out" / "module_definition.json",
        "eregulons": lambda s: (REPO / "paml/malignant" / "out"
                                / f"eRegulons_{s}_repair.tsv"),
        "vocab_tag": "repair",
        "expect_universe": 201,
        "expect_old": (34, 21),
    },
}

CATALOGUE = REPO / "analysis/network_comparison" / "out" / "regulon_catalogue.csv"

# strip the SCENIC+ "_(123g)" / "_(45r)" suffix, exactly as the upstream builders do
regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

_LOG_LINES = []

def log(m):
    line = f"[modules] {m}"
    print(line, flush=True)
    _LOG_LINES.append(line)

def die(m):
    log(f"ABORT -- {m}")
    (OUT / "_log_m_00.txt").write_text("\n".join(_LOG_LINES) + "\n", encoding="utf-8",
                                       newline="\n")
    sys.exit(1)

# Load

def load_arms():
    A = {}
    for arm, cfg in ARMS.items():
        fun = pd.read_csv(cfg["funnel"], index_col=0)
        if fun.index.duplicated().any():
            die(f"{arm}: duplicate regulon names in {cfg['funnel']}")
        for c in ("ax1_directional", "ax2_selective", "ax3_celltype_robust", "ax4_tier1",
                  "survives_4", "TF", "tier", "direction", "sign_mode_lead"):
            if c not in fun.columns:
                die(f"{arm}: {cfg['funnel']} has no column `{c}`")
        for c in ("ax1_directional", "ax2_selective", "ax3_celltype_robust", "ax4_tier1",
                  "survives_4"):
            if fun[c].dtype != bool:
                die(f"{arm}: funnel column `{c}` is {fun[c].dtype}, not bool -- refusing to "
                    f"coerce, because a string 'False' is truthy and would silently pass rows")
        committed = json.loads(cfg["committed"].read_text(encoding="utf-8"))
        A[arm] = {"funnel": fun, "committed": committed, "cfg": cfg}
        log(f"{arm}: funnel {len(fun)} rows from {cfg['funnel'].relative_to(REPO)}; "
            f"committed definition '{committed['name']}' "
            f"{committed['n_regulons']}r/{committed['n_TFs']}TF, "
            f"background {len(committed['background_regulons'])}")
    return A

def check_entry_condition(A):
    facts = {}
    for arm, d in A.items():
        cfg = d["cfg"]
        per_stage = {}
        for s in STAGES:
            p = cfg["eregulons"](s)
            if not p.exists():
                die(f"{arm}: missing {p}")
            e = pd.read_csv(p, sep="\t", usecols=["eRegulon_name"])
            per_stage[s] = set(regkey(v) for v in e.eRegulon_name)
        inter = per_stage["Diagnosis"] & per_stage["Relapse"]
        idx = set(d["funnel"].index)
        if inter != idx:
            die(f"{arm}: funnel universe is NOT the Dx n Rel eRegulon-name intersection -- "
                f"{len(inter - idx)} in intersection only, {len(idx - inter)} in funnel only")
        if len(idx) != cfg["expect_universe"]:
            die(f"{arm}: universe is {len(idx)}, expected {cfg['expect_universe']}")
        facts[arm] = {"n_Diagnosis": len(per_stage["Diagnosis"]),
                      "n_Relapse": len(per_stage["Relapse"]),
                      "n_intersection": len(inter)}
        log(f"{arm}: entry condition OK -- Dx {len(per_stage['Diagnosis'])} n Rel "
            f"{len(per_stage['Relapse'])} = {len(inter)}, set-identical to the funnel index "
            f"(0 exceptions)")
    return facts

# Rules

def rule_T3(fun):
    """ax1 & ax2 & ax3 -- the new per-arm module."""
    m = fun.ax1_directional & fun.ax2_selective & fun.ax3_celltype_robust
    return list(fun.index[m])

def rule_old(fun):
    """ax1 & ax2 & ax3 & ax4_tier1 -- the retired rule, reproduced for the parity check."""
    m = fun.ax1_directional & fun.ax2_selective & fun.ax3_celltype_robust & fun.ax4_tier1
    return list(fun.index[m])

def tfs_of(fun, regs):
    return sorted(fun.loc[regs, "TF"].unique().tolist())

# PARITY CHECK

def parity_gate(A):
    ''
    log("")
    log("=" * 78)
    log("PARITY CHECK -- old rule  ax1_directional & ax2_selective & ax3_celltype_robust "
        "& ax4_tier1")
    log("=" * 78)
    failures = 0
    for arm, d in A.items():
        fun, com = d["funnel"], d["committed"]
        derived = rule_old(fun)
        derived_tfs = tfs_of(fun, derived)
        # also confirm the funnel's own precomputed survives_4 agrees with the expression
        s4 = list(fun.index[fun.survives_4])
        exp_r, exp_tf = d["cfg"]["expect_old"]

        checks = [
            ("survives_4 == ax1&2&3&4 (ordered)", s4 == derived,
             f"derived {len(derived)} vs survives_4 {len(s4)}"),
            ("n_regulons", com["n_regulons"] == len(derived),
             f"committed {com['n_regulons']} vs derived {len(derived)}"),
            ("n_TFs", com["n_TFs"] == len(derived_tfs),
             f"committed {com['n_TFs']} vs derived {len(derived_tfs)}"),
            ("regulons (ordered list)", list(com["regulons"]) == derived, None),
            ("regulons (as sets)", set(com["regulons"]) == set(derived), None),
            ("TFs (ordered list)", list(com["TFs"]) == derived_tfs, None),
            ("background_regulons == funnel universe (as sets)",
             set(com["background_regulons"]) == set(fun.index),
             f"committed {len(com['background_regulons'])} vs universe {len(fun.index)}"),
            (f"expected headline {exp_r}r/{exp_tf}TF",
             (len(derived), len(derived_tfs)) == (exp_r, exp_tf),
             f"got {len(derived)}r/{len(derived_tfs)}TF"),
        ]
        log(f"-- arm {arm} ({d['cfg']['committed'].relative_to(REPO)})")
        for label, ok, extra in checks:
            log(f"   [{'PASS' if ok else 'FAIL'}] {label}"
                + (f"  ({extra})" if extra and not ok else ""))
            if not ok:
                failures += 1
        # per-field diff on the membership sets
        for field, derived_set in (("regulons", set(derived)), ("TFs", set(derived_tfs))):
            com_set = set(com[field])
            only_com, only_der = sorted(com_set - derived_set), sorted(derived_set - com_set)
            if only_com or only_der:
                log(f"   DIFF {field}: only in committed = {only_com}")
                log(f"   DIFF {field}: only in derived   = {only_der}")
            else:
                log(f"   DIFF {field}: none -- {len(com_set)} members identical")
        # T5 is a membership field too; it is `survives_all5` in the funnel
        t5_derived = list(fun.index[fun.survives_all5]) if "survives_all5" in fun.columns else []
        t5_ok = list(com.get("T5_subset_reference_only", [])) == t5_derived
        log(f"   [{'PASS' if t5_ok else 'FAIL'}] T5_subset_reference_only (survives_all5): "
            f"committed {len(com.get('T5_subset_reference_only', []))} vs derived "
            f"{len(t5_derived)}")
        if not t5_ok:
            failures += 1

    # Known non-membership defect: reported, not counted as a parity failure.
    rp = A["repair"]["committed"]["provenance"]
    if "longitudinal/published" in rp:
        log("")
        log("NOTE (not a parity failure -- provenance is not a membership field): "
            "adult/corrected/out/module_definition.json carries the 34/21 CORRECTED core but "
            "its `provenance` names longitudinal/published/out/nomination_funnel.csv, which is the "
            "PUBLISHED arm's funnel and yields 26/14. Half-updated copy of "
            "adult/published/scripts/build_module.py. Not propagated by this script.")

    log("")
    log(f"PARITY CHECK: {'PASS' if failures == 0 else 'FAIL'} "
        f"({failures} membership-field failure(s))")
    log("=" * 78)
    return failures

# Carry-over with a hard error on any missing member

def carry_block(arm, com, regs):
    ''
    bg = list(com["background_regulons"])
    missing_bg = [r for r in regs if r not in set(bg)]
    if missing_bg:
        die(f"{arm}: {len(missing_bg)} module regulons absent from the carried-over "
            f"background_regulons: {missing_bg}")
    orient, targets = {}, {}
    for s in STAGES:
        o, t = com["orientation"][s], com["targets"][s]
        mo = [r for r in regs if r not in o]
        mt = [r for r in regs if r not in t]
        if mo:
            die(f"{arm}/{s}: {len(mo)} module regulons unsigned in the carried-over "
                f"orientation block: {mo}")
        if mt:
            die(f"{arm}/{s}: {len(mt)} module regulons have no target set in the carried-over "
                f"targets block: {mt}")
        orient[s] = o
        targets[s] = t
    log(f"{arm}: carry-over verified -- all {len(regs)} module regulons present in "
        f"background ({len(bg)}), orientation and targets under both vocabularies")
    return bg, orient, targets

def set_sizes(name, targets, regs):
    out = {}
    for s in STAGES:
        n = [len(targets[s][r]) for r in regs]
        out[s] = {f"{name}_median": float(np.median(n)),
                  f"{name}_min": int(min(n)),
                  f"{name}_max": int(max(n)),
                  f"{name}_n_regulons_measured": len(n)}
    return out

# Build

T5_NOTE = (
    "RETIRED AXIS -- empty by design. The old link definition's "
    "`T5_subset_reference_only` held the 10 regulons that passed the ENCODE/K562 ChIP axis "
    "(`survives_all5` = ax1&2&3&4&5); they had already FAILED the internal strictest null "
    "(null C). Two independent reasons not to carry that membership into a post-ax4 module: "
    "(i) the ChIP axis was dropped from the manuscript (2026-08-13), and it was never computable "
    "on the repair arm (no gate_malig_{stage}_repair_K562_w5000.csv exists, so "
    "nomination.py sets ax5=False for every repair regulon by construction); "
    "(ii) `survives_all5` is ax4-CONDITIONED, so importing it into a module defined to be "
    "ax4-free would smuggle the retired axis back in. The key is retained with an empty list so "
    "the JSON schema stays identical for existing readers. The historical 10 remain recoverable "
    "from adult/published/out/module_definition.json and from "
    "longitudinal/published/out/nomination_funnel.csv (`survives_all5`)."
)

def build(A, entry_facts):
    fun_l, fun_r = A["link"]["funnel"], A["repair"]["funnel"]
    t3 = {"link": rule_T3(fun_l), "repair": rule_T3(fun_r)}
    for arm in ("link", "repair"):
        log(f"{arm}: T3 = {len(t3[arm])} regulons / {len(tfs_of(A[arm]['funnel'], t3[arm]))} TFs")

    # ---- T4 = regulon-level intersection --------------------------------- #
    set_l, set_r = set(t3["link"]), set(t3["repair"])
    t4 = [r for r in t3["link"] if r in set_r]          # link funnel order, deterministic
    tf_t4 = tfs_of(fun_l, t4)
    # cross-arm TF-label sanity: the same regulon name must carry the same TF on both arms
    bad = [r for r in t4 if fun_l.loc[r, "TF"] != fun_r.loc[r, "TF"]]
    if bad:
        die(f"TF label disagrees across arms for shared regulons: {bad}")
    tf_l, tf_r = set(tfs_of(fun_l, t3["link"])), set(tfs_of(fun_r, t3["repair"]))
    tf_level = sorted(tf_l & tf_r)
    tf_only = sorted(set(tf_level) - set(tf_t4))
    log(f"T4_intersection (REGULON level) = {len(t4)} regulons carrying {len(tf_t4)} distinct TFs: "
        f"{tf_t4}")
    log(f"TF-level intersection of the two T3s = {len(tf_level)} TFs (reported statistic, NOT the "
        f"module's TF list)")
    log(f"shared only through non-matching sign modes ({len(tf_only)}): {tf_only}")
    for tf in tf_only:
        log(f"   {tf}: link {sorted(r for r in set_l if fun_l.loc[r, 'TF'] == tf)} | "
            f"repair {sorted(r for r in set_r if fun_r.loc[r, 'TF'] == tf)}")

    entry_txt = {
        arm: (f"entry condition: this arm's universe is the Diagnosis n Relapse eRegulon-name "
              f"intersection of the {ARMS[arm]['vocab_tag']} vocabulary "
              f"({entry_facts[arm]['n_Diagnosis']} n {entry_facts[arm]['n_Relapse']} = "
              f"{entry_facts[arm]['n_intersection']}), re-verified set-identical to the funnel "
              f"index on {BUILD_DATE} with 0 exceptions")
        for arm in ARMS}

    written = []

    # ---- the two per-arm T3 definitions ---------------------------------- #
    for arm, src_note in (
        ("link", "adult/published/out/module_definition.json (frozen 2026-08-06), which took "
                 "orientation from longitudinal/published/out/e2_paired_wilcoxon_reconciliation.csv "
                 "and target sets from paml/malignant/out/eRegulons_{stage}_link.tsv"),
        ("repair", "adult/corrected/out/module_definition.json (frozen 2026-08-06), which took "
                   "orientation from "
                   "longitudinal/corrected/out/e2_paired_wilcoxon_reconciliation.csv and target "
                   "sets from paml/malignant/out/eRegulons_{stage}_repair.tsv. NOTE: that "
                   "file's own `provenance` string wrongly names "
                   "longitudinal/published/out/nomination_funnel.csv; the membership it carries "
                   "(34/21) is the repair arm's, so the defect is in the string only. Not "
                   "propagated here."),
    ):
        name = f"T3_{arm}"
        fun, com = A[arm]["funnel"], A[arm]["committed"]
        regs = t3[arm]
        bg, orient, targets = carry_block(arm, com, regs)
        mod = {
            "name": name,
            "provenance": (
                f"{ARMS[arm]['funnel'].relative_to(REPO).as_posix()}; rule = "
                f"ax1_directional & ax2_selective & ax3_celltype_robust  (the fourth axis "
                f"ax4_tier1 is RETIRED as of {BUILD_DATE} and is NOT applied); "
                f"{entry_txt[arm]}; background_regulons, orientation and targets carried over "
                f"VERBATIM from {src_note}; membership fields and target_set_sizes recomputed by "
                f"modules/scripts/module_definitions.py on {BUILD_DATE}."),
            "frozen": (f"{BUILD_DATE}, frozen before any downstream statistic (replication "
                       f"scoring, external-cohort scoring, figures, manuscript numbers) was "
                       f"recomputed under the retired-ax4 redefinition. Supersedes the "
                       f"2026-08-06 freeze of {ARMS[arm]['committed'].relative_to(REPO).as_posix()}"
                       f" ('T4_tier1_core', {com['n_regulons']}r/{com['n_TFs']}TF), which is "
                       f"retired as a module but retained unmodified for audit."),
            "membership_level": "regulon",
            "arm_keyed": False,
            "vocabularies": [ARMS[arm]["vocab_tag"]],
            "n_regulons": len(regs),
            "n_TFs": len(tfs_of(fun, regs)),
            "TFs": tfs_of(fun, regs),
            "regulons": regs,
            "T5_subset_reference_only": [],
            "T5_note": T5_NOTE,
            "background_regulons": bg,
            "orientation": orient,
            "targets": targets,
            "target_set_sizes": set_sizes(name, targets, regs),
        }
        p = OUT / f"module_definition_{name}.json"
        p.write_text(json.dumps(mod, indent=1), encoding="utf-8", newline="\n")
        written.append(p)
        log(f"wrote {p.relative_to(REPO)} -- {mod['n_regulons']} regulons / {mod['n_TFs']} TFs, "
            f"background {len(bg)}")
        for s in STAGES:
            v = mod["target_set_sizes"][s]
            log(f"   {s} target-set size: median {v[f'{name}_median']:.0f} "
                f"(range {v[f'{name}_min']}-{v[f'{name}_max']})")

    # ---- T4_intersection, arm-keyed -------------------------------------- #
    # The large blocks are keyed by arm first (deliberate): code that uses the flat
    # `mod["orientation"]["Diagnosis"]` idiom raises KeyError instead of silently scoring
    # one arm's vocabulary as the other's.
    bgs, orients, targs, sizes = {}, {}, {}, {}
    for arm in ("link", "repair"):
        b, o, t = carry_block(arm, A[arm]["committed"], t4)
        bgs[arm], orients[arm], targs[arm] = b, o, t
        sizes[arm] = set_sizes("T4_intersection", t, t4)
    t4mod = {
        "name": "T4_intersection",
        "provenance": (
            "REGULON-LEVEL intersection of T3_link and T3_repair. T3_link = "
            "longitudinal/published/out/nomination_funnel.csv with rule ax1_directional & "
            "ax2_selective & ax3_celltype_robust; T3_repair = "
            "longitudinal/corrected/out/nomination_funnel.csv with the same rule; the retired "
            "ax4_tier1 is not applied on either arm. "
            + entry_txt["link"] + "; " + entry_txt["repair"] + ". "
            "background_regulons, orientation and targets carried over VERBATIM from "
            "adult/published/out/module_definition.json (link vocabulary) and "
            "adult/corrected/out/module_definition.json (repair vocabulary), kept SEPARATE and "
            "labelled by arm -- never averaged or merged. Membership fields and target_set_sizes "
            f"recomputed by modules/scripts/module_definitions.py on {BUILD_DATE}."),
        "frozen": (f"{BUILD_DATE}, frozen before any downstream statistic (replication scoring, "
                   f"external-cohort scoring, figures, manuscript numbers) was recomputed under "
                   f"the retired-ax4 redefinition. Replaces, and does not extend, the retired "
                   f"'T4_tier1_core' modules in adult/published/out/module_definition.json (26r/14TF) "
                   f"and adult/corrected/out/module_definition.json (34r/21TF), both retained "
                   f"unmodified for audit."),
        "membership_level": "regulon",
        "arm_keyed": True,
        "vocabularies": ["link", "repair"],
        "n_regulons": len(t4),
        "n_TFs": len(tf_t4),
        "TFs": tf_t4,
        "regulons": t4,
        "tf_level_note": (
            f"LEVEL MATTERS -- DO NOT CONFLATE. This module is defined at the REGULON level: "
            f"{len(t4)} regulons, which happen to carry {len(tf_t4)} distinct TFs. `TFs` above is "
            f"that {len(tf_t4)}-TF list and is the module's TF list. Separately, the TF-LEVEL "
            f"intersection of the two T3 TF lists is {len(tf_level)} TFs; it is recorded below as "
            f"a reported statistic ONLY and is NOT the module. The two differ because "
            f"{len(tf_only)} TFs ({', '.join(tf_only)}) appear in both T3s only through "
            f"NON-MATCHING sign modes, i.e. no single regulon name of theirs is in both T3s."),
        "n_TFs_tf_level_intersection_reference_only": len(tf_level),
        "tf_level_intersection_reference_only": tf_level,
        "tf_level_only_via_nonmatching_sign_modes": tf_only,
        "T5_subset_reference_only": [],
        "T5_note": T5_NOTE,
        "background_regulons": bgs,
        "orientation": orients,
        "targets": targs,
        "target_set_sizes": sizes,
    }
    p = OUT / "shared_module.json"
    p.write_text(json.dumps(t4mod, indent=1), encoding="utf-8", newline="\n")
    written.append(p)
    log(f"wrote {p.relative_to(REPO)} -- {len(t4)} regulons / {len(tf_t4)} TFs (regulon level); "
        f"TF-level intersection {len(tf_level)} recorded separately; backgrounds "
        f"link {len(bgs['link'])} / repair {len(bgs['repair'])}")
    for arm in ("link", "repair"):
        for s in STAGES:
            v = sizes[arm][s]
            log(f"   {arm}/{s} target-set size: median "
                f"{v['T4_intersection_median']:.0f} (range {v['T4_intersection_min']}-"
                f"{v['T4_intersection_max']})")

    # ---- crosswalk -------------------------------------------------------- #
    p = crosswalk(A, t3, t4)
    written.append(p)
    return t3, t4, tf_t4, tf_level, tf_only, written

# Crosswalk

def crosswalk(A, t3, t4):
    fun_l, fun_r = A["link"]["funnel"], A["repair"]["funnel"]
    old_l = rule_old(fun_l)
    old_r = rule_old(fun_r)

    cat = pd.read_csv(CATALOGUE)
    cat_tier = dict(zip(cat.regulon, cat.tier))
    # short tier, exactly as both nomination scripts derive it
    cat_short = {r: str(t).split()[0] for r, t in cat_tier.items()}

    members = {
        "in_old_T4_link_26": set(old_l),
        "in_old_core_repair_34": set(old_r),
        "in_T3_link_98": set(t3["link"]),
        "in_T3_repair_55": set(t3["repair"]),
        "in_T4_new_13": set(t4),
    }
    union = sorted(set().union(*members.values()))

    rows = []
    for r in union:
        in_l, in_r = r in fun_l.index, r in fun_r.index
        tf = fun_l.loc[r, "TF"] if in_l else fun_r.loc[r, "TF"]
        m = re.search(r"_([+-]/[+-])$", str(r))
        sign_mode = m.group(1) if m else "?"
        row = {
            "regulon": r,
            "TF": tf,
            "sign_mode": sign_mode,
            "sign_mode_lead": sign_mode[0] if sign_mode != "?" else "?",
            "in_universe_link": in_l,
            "in_universe_repair": in_r,
            "tier_link": fun_l.loc[r, "tier"] if in_l else "",
            "tier_repair": fun_r.loc[r, "tier"] if in_r else "",
            "catalogue_tier_raw": cat_tier.get(r, "not in catalogue"),
            "direction_link": fun_l.loc[r, "direction"] if in_l else "",
            "direction_repair": fun_r.loc[r, "direction"] if in_r else "",
        }
        for k, s in members.items():
            row[k] = r in s
        rows.append(row)
    df = pd.DataFrame(rows)

    # -- invariants asserted -------------------------- #
    # (a) the funnels' `tier` must be the catalogue short tier, default "none"
    for arm, fun in (("link", fun_l), ("repair", fun_r)):
        exp = pd.Series([cat_short.get(i, "none") for i in fun.index], index=fun.index)
        bad = fun.index[fun["tier"].astype(str) != exp]
        if len(bad):
            die(f"{arm}: funnel `tier` disagrees with the catalogue for {len(bad)} regulons "
                f"(e.g. {list(bad[:5])})")
    # (b) sign_mode parsed from the name must match the funnel's own sign_mode_lead
    for arm, fun in (("link", fun_l), ("repair", fun_r)):
        sub = df[df[f"in_universe_{arm}"]]
        got = fun.loc[sub.regulon, "sign_mode_lead"].astype(str).values
        if not (got == sub.sign_mode_lead.values).all():
            die(f"{arm}: sign_mode parsed from the regulon name disagrees with the funnel's "
                f"sign_mode_lead column")
    # (c) the redefinition is a relaxation: dropping ax4 can only add members
    if not set(old_l) <= set(t3["link"]):
        die(f"old link T4 is NOT a subset of T3_link -- {sorted(set(old_l) - set(t3['link']))}")
    if not set(old_r) <= set(t3["repair"]):
        die(f"old repair core is NOT a subset of T3_repair -- "
            f"{sorted(set(old_r) - set(t3['repair']))}")
    log("crosswalk invariants OK: funnel tier == catalogue tier on both arms; parsed sign mode "
        "== funnel sign_mode_lead on both arms; both old modules are subsets of their arm's T3")

    p = OUT / "module_membership_crosswalk.csv"
    df.to_csv(p, index=False, lineterminator="\n")
    log("")
    log(f"wrote {p.relative_to(REPO)} -- {len(df)} rows (union of the five membership sets)")
    log("-- crosswalk headline counts --")
    log(f"   old link T4        : {len(old_l)}")
    log(f"   old repair core    : {len(old_r)}")
    log(f"   T3_link            : {len(t3['link'])}  "
        f"(contains {len(set(old_l) & set(t3['link']))}/{len(old_l)} of the old link T4; "
        f"gains {len(set(t3['link']) - set(old_l))} new members)")
    log(f"   T3_repair          : {len(t3['repair'])}  "
        f"(contains {len(set(old_r) & set(t3['repair']))}/{len(old_r)} of the old repair core; "
        f"gains {len(set(t3['repair']) - set(old_r))} new members)")
    log(f"   T4_new             : {len(t4)}  "
        f"(contains {len(set(old_l) & set(t4))}/{len(old_l)} of the old link T4 and "
        f"{len(set(old_r) & set(t4))}/{len(old_r)} of the old repair core; "
        f"{len(set(t4) - set(old_l) - set(old_r))} members in neither old module)")
    log(f"   old link T4 also in T3_repair : {len(set(old_l) & set(t3['repair']))}")
    log(f"   old link T4 evaluable on repair arm (in its universe) : "
        f"{len(set(old_l) & set(fun_r.index))}")
    log(f"   old repair core also in T3_link : {len(set(old_r) & set(t3['link']))}")
    log(f"   union size         : {len(union)}")
    return p

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("parity", "build"), default="build")
    args = ap.parse_args()

    A = load_arms()
    entry_facts = check_entry_condition(A)
    failures = parity_gate(A)
    if failures:
        die(f"parity check failed with {failures} membership-field difference(s) -- refusing to "
            f"emit any new definition. Nothing was written.")
    if args.mode == "parity":
        log("mode=parity -- check passed, nothing written.")
    else:
        build(A, entry_facts)
    (OUT / "_log_m_00.txt").write_text("\n".join(_LOG_LINES) + "\n", encoding="utf-8",
                                       newline="\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())
