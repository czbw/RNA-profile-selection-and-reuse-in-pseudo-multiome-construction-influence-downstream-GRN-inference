#!/usr/bin/env bash
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -uo pipefail

ROOT="${AMLGRN_ROOT_DIR:-$AMLGRN_PROJECT}"
BUNDLE="${AMLGRN_BUNDLE:-$ROOT/benchmarks/bmmc}"
SIF="$ROOT/miniforge.sif"
STAGE="${AMLGRN_STAGE:-bmmc}"
SEED="${AMLGRN_SEED:-555}"
ARMS=("$@")
[ "${#ARMS[@]}" -gt 0 ] || ARMS=(pseudomultiome pseudomultiome_k10)

PASS=0; FAILN=0
ok()  { PASS=$((PASS+1));  printf '  ok    %s\n' "$*"; }
bad() { FAILN=$((FAILN+1)); printf '  FAIL  %s\n' "$*" >&2; }

echo "=================================================================="
echo " BMMC metacell preflight   $(date)"
echo " bundle : $BUNDLE"
echo " stage  : $STAGE   seed: $SEED"
echo " arms   : ${ARMS[*]}"
echo "=================================================================="

[ -d "$BUNDLE" ] || { echo "FATAL: no bundle at $BUNDLE"; exit 2; }
cd "$BUNDLE" || { echo "FATAL: cannot cd $BUNDLE"; exit 2; }

echo; echo "---- [1] the k each arm resolves to, parsed the way run_grn_bmmc.sh parses it ----"
# run_grn_bmmc.sh reads the metacell size from the arm name, so two metacell arms cannot silently
# share a setting. Same parse here, so a mistyped arm name fails now rather than hours into the run.
declare -a KS=()
for ARM in "${ARMS[@]}"; do
  case "$ARM" in
    pseudomultiome)    K=5 ;;
    pseudomultiome_k*) K="${ARM#pseudomultiome_k}" ;;
    *) bad "'$ARM' is not a metacell arm -- this preflight covers pseudomultiome / pseudomultiome_k<N> only"; continue ;;
  esac
  case "$K" in ''|*[!0-9]*) bad "cannot parse metacell size from arm '$ARM'"; continue ;; esac
  if [ "$K" -lt 2 ]; then bad "arm '$ARM' -> k=$K; run_grn_bmmc.sh requires k >= 2"; continue; fi
  KS+=("$K")
  ok "arm '$ARM' -> k=$K"
done

echo; echo "---- [2] scripts these arms execute ----"
for f in "$BUNDLE/submit_bmmc.sh" \
         "$BUNDLE/submit_bmmc_arm.sh" \
         "$BUNDLE/run_bmmc.sh" \
         "$BUNDLE/scripts/run_grn_bmmc.sh" \
         "$BUNDLE/scripts/chain_guard.sh" \
         "$BUNDLE/scripts/egrn_guard.py" \
         "$BUNDLE/scripts/metacell_floor.py" \
         "$BUNDLE/scripts/config.yaml"; do
  [ -s "$f" ] && ok "present: ${f#$BUNDLE/}" || bad "MISSING: $f"
done
# submit_bmmc.sh asserts these three files unconditionally and exits 2 without them, for every arm
# set including metacell arms. A metacell arm does not use them, but the submission still requires them.
for f in "$BUNDLE/scripts/atac_chain.py" \
         "$BUNDLE/scripts/region_sets.py" \
         "$BUNDLE/scripts/build_arm_rna.py"; do
  [ -s "$f" ] && ok "present (submit_bmmc.sh gates on it): ${f#$BUNDLE/}" \
              || bad "MISSING: $f -- submit_bmmc.sh exits 2 before entering the container"
done

echo; echo "---- [3] the floor assertion is WIRED into run_grn_bmmc.sh ----"
if grep -Fq "metacell_floor.py" "$BUNDLE/scripts/run_grn_bmmc.sh"; then
  ok "run_grn_bmmc.sh calls metacell_floor.py"
else
  bad "run_grn_bmmc.sh never mentions metacell_floor.py -- the metacell floor would not be asserted at the point of consumption"
fi

echo; echo "---- [4] 'no assignment table' -- by EXECUTING the dispatch, not by grepping it ----"
NA_SRC="$(sed -n '/^needs_assignment() {/,/^}/p' "$BUNDLE/run_bmmc.sh")"
IM_SRC="$(sed -n '/^is_metacell() {/,/^}/p' "$BUNDLE/run_bmmc.sh")"
if [ -z "$NA_SRC" ] || [ -z "$IM_SRC" ]; then
  bad "could not extract needs_assignment()/is_metacell() from run_bmmc.sh -- refusing to report a pass on an empty extraction"
elif ! printf '%s' "$NA_SRC" | grep -Fq 'return'; then
  bad "the extracted needs_assignment() has no 'return' -- the extraction is wrong, so any verdict from it would be meaningless"
else
  eval "$NA_SRC"
  eval "$IM_SRC"
  for ARM in "${ARMS[@]}"; do
    if is_metacell "$ARM"; then ok "run_bmmc.sh is_metacell('$ARM') -> true"
    else bad "run_bmmc.sh is_metacell('$ARM') -> FALSE; it would be treated as a matching arm"; fi
    if needs_assignment "$ARM"; then
      bad "run_bmmc.sh needs_assignment('$ARM') -> TRUE; it would demand uploads/arm_${ARM}_${STAGE}_seed${SEED}_full.csv"
    else
      ok "run_bmmc.sh needs_assignment('$ARM') -> false (exempt by design)"
    fi
  done
fi

SB_SRC="$(sed -n '/^SEED_USED=/,/^done$/p' "$BUNDLE/submit_bmmc.sh")"
if [ -z "$SB_SRC" ] || ! printf '%s' "$SB_SRC" | grep -Fq 'MISSING assignment table'; then
  bad "could not extract submit_bmmc.sh's assignment-table loop (or the extraction lost its error branch) -- not reporting a pass on it"
else
  if ( BUNDLE="$BUNDLE"; AMLGRN_SEED="$SEED"; set -- "${ARMS[@]}"; eval "$SB_SRC" ) >/dev/null 2>&1; then
    ok "submit_bmmc.sh's preflight loop passes for ${ARMS[*]} with no assignment table on disk"
  else
    bad "submit_bmmc.sh's preflight loop DEMANDS an assignment table for one of: ${ARMS[*]}"
  fi
fi

echo; echo "---- [5] the frozen chromatin chain (built by phase [2]; these arms only consume it) ----"
for f in "$BUNDLE/out/cistopic_obj_${STAGE}.pkl" \
         "$BUNDLE/out/menr_cistarget_${STAGE}.hdf5" \
         "$BUNDLE/out/atac_chain_${STAGE}.sha256" \
         "$BUNDLE/out/cistopic_cell_names_${STAGE}.txt"; do
  [ -s "$f" ] && ok "present: ${f#$BUNDLE/} ($(du -h "$f" 2>/dev/null | cut -f1))" \
              || bad "MISSING: $f -- run 'AMLGRN_STOP_AFTER=chain sbatch submit_bmmc.sh' first"
done
[ -d "$BUNDLE/out/region_sets_${STAGE}" ] && ok "present: out/region_sets_${STAGE}/" \
  || bad "MISSING: $BUNDLE/out/region_sets_${STAGE} -- run_bmmc.sh phase [2b] builds it"

echo; echo "---- [6] the metacell GEX ----"
GEX="$BUNDLE/out/native_rna_${STAGE}.metacell.h5ad"
if [ -s "$GEX" ]; then
  ok "present: out/$(basename "$GEX") ($(du -h "$GEX" | cut -f1))"
else
  bad "MISSING: $GEX -- built by scripts/native_rna.py + scripts/filter_gex_for_metacells.py and scp'd into out/; see scripts/UPLOAD_MANIFEST.md sec.3"
fi
[ -s "$BUNDLE/uploads/b00_cell_metadata.csv" ] && ok "present: uploads/b00_cell_metadata.csv" \
  || bad "MISSING: $BUNDLE/uploads/b00_cell_metadata.csv"
[ -s "$BUNDLE/data/multiome_BMMC.h5ad" ] && ok "present: data/multiome_BMMC.h5ad" \
  || bad "MISSING: $BUNDLE/data/multiome_BMMC.h5ad -- run_bmmc.sh phase [1] asserts it for every arm"

echo; echo "---- [7] the GRCh38 cisTarget databases (reused from Arm A, never staged here) ----"
for f in "$ROOT/databases/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather" \
         "$ROOT/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
         "$ROOT/databases/hg38_genome_annotation.tsv" \
         "$ROOT/databases/hg38_chromsizes_3col.tsv"; do
  [ -s "$f" ] && ok "present: ${f#$ROOT/}" || bad "MISSING DATABASE: $f"
done

echo; echo "---- [8] CR bytes and bash -n on what these arms run ----"
for f in "$BUNDLE/submit_bmmc.sh" "$BUNDLE/submit_bmmc_arm.sh" "$BUNDLE/run_bmmc.sh" \
         "$BUNDLE/preflight_metacell_bmmc.sh" \
         "$BUNDLE/scripts/run_grn_bmmc.sh" "$BUNDLE/scripts/chain_guard.sh" \
         "$BUNDLE/scripts/metacell_floor.py"; do
  [ -e "$f" ] || continue
  n="$(tr -cd '\r' < "$f" | wc -c | tr -d ' ')"
  [ "$n" -eq 0 ] && ok "LF-only: ${f#$BUNDLE/}" \
                 || bad "$n CR byte(s) in ${f#$BUNDLE/} -- fix with: sed -i 's/\r\$//' $f"
done
for f in "$BUNDLE/submit_bmmc.sh" "$BUNDLE/submit_bmmc_arm.sh" "$BUNDLE/run_bmmc.sh" \
         "$BUNDLE/preflight_metacell_bmmc.sh" "$BUNDLE/scripts/run_grn_bmmc.sh"; do
  [ -e "$f" ] || continue
  bash -n "$f" 2>/dev/null && ok "bash -n: ${f#$BUNDLE/}" || bad "bash -n FAILED: ${f#$BUNDLE/}"
done

echo; echo "---- [9] the floor + cell-id assertion, run for real, once per k ----"
# This check decides whether prepare_GEX_ACC will survive. It needs anndata, so it runs inside the
# same container image the job uses, on the login node, for a few seconds.
# AMLGRN_PREFLIGHT_NO_CONTAINER=1 skips it and prints the command instead.
if [ ! -s "$GEX" ]; then
  bad "[9] skipped: the metacell GEX is missing (see [6]) -- the floor cannot be asserted against a file that is not there"
elif [ "${AMLGRN_PREFLIGHT_NO_CONTAINER:-0}" = "1" ]; then
  bad "[9] SKIPPED by AMLGRN_PREFLIGHT_NO_CONTAINER=1 -- the floor and the cell-id convention are UNCHECKED. Run this before submitting."
  for K in "${KS[@]}"; do
    echo "        apptainer exec --bind $ROOT $SIF bash -c 'source \$(conda info --base)/etc/profile.d/conda.sh && conda activate scenicplus && cd $BUNDLE/scripts && python metacell_floor.py --gex $GEX --stage $STAGE --k $K --analysis-set $BUNDLE/out/cistopic_cell_names_${STAGE}.txt --metadata $BUNDLE/uploads/b00_cell_metadata.csv --require-analysis-set'"
  done
elif [ ! -s "$SIF" ]; then
  bad "[9] cannot run: $SIF is missing, so the floor and the cell-id convention are UNCHECKED"
else
  export APPTAINERENV_CONDA_ENVS_DIRS="$ROOT/envs"
  export APPTAINERENV_AMLGRN_ROOT="$BUNDLE"
  export APPTAINERENV_AMLGRN_STAGE="$STAGE"
  export APPTAINERENV_AMLGRN_CONFIG="config.yaml"
  for K in "${KS[@]}"; do
    echo "  --- metacell_floor.py --k $K ---"
    apptainer exec --bind "$ROOT" "$SIF" bash -c \
      "source \$(conda info --base)/etc/profile.d/conda.sh && conda activate scenicplus && cd '$BUNDLE/scripts' && python metacell_floor.py --gex '$GEX' --stage '$STAGE' --k $K --analysis-set '$BUNDLE/out/cistopic_cell_names_${STAGE}.txt' --metadata '$BUNDLE/uploads/b00_cell_metadata.csv' --require-analysis-set"
    rc=$?
    case "$rc" in
      0) ok "k=$K: every shared cell type clears the floor AND the GEX cell ids are the frozen analysis set" ;;
      3) bad "k=$K: FLOOR VIOLATED -- prepare_GEX_ACC would raise 'Sample larger than population'. Rebuild the metacell GEX with filter_gex_for_metacells.py --k $K, or drop this arm." ;;
      4) bad "k=$K: CELL-ID CONVENTION VIOLATED -- the GEX ids are not the frozen chain's analysis set " ;;
      2) bad "k=$K: b_12 could not run (missing/unusable input, or --require-analysis-set with no sidecar)" ;;
      *) bad "k=$K: metacell_floor.py exited $rc" ;;
    esac
  done
fi

echo
echo "=================================================================="
if [ "$FAILN" -ne 0 ]; then
  echo " METACELL PREFLIGHT FAILED: $FAILN problem(s), $PASS check(s) passed"
  echo " Nothing was submitted and nothing was written."
  echo "=================================================================="
  exit 1
fi
echo " METACELL PREFLIGHT OK: $PASS check(s) passed"
echo
echo " Submit with (GRN-half resources; --propagate=NONE and the mail directives are in"
echo " submit_bmmc_arm.sh's #SBATCH header):"
echo
echo "   cd $BUNDLE && sbatch --propagate=NONE submit_bmmc_arm.sh ${ARMS[*]}"
echo "=================================================================="
