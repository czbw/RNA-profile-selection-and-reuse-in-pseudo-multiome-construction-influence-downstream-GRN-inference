#!/usr/bin/env bash
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus
command -v scenicplus >/dev/null 2>&1 || {
  echo "[seedrep] ABORT — scenicplus is not on PATH after activating the env. CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS:-unset}" >&2
  exit 1; }
echo ">> env: $(command -v python) | scenicplus: $(command -v scenicplus)"

STAGE="${AMLGRN_STAGE:-pbmc10k}"
ARM="${AMLGRN_ARM:-truth}"
REP="${AMLGRN_REP:-rep2}"

ROOT="${AMLGRN_SRC_BUNDLE:?export AMLGRN_SRC_BUNDLE (the benchmarks/pbmc bundle)}"
# Note: avoid apostrophes inside ${VAR:?msg} messages; an unmatched quote makes `bash -n` report
# "unexpected EOF" at the last line of the file.
REPOUT="${AMLGRN_REP_OUT:?export AMLGRN_REP_OUT (the out directory of this bundle)}"
DB="${AMLGRN_DB:?export AMLGRN_DB}"
NCPU="${AMLGRN_NCPU:-16}"

OUT="$ROOT/out"
SRC="${STAGE}_${ARM}"
TAG="${STAGE}_${ARM}_${REP}"
TMP="${AMLGRN_REP_TMP:?export AMLGRN_REP_TMP}"
RANK="$DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"

mkdir -p "$TMP" "$REPOUT"

die()  { echo "[seedrep] ABORT — $*" >&2; exit 1; }
step() { echo; echo "===== [$TAG] $* ====="; date '+%Y-%m-%d %H:%M:%S'; }

echo "=================================================================="
echo " PBMC RE-RUN REPLICATE   $(date)"
echo " source arm : $SRC        (read-only)"
echo " replicate  : $TAG        -> $REPOUT"
echo " n_cpu      : $NCPU       tmp: $TMP"
echo "=================================================================="

# --- the reused inputs: checked, never rebuilt --------------------------------------------------
MUD="$OUT/mudata_${SRC}.h5mu"
SS="$OUT/search_space_${SRC}.tsv"
TFN="$OUT/tf_names_${SRC}.txt"
CIS="$OUT/cistromes_direct_${SRC}.h5ad"
for p in "$MUD" "$SS" "$TFN" "$CIS" "$RANK"; do
  [ -s "$p" ] || die "missing/empty reused input: $p"
done
echo ">> reused inputs present (mudata, search_space, tf_names, cistromes, ranking db)"

# The chain must not have changed since the source arm ran (the same check run_grn.sh makes).
[ -f "$ROOT/scripts/chain_guard.sh" ] \
  || die "chain_guard.sh not found at $ROOT/scripts/ — refusing to run a replicate that cannot prove the chain is unchanged"

# chain_guard resolves the chain as "$AMLGRN_ROOT/out", so AMLGRN_ROOT must be the source bundle.
# Set it explicitly so the script also runs standalone.
AMLGRN_ROOT="$ROOT" bash "$ROOT/scripts/chain_guard.sh" verify "$STAGE" "seedrep:$ARM"
cg=$?
# Distinguish "the guard could not run" from "the guard ran and the chain changed", so that a
# failure to check is not reported as a mismatch. chain_guard exit codes: 7 bad/missing AMLGRN_ROOT
# or out/, 3 chain artifact missing, 5 no reference recorded, 6 genuine mismatch, 2 usage.
case "$cg" in
  0) echo ">> chain_guard verify OK — the chain is byte-identical to the one $SRC consumed" ;;
  6) die "chain_guard reports a MISMATCH: the chromatin chain has MOVED since $SRC was built. A replicate against a different chain measures nothing. Do not re-run until this is explained." ;;
  7) die "chain_guard could not run: AMLGRN_ROOT is unset, empty, or not a directory (it must be the benchmarks/pbmc bundle). THE CHAIN WAS NOT CHECKED — this is not evidence that it moved." ;;
  3) die "chain_guard could not run: a chain artifact (cistopic_obj / region_sets / menr) is missing or empty under $ROOT/out. THE CHAIN WAS NOT COMPARED." ;;
  5) die "chain_guard could not run: no atac_chain_${STAGE}.sha256 reference exists, so there is nothing to verify against. THE CHAIN WAS NOT COMPARED." ;;
  *) die "chain_guard exited $cg (usage or unexpected error). THE CHAIN WAS NOT VERIFIED; refusing to proceed." ;;
esac

# --- do not overwrite a completed replicate -----------------------------------------------------
if [ -s "$REPOUT/eRegulons_${TAG}.tsv" ] && [ -s "$REPOUT/eRegulon_AUC_${TAG}.h5mu" ]; then
  echo ">> replicate already complete — nothing to do (delete the outputs to force a re-run)"
  exit 0
fi

# --- 1) region_to_gene --------------------------------------------------------------------------
step "region_to_gene (GBM — stochastic)"
[ -s "$REPOUT/r2g_${TAG}.tsv" ] || scenicplus grn_inference region_to_gene \
  --multiome_mudata_fname "$MUD" --search_space_fname "$SS" \
  --temp_dir "$TMP" --out_region_to_gene_adjacencies "$REPOUT/r2g_${TAG}.tsv" \
  --importance_scoring_method GBM --correlation_scoring_method SR \
  --n_cpu "$NCPU" || die "region_to_gene"

# --- 2) TF_to_gene ------------------------------------------------------------------------------
step "TF_to_gene (GBM — stochastic)"
[ -s "$REPOUT/tf2g_${TAG}.tsv" ] || scenicplus grn_inference TF_to_gene \
  --multiome_mudata_fname "$MUD" --tf_names "$TFN" \
  --temp_dir "$TMP" --out_tf_to_gene_adjacencies "$REPOUT/tf2g_${TAG}.tsv" \
  --method GBM --n_cpu "$NCPU" || die "TF_to_gene"

step "eGRN"
if [ ! -s "$REPOUT/eRegulons_${TAG}.tsv" ]; then
  scenicplus grn_inference eGRN \
    --TF_to_gene_adj_fname "$REPOUT/tf2g_${TAG}.tsv" \
    --region_to_gene_adj_fname "$REPOUT/r2g_${TAG}.tsv" \
    --cistromes_fname "$CIS" \
    --ranking_db_fname "$RANK" \
    --eRegulon_out_fname "$REPOUT/eRegulons_${TAG}.tsv" \
    --temp_dir "$TMP"
  rc=$?
  # The source arm `truth` produced an eGRN, so a failure here is a defect, not the structurally
  # impossible result that egrn_guard.py certifies for pairings that destroy RNA-ATAC covariation.
  if [ "$rc" -ne 0 ]; then
    die "eGRN returned $rc. The source arm '$ARM' HAS an eRegulon table, so an empty adjacency is not a credible explanation and this is not routed to egrn_guard.py."
  fi
fi

# --- 4) AUCell ----------------------------------------------------------------------------------
step "AUCell"
[ -s "$REPOUT/eRegulon_AUC_${TAG}.h5mu" ] || scenicplus grn_inference AUCell \
  --eRegulon_fname "$REPOUT/eRegulons_${TAG}.tsv" \
  --multiome_mudata_fname "$MUD" \
  --aucell_out_fname "$REPOUT/eRegulon_AUC_${TAG}.h5mu" --n_cpu 1 || die "AUCell"

echo
echo "=== [$TAG] REPLICATE COMPLETE ==="
wc -l "$REPOUT/eRegulons_${TAG}.tsv" 2>/dev/null
ls -lh "$REPOUT"/{r2g,tf2g}_"${TAG}".tsv "$REPOUT/eRegulons_${TAG}.tsv" "$REPOUT/eRegulon_AUC_${TAG}.h5mu" 2>/dev/null

step "computing the reproducibility ceiling"
python "$(dirname "$0")/compare_replicates.py" \
  --src "$OUT" --rep "$REPOUT" --stage "$STAGE" --arm "$ARM" --reptag "$REP" \
  || die "compare_replicates.py"

echo
echo "=== [$TAG] CEILING COMPUTED — see $REPOUT/seedrep_ceiling.json ==="
