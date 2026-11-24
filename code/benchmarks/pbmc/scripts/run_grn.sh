#!/usr/bin/env bash
set -uo pipefail
STAGE="${1:?stage}"; PAIRING="${2:?truth|link|pseudomultiome}"
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT}"
# [PBMC] env override for DB: AMLGRN_ROOT is the bundle, but the feathers stay one level up.
OUT="$ROOT/out"; DB="${AMLGRN_DB:-$ROOT/databases}"; TMP="$ROOT/work/grn_tmp_${STAGE}_${PAIRING}"
mkdir -p "$TMP"
RANK="$DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
MOTIF_TBL="$DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TAG="${STAGE}_${PAIRING}"
KMETA=""          # non-empty ONLY for metacell (non-multiome) arms; also selects the branch below
step() { echo "===== [$TAG] $* ====="; }
die()  { echo "FAILED [$TAG]: $*"; exit 1; }

case "$PAIRING" in
  truth|link|link_1to1|link_subcluster_1to1|rand_subcluster_matched|rand_subcluster_1to1|rand_lineage_matched|\
rand_compartment_matched|rand_global)
                  GEX="$OUT/paired_rna_${STAGE}.${PAIRING}.h5ad" ;;
  pseudomultiome|pseudomultiome_k*)                          # Option B (non-multiome metacells)
    case "$PAIRING" in
      pseudomultiome)    KMETA=5 ;;
      pseudomultiome_k*) KMETA="${PAIRING#pseudomultiome_k}" ;;
    esac
    case "$KMETA" in ''|*[!0-9]*) die "cannot parse metacell size from arm '$PAIRING'" ;; esac
    [ "$KMETA" -ge 2 ] || die "metacell size must be >= 2 (got $KMETA)"

    # The metacell GEX is built once per stage and does not depend on k: it is a cell-type filter
    # with floor max(min_cells_per_celltype=25, k), i.e. 25 for every k <= 25. It is shared by all
    # metacell arms and not rebuilt here; a missing file is an error.
    GEX="$OUT/native_rna_${STAGE}.metacell.h5ad"
    # Input checks and the floor assertion are skipped in resolve-only mode, which only answers
    # "does this script resolve this arm name?". preflight_arm.sh checks these inputs separately.
    if [ "${AMLGRN_RESOLVE_ONLY:-0}" != "1" ]; then
      [ -s "$GEX" ] || die "missing $GEX (built once per stage; k-independent, not rebuilt here)"

      # Check the floor: SCENIC+ draws each metacell with random.sample(cells_of_that_type, k)
      # without replacement, so a cell type with fewer than k cells makes prepare_GEX_ACC fail
      # ("Sample larger than population").
      python "$(dirname "$0")/metacell_floor.py" --gex "$GEX" --stage "$STAGE" --k "$KMETA" \
        || die "metacell floor check failed for k=$KMETA"
    fi
    ;;
  *) die "unknown pairing $PAIRING" ;;
esac

if [ "${AMLGRN_RESOLVE_ONLY:-0}" = "1" ]; then
  echo "resolve: arm='$PAIRING' -> GEX=$(basename "$GEX")${KMETA:+  metacell_k=$KMETA}"
  exit 0
fi

[ -s "$GEX" ] || die "missing GEX $GEX (run build_paired_rna.py / build_native_rna.py)"

MENR="$OUT/menr_cistarget_${STAGE}.hdf5"
[ -s "$MENR" ] || die "missing/empty $MENR -- run_full.sh step [5b] builds it; run_grn.sh no longer does"
bash "$(dirname "$0")/chain_guard.sh" verify "$STAGE" "arm:$PAIRING" || die "chain_guard verify"

# --- B) prepare_GEX_ACC (pairing-specific) -----------------------------------------------
MUD="$OUT/mudata_${TAG}.h5mu"
if [ ! -s "$MUD" ]; then
  step "prepare_GEX_ACC ($PAIRING)"
  if [ -z "$KMETA" ]; then   # [PBMC] truth, link, link_1to1 and every rand_* arm are multiome-mode
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$OUT/cistopic_obj_${STAGE}.pkl" \
      --GEX_anndata_fname "$GEX" --out_file "$MUD" \
      --do_not_use_raw_for_GEX_anndata || die "prepare_GEX_ACC(link)"
  else
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$OUT/cistopic_obj_${STAGE}.pkl" \
      --GEX_anndata_fname "$GEX" --out_file "$MUD" \
      --is_not_multiome --key_to_group_by Classified_Celltype \
      --nr_cells_per_metacells "$KMETA" --do_not_use_raw_for_GEX_anndata \
      || die "prepare_GEX_ACC($PAIRING, k=$KMETA)"
  fi
fi

step "search_spance"
[ -s "$OUT/search_space_${TAG}.tsv" ] || scenicplus prepare_data search_spance \
  --multiome_mudata_fname "$MUD" \
  --gene_annotation_fname "$DB/hg38_genome_annotation.tsv" \
  --chromsizes_fname "$DB/hg38_chromsizes_3col.tsv" \
  --out_fname "$OUT/search_space_${TAG}.tsv" || die "search_spance"

step "prepare_menr (reuse $MENR, re-filter TFs to this GEX)"
[ -s "$OUT/tf_names_${TAG}.txt" ] || scenicplus prepare_data prepare_menr \
  --paths_to_motif_enrichment_results "$MENR" \
  --multiome_mudata_fname "$MUD" \
  --out_file_tf_names "$OUT/tf_names_${TAG}.txt" \
  --out_file_direct_annotation "$OUT/cistromes_direct_${TAG}.h5ad" \
  --out_file_extended_annotation "$OUT/cistromes_extended_${TAG}.h5ad" \
  --direct_annotation Direct_annot \
  --extended_annotation Orthology_annot Motif_similarity_annot || die "prepare_menr"

step "region_to_gene"
[ -s "$OUT/r2g_${TAG}.tsv" ] || scenicplus grn_inference region_to_gene \
  --multiome_mudata_fname "$MUD" --search_space_fname "$OUT/search_space_${TAG}.tsv" \
  --temp_dir "$TMP" --out_region_to_gene_adjacencies "$OUT/r2g_${TAG}.tsv" \
  --importance_scoring_method GBM --correlation_scoring_method SR \
  --n_cpu "${AMLGRN_NCPU:-8}" || die "region_to_gene"

step "TF_to_gene  (the gate — heaviest for Option A on individual cells)"
[ -s "$OUT/tf2g_${TAG}.tsv" ] || scenicplus grn_inference TF_to_gene \
  --multiome_mudata_fname "$MUD" --tf_names "$OUT/tf_names_${TAG}.txt" \
  --temp_dir "$TMP" --out_tf_to_gene_adjacencies "$OUT/tf2g_${TAG}.tsv" \
  --method GBM --n_cpu "${AMLGRN_NCPU:-8}" || die "TF_to_gene"

step "eGRN"
if [ ! -s "$OUT/eRegulons_${TAG}.tsv" ]; then
  scenicplus grn_inference eGRN \
    --TF_to_gene_adj_fname "$OUT/tf2g_${TAG}.tsv" \
    --region_to_gene_adj_fname "$OUT/r2g_${TAG}.tsv" \
    --cistromes_fname "$OUT/cistromes_direct_${TAG}.h5ad" \
    --ranking_db_fname "$RANK" \
    --eRegulon_out_fname "$OUT/eRegulons_${TAG}.tsv" \
    --temp_dir "$TMP"
  if [ $? -ne 0 ]; then
    echo "--- eGRN returned non-zero; asking egrn_guard.py whether this is a result or a defect ---"
    python "$(dirname "$0")/egrn_guard.py" \
      --r2g "$OUT/r2g_${TAG}.tsv" --tag "$TAG" \
      --out "$OUT/eGRN_IMPOSSIBLE_${TAG}.json"
    if [ $? -eq 0 ]; then
      echo "=== [$TAG] eGRN STRUCTURALLY IMPOSSIBLE (confirmed from r2g) -- arm terminates here ==="
      echo "=== [$TAG] L1/L2 remain scorable from r2g/tf2g; L3-L5 do not exist for this arm ==="
      exit 45
    fi
    die "eGRN"
  fi
fi

step "AUCell"
[ -s "$OUT/eRegulon_AUC_${TAG}.h5mu" ] || scenicplus grn_inference AUCell \
  --eRegulon_fname "$OUT/eRegulons_${TAG}.tsv" \
  --multiome_mudata_fname "$MUD" \
  --aucell_out_fname "$OUT/eRegulon_AUC_${TAG}.h5mu" --n_cpu 1 || die "AUCell"

echo "=== [$TAG] GRN COMPLETE ==="; wc -l "$OUT/eRegulons_${TAG}.tsv" 2>/dev/null
