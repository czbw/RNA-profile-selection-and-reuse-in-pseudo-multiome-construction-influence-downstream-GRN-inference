#!/usr/bin/env bash
set -uo pipefail
STAGE="${1:?stage}"; PAIRING="${2:?link|pseudomultiome}"
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT}"
OUT="$ROOT/out"; DB="$ROOT/databases"; TMP="$ROOT/work/grn_tmp_${STAGE}_${PAIRING}"
mkdir -p "$TMP"
RANK="$DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
MOTIF_TBL="$DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TAG="${STAGE}_${PAIRING}"
step() { echo "===== [$TAG] $* ====="; }
die()  { echo "FAILED [$TAG]: $*"; exit 1; }

case "$PAIRING" in
  link)           GEX="$OUT/paired_rna_${STAGE}.h5ad" ;;     # Option A (multiome; cell_id_sep '___' -> names match cisTopic)
  pseudomultiome)                                            # Option B (non-multiome metacells)
    # SCENIC+'s metacell generator fails on very rare cell types: nr_metacells = round(n/5)*2 is
    # 2 even for 3-4 cells, so random.sample() is asked for 5 cells from 3 ("Sample larger than
    # population"). Drop cell types too thin in either modality first.
    GEX="$OUT/native_rna_${STAGE}.metacell.h5ad"
    [ -f "$GEX" ] || python filter_gex_for_metacells.py "$STAGE" || die "filter_gex_for_metacells"
    ;;
  *) die "unknown pairing $PAIRING" ;;
esac
[ -f "$GEX" ] || die "missing GEX $GEX (run build_paired_rna.py / build_native_rna.py)"

# --- A) motif enrichment (per STAGE, shared by both pairings) -----------------------------
MENR="$OUT/menr_cistarget_${STAGE}.hdf5"
if [ ! -f "$MENR" ]; then
  step "motif_enrichment_cistarget (DB)"
  scenicplus grn_inference motif_enrichment_cistarget \
    --region_set_folder "$OUT/region_sets_${STAGE}" \
    --cistarget_db_fname "$RANK" \
    --output_fname_cistarget_result "$MENR" \
    --temp_dir "$TMP" --species homo_sapiens \
    --path_to_motif_annotations "$MOTIF_TBL" \
    --annotation_version v10nr_clust --n_cpu 1 || die "motif_enrichment_cistarget"
fi

# --- B) prepare_GEX_ACC (pairing-specific) -----------------------------------------------
MUD="$OUT/mudata_${TAG}.h5mu"
if [ ! -f "$MUD" ]; then
  step "prepare_GEX_ACC ($PAIRING)"
  if [ "$PAIRING" = "link" ]; then
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$OUT/cistopic_obj_${STAGE}.pkl" \
      --GEX_anndata_fname "$GEX" --out_file "$MUD" \
      --do_not_use_raw_for_GEX_anndata || die "prepare_GEX_ACC(link)"
  else
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$OUT/cistopic_obj_${STAGE}.pkl" \
      --GEX_anndata_fname "$GEX" --out_file "$MUD" \
      --is_not_multiome --key_to_group_by Classified_Celltype \
      --nr_cells_per_metacells 5 --do_not_use_raw_for_GEX_anndata || die "prepare_GEX_ACC(pseudomultiome)"
  fi
fi

step "search_spance"
[ -f "$OUT/search_space_${TAG}.tsv" ] || scenicplus prepare_data search_spance \
  --multiome_mudata_fname "$MUD" \
  --gene_annotation_fname "$DB/hg38_genome_annotation.tsv" \
  --chromsizes_fname "$DB/hg38_chromsizes_3col.tsv" \
  --out_fname "$OUT/search_space_${TAG}.tsv" || die "search_spance"

step "prepare_menr (reuse $MENR, re-filter TFs to this GEX)"
[ -f "$OUT/tf_names_${TAG}.txt" ] || scenicplus prepare_data prepare_menr \
  --paths_to_motif_enrichment_results "$MENR" \
  --multiome_mudata_fname "$MUD" \
  --out_file_tf_names "$OUT/tf_names_${TAG}.txt" \
  --out_file_direct_annotation "$OUT/cistromes_direct_${TAG}.h5ad" \
  --out_file_extended_annotation "$OUT/cistromes_extended_${TAG}.h5ad" \
  --direct_annotation Direct_annot \
  --extended_annotation Orthology_annot Motif_similarity_annot || die "prepare_menr"

step "region_to_gene"
[ -f "$OUT/r2g_${TAG}.tsv" ] || scenicplus grn_inference region_to_gene \
  --multiome_mudata_fname "$MUD" --search_space_fname "$OUT/search_space_${TAG}.tsv" \
  --temp_dir "$TMP" --out_region_to_gene_adjacencies "$OUT/r2g_${TAG}.tsv" \
  --importance_scoring_method GBM --correlation_scoring_method SR \
  --n_cpu "${AMLGRN_NCPU:-8}" || die "region_to_gene"

step "TF_to_gene  (heaviest step for Option A on individual cells)"
[ -f "$OUT/tf2g_${TAG}.tsv" ] || scenicplus grn_inference TF_to_gene \
  --multiome_mudata_fname "$MUD" --tf_names "$OUT/tf_names_${TAG}.txt" \
  --temp_dir "$TMP" --out_tf_to_gene_adjacencies "$OUT/tf2g_${TAG}.tsv" \
  --method GBM --n_cpu "${AMLGRN_NCPU:-8}" || die "TF_to_gene"

step "eGRN"
[ -f "$OUT/eRegulons_${TAG}.tsv" ] || scenicplus grn_inference eGRN \
  --TF_to_gene_adj_fname "$OUT/tf2g_${TAG}.tsv" \
  --region_to_gene_adj_fname "$OUT/r2g_${TAG}.tsv" \
  --cistromes_fname "$OUT/cistromes_direct_${TAG}.h5ad" \
  --ranking_db_fname "$RANK" \
  --eRegulon_out_fname "$OUT/eRegulons_${TAG}.tsv" \
  --temp_dir "$TMP" || die "eGRN"

step "AUCell"
[ -f "$OUT/eRegulon_AUC_${TAG}.h5mu" ] || scenicplus grn_inference AUCell \
  --eRegulon_fname "$OUT/eRegulons_${TAG}.tsv" \
  --multiome_mudata_fname "$MUD" \
  --aucell_out_fname "$OUT/eRegulon_AUC_${TAG}.h5mu" --n_cpu 1 || die "AUCell"

echo "=== [$TAG] GRN COMPLETE ==="; wc -l "$OUT/eRegulons_${TAG}.tsv" 2>/dev/null
