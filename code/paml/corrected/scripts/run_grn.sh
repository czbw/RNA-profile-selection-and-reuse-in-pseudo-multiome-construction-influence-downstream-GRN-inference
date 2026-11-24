#!/usr/bin/env bash
set -uo pipefail
STAGE="${1:?stage}"; PAIRING="${2:?link|pseudomultiome|repair}"
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT}"
OUT="$ROOT/out"; DB="$ROOT/databases"; TMP="$ROOT/work/grn_tmp_${STAGE}_${PAIRING}"
mkdir -p "$TMP"
RANK="$DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
MOTIF_TBL="$DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TAG="${STAGE}_${PAIRING}"
step() { echo "===== [$TAG] $* ====="; }
die()  { echo "FAILED [$TAG]: $*"; exit 1; }

# path/to/x.h5mu -> path/to/x.part.h5mu   (extension preserved: anndata/mudata/polars all
# dispatch on it, so the sidecar must not be a bare ".part").
partof() { local f="$1"; printf '%s.part.%s' "${f%.*}" "${f##*.}"; }
# Move every .part sidecar into place; call only after the writing command returned 0.
commit() { local f; for f in "$@"; do mv -f "$(partof "$f")" "$f" || return 1; done; }

case "$PAIRING" in
  link)           GEX="$OUT/paired_rna_${STAGE}.h5ad" ;;
  repair)         GEX="$OUT/paired_rna_${STAGE}.repair.h5ad" ;;
                  # within-stratum pairing (p_a01). TAG becomes <STAGE>_repair, so outputs do not
                  # collide with <STAGE>_link. Reuses the same cistopic_obj / region_sets / menr;
                  # the ATAC-only topic model is not re-run.
  pseudomultiome)                                            # Option B (non-multiome metacells)
    # SCENIC+'s metacell generator fails on very rare cell types (it samples 5 cells from a type
    # with only 3-4). Drop cell types too thin in either modality first.
    GEX="$OUT/native_rna_${STAGE}.metacell.h5ad"
    [ -s "$GEX" ] || python filter_gex_for_metacells.py "$STAGE" || die "filter_gex_for_metacells"
    ;;
  *) die "unknown pairing $PAIRING" ;;
esac
[ -s "$GEX" ] || die "missing GEX $GEX (run build_paired_rna.py / build_native_rna.py)"

# --- A) motif enrichment (per STAGE, shared by both pairings) -----------------------------
MENR="$OUT/menr_cistarget_${STAGE}.hdf5"
if [ ! -s "$MENR" ]; then
  step "motif_enrichment_cistarget (DB)"
  scenicplus grn_inference motif_enrichment_cistarget \
    --region_set_folder "$OUT/region_sets_${STAGE}" \
    --cistarget_db_fname "$RANK" \
    --output_fname_cistarget_result "$(partof "$MENR")" \
    --temp_dir "$TMP" --species homo_sapiens \
    --path_to_motif_annotations "$MOTIF_TBL" \
    --annotation_version v10nr_clust --n_cpu 1 || die "motif_enrichment_cistarget"
  commit "$MENR" || die "motif_enrichment_cistarget (commit)"
fi

# --- B) prepare_GEX_ACC (pairing-specific) -----------------------------------------------
# Memory: this step and AUCell below both materialise a dense (cells x consensus regions) matrix.
# Diagnosis/malignant has 624,963 regions, so ~107 GB (float32) per copy at the pseudomultiome's
# 42,622 metacells and ~2.5x that at link's ~100k cells. AMLGRN_IMPUTE_BUDGET_GB does not limit
# these steps, so Diagnosis needs a large-memory node.
MUD="$OUT/mudata_${TAG}.h5mu"
if [ ! -s "$MUD" ]; then
  step "prepare_GEX_ACC ($PAIRING)"
  if [ "$PAIRING" = "link" ] || [ "$PAIRING" = "repair" ]; then
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$OUT/cistopic_obj_${STAGE}.pkl" \
      --GEX_anndata_fname "$GEX" --out_file "$(partof "$MUD")" \
      --do_not_use_raw_for_GEX_anndata || die "prepare_GEX_ACC(link)"
  else
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$OUT/cistopic_obj_${STAGE}.pkl" \
      --GEX_anndata_fname "$GEX" --out_file "$(partof "$MUD")" \
      --is_not_multiome --key_to_group_by Classified_Celltype \
      --nr_cells_per_metacells 5 --do_not_use_raw_for_GEX_anndata || die "prepare_GEX_ACC(pseudomultiome)"
  fi
  commit "$MUD" || die "prepare_GEX_ACC (commit)"
fi

step "search_spance"
SS="$OUT/search_space_${TAG}.tsv"
if [ ! -s "$SS" ]; then
  scenicplus prepare_data search_spance \
    --multiome_mudata_fname "$MUD" \
    --gene_annotation_fname "$DB/hg38_genome_annotation.tsv" \
    --chromsizes_fname "$DB/hg38_chromsizes_3col.tsv" \
    --out_fname "$(partof "$SS")" || die "search_spance"
  commit "$SS" || die "search_spance (commit)"
fi

# prepare_menr writes three outputs in order: tf_names, direct cistromes, extended cistromes.
# The guard tests the last one written, so a run that died part-way is not accepted.
step "prepare_menr (reuse $MENR, re-filter TFs to this GEX)"
TFN="$OUT/tf_names_${TAG}.txt"
CIS_D="$OUT/cistromes_direct_${TAG}.h5ad"
CIS_E="$OUT/cistromes_extended_${TAG}.h5ad"
if [ ! -s "$CIS_E" ]; then
  scenicplus prepare_data prepare_menr \
    --paths_to_motif_enrichment_results "$MENR" \
    --multiome_mudata_fname "$MUD" \
    --out_file_tf_names "$(partof "$TFN")" \
    --out_file_direct_annotation "$(partof "$CIS_D")" \
    --out_file_extended_annotation "$(partof "$CIS_E")" \
    --direct_annotation Direct_annot \
    --extended_annotation Orthology_annot Motif_similarity_annot || die "prepare_menr"
  commit "$TFN" "$CIS_D" "$CIS_E" || die "prepare_menr (commit)"
fi

step "region_to_gene"
R2G="$OUT/r2g_${TAG}.tsv"
if [ ! -s "$R2G" ]; then
  scenicplus grn_inference region_to_gene \
    --multiome_mudata_fname "$MUD" --search_space_fname "$SS" \
    --temp_dir "$TMP" --out_region_to_gene_adjacencies "$(partof "$R2G")" \
    --importance_scoring_method GBM --correlation_scoring_method SR \
    --n_cpu "${AMLGRN_NCPU:-8}" || die "region_to_gene"
  commit "$R2G" || die "region_to_gene (commit)"
fi

step "TF_to_gene  (the gate — heaviest for Option A on individual cells)"
TF2G="$OUT/tf2g_${TAG}.tsv"
if [ ! -s "$TF2G" ]; then
  scenicplus grn_inference TF_to_gene \
    --multiome_mudata_fname "$MUD" --tf_names "$TFN" \
    --temp_dir "$TMP" --out_tf_to_gene_adjacencies "$(partof "$TF2G")" \
    --method GBM --n_cpu "${AMLGRN_NCPU:-8}" || die "TF_to_gene"
  commit "$TF2G" || die "TF_to_gene (commit)"
fi

step "eGRN"
EREG="$OUT/eRegulons_${TAG}.tsv"
if [ ! -s "$EREG" ]; then
  scenicplus grn_inference eGRN \
    --TF_to_gene_adj_fname "$TF2G" \
    --region_to_gene_adj_fname "$R2G" \
    --cistromes_fname "$CIS_D" \
    --ranking_db_fname "$RANK" \
    --eRegulon_out_fname "$(partof "$EREG")" \
    --temp_dir "$TMP" || die "eGRN"
  commit "$EREG" || die "eGRN (commit)"
fi

# AUCell is the memory peak of the GRN half: it re-reads the dense (cells x 624,963 regions)
# accessibility matrix and builds a same-shaped rank matrix. --n_cpu is already 1, so the only
# lever is the memory request (run Diagnosis on a large-memory node).
step "AUCell"
AUC="$OUT/eRegulon_AUC_${TAG}.h5mu"
if [ ! -s "$AUC" ]; then
  scenicplus grn_inference AUCell \
    --eRegulon_fname "$EREG" \
    --multiome_mudata_fname "$MUD" \
    --aucell_out_fname "$(partof "$AUC")" --n_cpu 1 || die "AUCell"
  commit "$AUC" || die "AUCell (commit)"
fi

echo "=== [$TAG] GRN COMPLETE ==="; wc -l "$EREG" 2>/dev/null
