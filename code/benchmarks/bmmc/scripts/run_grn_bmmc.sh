#!/usr/bin/env bash
# Per-arm GRN half. The frozen shared chromatin chain is re-verified at entry, so the only thing
# that differs between arms is the cell pairing.
set -uo pipefail

STAGE="${1:?stage}"; ARM="${2:?arm name, as dispatched by run_bmmc.sh}"
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT (run_bmmc.sh does)}"
OUT="$ROOT/out"
DB="${AMLGRN_DB:-$ROOT/databases}"          # the GRCh38 feathers live one level up on CUBIC
TMP="$ROOT/work/grn_tmp_${STAGE}_${ARM}"
mkdir -p "$TMP" || { echo "FATAL: cannot create $TMP"; exit 1; }

RANK="$DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
MOTIF_TBL="$DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TAG="${STAGE}_${ARM}"
NCPU="${AMLGRN_NCPU:-8}"

# non-empty only for metacell (non-multiome) arms; it also selects the prepare_GEX_ACC branch
KMETA=""
# The obs column SCENIC+ groups metacells by. Both reference bundles hardcode Classified_Celltype
# and config.yaml declares the same name under pseudobulk.group_by, so the arm builder must
# emit that column (it copies its own label column into Classified_Celltype, because the SCENIC+
# call passes --key_to_group_by Classified_Celltype).
GROUP_BY="${AMLGRN_GROUP_BY:-Classified_Celltype}"

step() { echo "===== [$TAG] $* ====="; }
die()  { echo "FAILED [$TAG]: $*" >&2; exit 1; }

# path/to/x.h5mu -> path/to/x.part.h5mu   (extension preserved; see the resumability note above)
partof() { local f="$1"; printf '%s.part.%s' "${f%.*}" "${f##*.}"; }
# Move every .part sidecar into place; call only after the writing command returned 0.
commit() { local f; for f in "$@"; do mv -f "$(partof "$f")" "$f" || return 1; done; }

echo "=================================================================="
echo " [$TAG] GRN half   $(date)"
echo " root : $ROOT"
echo " db   : $DB"
echo " tmp  : $TMP     ncpu: $NCPU   seed: ${AMLGRN_SEED:-unset}"
echo "=================================================================="

CIS="$OUT/cistopic_obj_${STAGE}.pkl"
MENR="$OUT/menr_cistarget_${STAGE}.hdf5"
[ -s "$CIS" ]  || die "missing/empty $CIS -- run_bmmc.sh phase [2a] builds it; this script never does"
[ -s "$MENR" ] || die "missing/empty $MENR -- run_bmmc.sh phase [2c] builds it; this script never does"
[ -d "$OUT/region_sets_${STAGE}" ] \
  || die "missing $OUT/region_sets_${STAGE} -- run_bmmc.sh phase [2b] builds it"
step "chain_guard verify (arm entry)"
bash chain_guard.sh verify "$STAGE" "arm-entry:$ARM" \
  || die "chain_guard verify -- the frozen chromatin chain has MOVED; every cross-arm number computed against it would be invalid"

# --- resolve this arm's GEX ---------------------------------------------------------------
# run_bmmc.sh dispatches an arm name and owns arm construction; phase [3] has already built the
# object this resolves to, and phase [1] already asserted every input the name implies. So this
# script keeps no arm whitelist: an unknown name has no RNA object and dies here on a named
# missing file, having written nothing.
case "$ARM" in
  pseudomultiome|pseudomultiome_k*)
    # Metacell (non-multiome) arms. SCENIC+ builds the metacells itself inside prepare_GEX_ACC, so
    # there is no per-cell assignment table; the metacell size is read from the arm name so two
    # metacell arms cannot silently share a setting and k lands in the output TAG.
    case "$ARM" in
      pseudomultiome)    KMETA=5 ;;                       # the field's common value
      pseudomultiome_k*) KMETA="${ARM#pseudomultiome_k}" ;;
    esac
    case "$KMETA" in ''|*[!0-9]*) die "cannot parse metacell size from arm '$ARM'" ;; esac
    [ "$KMETA" -ge 2 ] || die "metacell size must be >= 2 (got $KMETA)"
    # The metacell GEX is built once per stage and is independent of k: it is a cell-type filter
    # whose floor is max(min_cells_per_celltype=25, k), i.e. 25 for every k <= 25. It is shared by
    # both metacell arms and never rebuilt here. native_rna.py + filter_gex_for_metacells.py
    # build it (locally and copied into out/, or here); a missing file is a hard error rather than
    # a silent rebuild, because rebuilding inside a per-arm script would make a stage-scoped shared
    # artifact a side effect of one arm.
    GEX="$OUT/native_rna_${STAGE}.metacell.h5ad"
    ;;
  *)
    GEX="$(printf '%s/paired_rna_%s.%s.h5ad' "$OUT" "$STAGE" "$ARM")"
    ;;
esac
if [ -n "$KMETA" ]; then
  [ -s "$GEX" ] || die "missing GEX $GEX -- the metacell arms consume a NATIVE (unpaired) RNA object, which native_rna.py builds and filter_gex_for_metacells.py filters; build_arm_rna.py refuses metacell arms by design. See scripts/UPLOAD_MANIFEST.md sec.3"
else
  [ -s "$GEX" ] || die "missing GEX $GEX -- run_bmmc.sh phase [3/4] builds it (build_arm_rna.py)"
fi
echo "  arm '$ARM' -> GEX=$(basename "$GEX")${KMETA:+   metacell_k=$KMETA}"

if [ -n "$KMETA" ]; then
  step "metacell floor + cell-id assertion (b_12, k=$KMETA)"
  python metacell_floor.py --gex "$GEX" --stage "$STAGE" --k "$KMETA" \
    --analysis-set "$OUT/cistopic_cell_names_${STAGE}.txt" \
    --metadata "$ROOT/uploads/b00_cell_metadata.csv" --require-analysis-set \
    || die "metacell_floor.py failed for k=$KMETA (3 = a cell type is below the k floor and prepare_GEX_ACC would crash; 4 = the GEX cell ids are not the frozen analysis set; 2 = an input is missing)"
fi

# --- B) prepare_GEX_ACC (pairing-specific) ------------------------------------------------
# Memory: this step and AUCell below both materialise a dense (cells x regions) matrix. At BMMC's
# 69,249 cells x 116,490 peaks that is ~32 GB per copy in float32, comfortably inside the 1200 G
# submit_bmmc.sh requests. Neither step honours AMLGRN_IMPUTE_BUDGET_GB; that guard lives in the
# region-sets step and only protects the DAR branch.
MUD="$OUT/mudata_${TAG}.h5mu"
if [ ! -s "$MUD" ]; then
  step "prepare_GEX_ACC ($ARM)"
  if [ -z "$KMETA" ]; then   # truth and every rand_* arm are multiome-mode: they differ only in
                             # which RNA cell sits under each ATAC cell
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$CIS" \
      --GEX_anndata_fname "$GEX" --out_file "$(partof "$MUD")" \
      --do_not_use_raw_for_GEX_anndata || { rm -f "$(partof "$MUD")"; die "prepare_GEX_ACC($ARM)"; }
  else
    scenicplus prepare_data prepare_GEX_ACC \
      --cisTopic_obj_fname "$CIS" \
      --GEX_anndata_fname "$GEX" --out_file "$(partof "$MUD")" \
      --is_not_multiome --key_to_group_by "$GROUP_BY" \
      --nr_cells_per_metacells "$KMETA" --do_not_use_raw_for_GEX_anndata \
      || { rm -f "$(partof "$MUD")"; die "prepare_GEX_ACC($ARM, k=$KMETA)"; }
  fi
  commit "$MUD" || die "prepare_GEX_ACC (commit)"
fi

# --- C) search space ----------------------------------------------------------------------
# `search_spance` is not a typo: it is the subcommand name as shipped by the scenicplus 1.0a2 CLI.
# Do not "fix" it.
SS="$OUT/search_space_${TAG}.tsv"
if [ ! -s "$SS" ]; then
  step "search_spance"
  scenicplus prepare_data search_spance \
    --multiome_mudata_fname "$MUD" \
    --gene_annotation_fname "$DB/hg38_genome_annotation.tsv" \
    --chromsizes_fname "$DB/hg38_chromsizes_3col.tsv" \
    --out_fname "$(partof "$SS")" || { rm -f "$(partof "$SS")"; die "search_spance"; }
  commit "$SS" || die "search_spance (commit)"
fi

# --- D) prepare_menr: reuse the shared motif enrichment, re-filter TFs to this arm's GEX ---
# prepare_menr writes three outputs, in this order: tf_names -> direct cistromes -> extended
# cistromes. The guard must test the last one written; testing tf_names would accept a run that
# died between the first and third file.
TFN="$OUT/tf_names_${TAG}.txt"
CIS_D="$OUT/cistromes_direct_${TAG}.h5ad"
CIS_E="$OUT/cistromes_extended_${TAG}.h5ad"
if [ ! -s "$CIS_E" ]; then
  step "prepare_menr (reuse $MENR, re-filter TFs to this GEX)"
  scenicplus prepare_data prepare_menr \
    --paths_to_motif_enrichment_results "$MENR" \
    --multiome_mudata_fname "$MUD" \
    --out_file_tf_names "$(partof "$TFN")" \
    --out_file_direct_annotation "$(partof "$CIS_D")" \
    --out_file_extended_annotation "$(partof "$CIS_E")" \
    --direct_annotation Direct_annot \
    --extended_annotation Orthology_annot Motif_similarity_annot \
    || { rm -f "$(partof "$TFN")" "$(partof "$CIS_D")" "$(partof "$CIS_E")"; die "prepare_menr"; }
  commit "$TFN" "$CIS_D" "$CIS_E" || die "prepare_menr (commit)"
fi

R2G="$OUT/r2g_${TAG}.tsv"
if [ ! -s "$R2G" ]; then
  step "region_to_gene"
  scenicplus grn_inference region_to_gene \
    --multiome_mudata_fname "$MUD" --search_space_fname "$SS" \
    --temp_dir "$TMP" --out_region_to_gene_adjacencies "$(partof "$R2G")" \
    --importance_scoring_method GBM --correlation_scoring_method SR \
    --n_cpu "$NCPU" || { rm -f "$(partof "$R2G")"; die "region_to_gene"; }
  commit "$R2G" || die "region_to_gene (commit)"
fi

# --- F) TF_to_gene ------------------------------------------------------------------------
TF2G="$OUT/tf2g_${TAG}.tsv"
if [ ! -s "$TF2G" ]; then
  step "TF_to_gene  (the gate — heaviest step for a per-cell multiome arm)"
  scenicplus grn_inference TF_to_gene \
    --multiome_mudata_fname "$MUD" --tf_names "$TFN" \
    --temp_dir "$TMP" --out_tf_to_gene_adjacencies "$(partof "$TF2G")" \
    --method GBM --n_cpu "$NCPU" || { rm -f "$(partof "$TF2G")"; die "TF_to_gene"; }
  commit "$TF2G" || die "TF_to_gene (commit)"
fi

# --- G) eGRN, and the result-vs-defect decision (see the exit-45 header block) ---------------
EREG="$OUT/eRegulons_${TAG}.tsv"
if [ ! -s "$EREG" ]; then
  step "eGRN"
  scenicplus grn_inference eGRN \
    --TF_to_gene_adj_fname "$TF2G" \
    --region_to_gene_adj_fname "$R2G" \
    --cistromes_fname "$CIS_D" \
    --ranking_db_fname "$RANK" \
    --eRegulon_out_fname "$(partof "$EREG")" \
    --temp_dir "$TMP"
  rc=$?
  # Two ways to end up with no network: a crash (rc != 0), or a zero exit that wrote nothing. Both
  # mean "no eGRN table", so both go to the guard, which refuses anything that is not the specific
  # empty-adjacency signature.
  if [ "$rc" -ne 0 ] || [ ! -s "$(partof "$EREG")" ]; then
    rm -f "$(partof "$EREG")"
    echo "--- [$TAG] eGRN produced no eRegulon table (rc=$rc); asking egrn_guard.py whether that"
    echo "--- is a RESULT (this pairing supports no network) or a DEFECT (anything else) ---"
    python egrn_guard.py --r2g "$R2G" --tag "$TAG" \
      --out "$OUT/eGRN_IMPOSSIBLE_${TAG}.json"
    grc=$?
    if [ "$grc" -eq 0 ]; then
      echo "=== [$TAG] eGRN STRUCTURALLY IMPOSSIBLE (confirmed from r2g) -- arm terminates here ==="
      echo "=== [$TAG] L1/L2 remain scorable from r2g/tf2g; L3-L5 do not exist for this arm ==="
      echo "=== [$TAG] exiting 45: a recorded RESULT, not a failure (run_bmmc.sh phase [4/4]) ==="
      exit 45
    fi
    die "eGRN (rc=$rc; egrn_guard.py did not confirm the empty-adjacency signature, so this failure is a DEFECT and stays fatal)"
  fi
  commit "$EREG" || die "eGRN (commit)"
fi

# --- H) AUCell ----------------------------------------------------------------------------
# AUCell, not eGRN, is the memory peak of the GRN half: it re-reads the dense (cells x regions)
# accessibility matrix and builds a same-shaped rank matrix beside it. --n_cpu is already 1, so the
# only lever is the memory request, hence the bigmem partition in submit_bmmc.sh.
AUC="$OUT/eRegulon_AUC_${TAG}.h5mu"
if [ ! -s "$AUC" ]; then
  step "AUCell"
  scenicplus grn_inference AUCell \
    --eRegulon_fname "$EREG" \
    --multiome_mudata_fname "$MUD" \
    --aucell_out_fname "$(partof "$AUC")" --n_cpu 1 \
    || { rm -f "$(partof "$AUC")"; die "AUCell"; }
  commit "$AUC" || die "AUCell (commit)"
fi

# --- I) the chain must not have moved during this arm either ---------------------------------
# Entry verification shows what the arm started from; this shows nothing shifted underneath it.
# run_bmmc.sh verifies once more after the last arm; this is the per-arm half of that check.
bash chain_guard.sh verify "$STAGE" "arm-exit:$ARM" \
  || die "chain_guard verify (arm exit) -- the chain moved DURING this arm; its outputs are not comparable"

echo "=== [$TAG] GRN COMPLETE ==="
wc -l "$EREG" 2>/dev/null
