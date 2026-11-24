#!/usr/bin/env bash
AMLGRN_PROJECT="${AMLGRN_PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
set -euo pipefail

WHICH="${1:-}"
case "$WHICH" in
  3k|10k|both) ;;
  *) echo "usage: $0 {3k|10k|both}" >&2; exit 2 ;;
esac

LOCAL_ROOT="$AMLGRN_PROJECT"
LOCAL_DATA="$LOCAL_ROOT/data/pbmc"
LOCAL_BUNDLE="$LOCAL_ROOT/benchmarks/pbmc"
REMOTE_ROOT="/path/to/AML-GRN"
REMOTE_BUNDLE="$REMOTE_ROOT/benchmarks/pbmc"
REMOTE_DATA="$REMOTE_BUNDLE/data"

say() { printf '\n== %s\n' "$*"; }

# 0. destinations must exist before scp (scp will not create them)
say "creating remote directories"
ssh cubic "mkdir -p '$REMOTE_DATA' '$REMOTE_BUNDLE/scripts' '$REMOTE_BUNDLE/work' '$REMOTE_BUNDLE/out' '$REMOTE_BUNDLE/logs'"

# helper: copy only if absent or a different size on the remote
put () {
  local src="$1" dstdir="$2" name
  name="$(basename "$src")"
  if [ ! -f "$src" ]; then echo "  MISSING LOCALLY: $src" >&2; return 1; fi
  local lsize rsize
  lsize=$(stat -c%s "$src")
  rsize=$(ssh cubic "stat -c%s '$dstdir/$name' 2>/dev/null || echo 0")
  if [ "$lsize" = "$rsize" ]; then
    printf '  skip (already there, %s bytes)  %s\n' "$lsize" "$name"
    return 0
  fi
  printf '  sending %s (%s bytes) ...\n' "$name" "$lsize"
  scp -q "$src" "cubic:$dstdir/"
  rsize=$(ssh cubic "stat -c%s '$dstdir/$name'")
  if [ "$lsize" != "$rsize" ]; then
    echo "  SIZE MISMATCH after transfer: local $lsize vs remote $rsize -- re-run" >&2
    return 1
  fi
  printf '  ok\n'
}

# 1. the bundle itself (small -- scripts, config, runner, README)
say "uploading bundle files"
for f in run_full.sh submit_full.sh preflight_arm.sh stage_big_files.sh; do
  [ -f "$LOCAL_BUNDLE/$f" ] && put "$LOCAL_BUNDLE/$f" "$REMOTE_BUNDLE" || true
done
for f in "$LOCAL_BUNDLE"/scripts/*; do
  [ -f "$f" ] && put "$f" "$REMOTE_BUNDLE/scripts" || true
done

# 2. the data payload
stage_set () {
  local s="$1"
  say "uploading ${s} payload"
  put "$LOCAL_DATA/pbmc_granulocyte_sorted_${s}_filtered_feature_bc_matrix.h5" "$REMOTE_DATA"
  put "$LOCAL_DATA/pbmc_granulocyte_sorted_${s}_atac_fragments.tsv.gz"          "$REMOTE_DATA"
  put "$LOCAL_DATA/pbmc_granulocyte_sorted_${s}_atac_fragments.tsv.gz.tbi"      "$REMOTE_DATA"
}
# Explicit `if` blocks are used because `[ A ] || [ B ] && cmd` parses as (A || B) && cmd, which
# returns non-zero when neither holds and would abort the script under `set -e`.
if [ "$WHICH" = "3k" ]  || [ "$WHICH" = "both" ]; then stage_set 3k;  fi
if [ "$WHICH" = "10k" ] || [ "$WHICH" = "both" ]; then stage_set 10k; fi

# 3. verify the shared assets this bundle reuses rather than ships
say "checking shared assets already on CUBIC (reused, not uploaded)"
ssh cubic "
  ok=1
  for p in '$REMOTE_ROOT/miniforge.sif' '$REMOTE_ROOT/envs/scenicplus' '$REMOTE_ROOT/envs/macs3env'; do
    if [ -e \"\$p\" ]; then echo \"  present : \$p\"; else echo \"  MISSING : \$p\"; ok=0; fi
  done
  echo '  cisTarget feathers (hg38, reused from Arm A):'
  ls -la '$REMOTE_ROOT/paml/healthy/databases/' 2>/dev/null | grep -i feather | awk '{printf \"    %s  %s\n\", \$5, \$9}' || echo '    NONE FOUND -- check the Arm A bundle'
  [ \"\$ok\" = 1 ] || echo '  ^ build the missing env(s) before submitting (see submit_container.sh header)'
"

say "chmod +x the runners"
ssh cubic "chmod +x '$REMOTE_BUNDLE'/*.sh '$REMOTE_BUNDLE'/scripts/*.sh 2>/dev/null || true"

cat <<EOF

== DONE. Next:
   ssh cubic
   cd $REMOTE_ROOT
   sbatch benchmarks/pbmc/submit_full.sh

EOF
