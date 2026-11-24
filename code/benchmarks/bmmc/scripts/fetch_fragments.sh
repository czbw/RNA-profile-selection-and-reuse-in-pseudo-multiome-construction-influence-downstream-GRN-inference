#!/usr/bin/env bash
set -uo pipefail

# --- 13 batches, from obs/batch; identical to the 13 S3 sample prefixes ---------------------
BATCHES="s1d1 s1d2 s1d3 s2d1 s2d4 s2d5 s3d3 s3d6 s3d7 s3d10 s4d1 s4d8 s4d9"

S3_PREFIX="s3://openproblems-bio/public/post_competition/multiome"
# anonymous HTTPS endpoint for the same objects — the curl fallback for machines with no AWS CLI
HTTP_PREFIX="https://openproblems-bio.s3.amazonaws.com/public/post_competition/multiome"

DEST=""
WITH_METRICS=0
for a in "$@"; do
  case "$a" in
    --with-metrics) WITH_METRICS=1 ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) DEST="$a" ;;
  esac
done
DEST="${DEST:-$(cd "$(dirname "$0")/.." && pwd)/data/fragments}"

FILES="atac_fragments.tsv.gz atac_fragments.tsv.gz.tbi"
[ "$WITH_METRICS" -eq 1 ] && FILES="$FILES per_barcode_metrics.csv"

HAVE_AWS=0
command -v aws >/dev/null 2>&1 && HAVE_AWS=1

echo "=================================================================="
echo " BMMC fragments   $(date)"
echo " dest    : $DEST"
echo " batches : 13"
echo " files   : $FILES"
echo " method  : $([ "$HAVE_AWS" -eq 1 ] && echo 'aws s3 cp --no-sign-request' || echo 'curl (no AWS CLI found)')"
echo "=================================================================="

mkdir -p "$DEST" || { echo "cannot create $DEST"; exit 2; }

FAIL=0
N_OK=0
N_SKIP=0
TOTAL=0

human() {  # bytes -> human, no dependency on numfmt (absent on some git-bash installs)
  awk -v b="$1" 'BEGIN{
    split("B KiB MiB GiB TiB",u," "); i=1
    while (b>=1024 && i<5) { b/=1024; i++ }
    printf (i==1 ? "%d %s" : "%.2f %s"), b, u[i]
  }'
}

size_of() { wc -c < "$1" 2>/dev/null | tr -d ' ' || echo 0; }

# Assert the Cell Ranger ARC header. gzip -dc on a BGZF file inflates the first block fine.
# `head -c` closes the pipe, so gzip takes SIGPIPE and returns non-zero under pipefail --
# hence the capture-with-|| true rather than a bare `| grep -q`.
verify_fragments() {
  local f="$1" head_txt
  head_txt="$(gzip -dc "$f" 2>/dev/null | head -c 8192 || true)"
  case "$head_txt" in
    *"pipeline_name=cellranger-arc"*) return 0 ;;
  esac
  echo
  echo "  *** HEADER CHECK FAILED: $f"
  echo "  *** expected the inflated first block to contain 'pipeline_name=cellranger-arc'."
  echo "  *** got (first 200 bytes, may be an S3 error document or a truncated transfer):"
  printf '%s\n' "${head_txt:0:200}" | sed 's/^/  *** | /'
  echo "  *** DELETE this file and re-run. Do NOT feed it to ArchR."
  return 1
}

get_one() {
  local batch="$1" name="$2" dst="$3" rc=0
  if [ "$HAVE_AWS" -eq 1 ]; then
    aws s3 cp --no-sign-request "$S3_PREFIX/$batch/$name" "$dst" --only-show-errors
    rc=$?
    [ "$rc" -eq 0 ] && return 0
    echo "    aws s3 cp returned $rc -- falling back to curl"
  fi
  # -f: fail on 4xx/5xx instead of writing the error body to disk (that is exactly the
  # non-zero-size garbage the header check exists to catch). -L: follow redirects.
  curl -fL --retry 3 --retry-delay 5 -o "$dst" "$HTTP_PREFIX/$batch/$name"
}

for BATCH in $BATCHES; do
  mkdir -p "$DEST/$BATCH"
  echo
  echo "---- $BATCH ----"
  for NAME in $FILES; do
    DST="$DEST/$BATCH/$NAME"

    if [ -s "$DST" ]; then
      SZ="$(size_of "$DST")"
      # a present-but-corrupt .tsv.gz must not be accepted just because it is non-empty
      if [ "$NAME" = "atac_fragments.tsv.gz" ] && ! verify_fragments "$DST"; then
        FAIL=$((FAIL+1)); continue
      fi
      printf '  skip  %-26s %12s  (%s)\n' "$NAME" "$SZ" "$(human "$SZ")"
      TOTAL=$((TOTAL+SZ)); N_SKIP=$((N_SKIP+1))
      continue
    fi

    rm -f "$DST"                      # a 0-byte leftover from a killed transfer
    printf '  get   %-26s ...\n' "$NAME"
    if ! get_one "$BATCH" "$NAME" "$DST"; then
      echo "  DOWNLOAD FAILED: $BATCH/$NAME"
      rm -f "$DST"                    # never leave a partial file for the next run to accept
      FAIL=$((FAIL+1)); continue
    fi

    SZ="$(size_of "$DST")"
    if [ ! -s "$DST" ]; then
      echo "  EMPTY AFTER DOWNLOAD: $BATCH/$NAME"; rm -f "$DST"; FAIL=$((FAIL+1)); continue
    fi
    if [ "$NAME" = "atac_fragments.tsv.gz" ] && ! verify_fragments "$DST"; then
      FAIL=$((FAIL+1)); continue
    fi
    printf '  ok    %-26s %12s  (%s)\n' "$NAME" "$SZ" "$(human "$SZ")"
    TOTAL=$((TOTAL+SZ)); N_OK=$((N_OK+1))
  done
done

echo
echo "==== FRAGMENT FETCH SUMMARY  $(date) ===="
echo "  downloaded : $N_OK"
echo "  already ok : $N_SKIP"
echo "  total bytes: $TOTAL  ($(human "$TOTAL"))"
echo "  dest       : $DEST"
if [ "$FAIL" -ne 0 ]; then
  echo "  $FAIL file(s) FAILED -- see above. Fix and re-run; completed files are skipped."
  echo "==== INCOMPLETE ===="
  exit 1
fi
echo "  every atac_fragments.tsv.gz inflated and carries pipeline_name=cellranger-arc"
echo
echo "  NEXT: these are inputs to the ArchR link construction ONLY, which runs LOCALLY."
echo "  Do NOT scp them to CUBIC -- the peak-matrix chain and the Seurat link arm do not read"
echo "  them, and the per-cell assignment table they produce is a small CSV."
echo "==== FRAGMENT FETCH COMPLETE ===="
