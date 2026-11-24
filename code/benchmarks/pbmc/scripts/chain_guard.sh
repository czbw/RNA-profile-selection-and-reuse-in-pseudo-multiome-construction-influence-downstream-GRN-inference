#!/usr/bin/env bash
set -euo pipefail

MODE="${1:?record|verify}"; STAGE="${2:?stage}"; WHO="${3:-unspecified}"
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT}"; OUT="$ROOT/out"

# ${VAR:?} catches only unset or empty values. A common bad value is a literal '$PWD' (from pasting
# an ssh-quoted command into an interactive shell), so validate that ROOT is a directory.
if [ ! -d "$ROOT" ]; then
  echo "chain_guard: AMLGRN_ROOT=$ROOT is not a directory." >&2
  case "$ROOT" in
    *'$'*) echo "  It still contains a literal '\$' -- you likely pasted an ssh-quoted command into"
           echo "  an interactive shell. Run:  export AMLGRN_ROOT=\$(pwd)   (no backslash)." >&2 ;;
  esac
  exit 7
fi
[ -d "$OUT" ] || { echo "chain_guard: $OUT does not exist -- wrong AMLGRN_ROOT?" >&2; exit 7; }
REF="$OUT/atac_chain_${STAGE}.sha256"
PIN="$(dirname "$0")/atac_chain_${STAGE}.sha256.pin"   # optional cross-run pin, version-controlled

CIS="$OUT/cistopic_obj_${STAGE}.pkl"
RGS="$OUT/region_sets_${STAGE}"
MEN="$OUT/menr_cistarget_${STAGE}.hdf5"

digest () {
  # emit three "sha256  label" lines, deterministic order
  [ -s "$CIS" ] || { echo "chain_guard: MISSING/EMPTY $CIS" >&2; exit 3; }
  [ -d "$RGS" ] || { echo "chain_guard: MISSING $RGS" >&2; exit 3; }
  [ -s "$MEN" ] || { echo "chain_guard: MISSING/EMPTY $MEN" >&2; exit 3; }
  echo "$(sha256sum "$CIS" | cut -d' ' -f1)  cistopic_obj"
  # Location-independent tree digest: sha256sum prints the file name beside each hash, so absolute
  # and relative paths would give different digests for identical files. cd into the directory and
  # hash './name' paths so the digest depends on file contents and names only. LC_ALL=C fixes the
  # sort order across locales.
  echo "$( (cd "$RGS" && find . -type f | LC_ALL=C sort | xargs sha256sum | sha256sum) \
           | cut -d' ' -f1 )  region_sets_tree"
  echo "$(sha256sum "$MEN" | cut -d' ' -f1)  menr_cistarget"
}

case "$MODE" in
  record)
    digest > "$REF"
    echo "[chain_guard] RECORDED $REF"
    sed 's/^/[chain_guard]   /' "$REF"
    if [ -s "$PIN" ]; then
      if diff -q "$REF" "$PIN" >/dev/null; then
        echo "[chain_guard] matches version-controlled pin: OK"
      else
        echo "[chain_guard] *** RECORDED CHAIN DIFFERS FROM THE PIN ***" >&2
        diff "$PIN" "$REF" >&2 || true
        echo "[chain_guard] the pin is the chain the published arms consumed. Aborting." >&2
        exit 4
      fi
    fi
    ;;
  verify)
    [ -s "$REF" ] || { echo "[chain_guard] no $REF -- run 'record' before the first arm" >&2; exit 5; }
    if digest | diff -q - "$REF" >/dev/null; then
      echo "[chain_guard] verify($WHO): chain unchanged -- OK"
    else
      echo "[chain_guard] *** CHAIN CHANGED between record and $WHO -- ABORTING ***" >&2
      digest | diff "$REF" - >&2 || true
      exit 6
    fi
    ;;
  *) echo "usage: chain_guard.sh record|verify <STAGE> [who]" >&2; exit 2 ;;
esac
