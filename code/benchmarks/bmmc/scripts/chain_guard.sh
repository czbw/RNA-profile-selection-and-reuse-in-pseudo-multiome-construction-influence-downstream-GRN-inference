#!/usr/bin/env bash
set -euo pipefail

MODE="${1:?record|verify}"; STAGE="${2:?stage}"; WHO="${3:-unspecified}"
ROOT="${AMLGRN_ROOT:?export AMLGRN_ROOT}"; OUT="$ROOT/out"

# ${VAR:?} only catches unset or empty values; it accepts a set-but-nonsense one, e.g. when an
# ssh-quoted command (export AMLGRN_ROOT=\$PWD) is pasted into an interactive shell and ROOT becomes
# the literal string '$PWD'. The failure would otherwise surface later as an opaque redirect
# error, so validate it here.
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
  local tree_hash
  tree_hash="$( cd "$RGS" && find . -type f -print0 | LC_ALL=C sort -z \
                | xargs -0 sha256sum | sha256sum | cut -d' ' -f1 )" \
    || { echo "chain_guard: hashing $RGS failed -- refusing a partial tree digest" >&2; exit 3; }
  echo "$tree_hash  region_sets_tree"
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
