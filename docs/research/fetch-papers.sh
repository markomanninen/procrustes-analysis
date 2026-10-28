#!/bin/sh
# fetch-papers.sh: download the open-access PDFs listed in papers-manifest.tsv into papers/.
#
# Usage:  ./fetch-papers.sh                    (every paper in papers-manifest.tsv)
#         ./fetch-papers.sh papers-core.tsv    (only the core set)
#         ARXIV_DELAY=5 ./fetch-papers.sh      (seconds between arXiv requests; default 3)
#
# - Paths are resolved from this script's own folder, so it can be run from anywhere.
# - Files that already exist in papers/ are skipped, so a second run only retries what failed.
# - Each download goes to a .part file and is kept only if it starts with %PDF. Some hosts
#   answer curl with an HTML page (publisher bot checks, PMC's browser check); those are
#   reported as failures and can be fetched by browser from the URL printed in the summary.
# - papers/ is git-ignored: not every paper's licence allows redistribution.

set -u

SCRIPT_DIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd -P) || exit 1
MANIFEST=${1:-"$SCRIPT_DIR/papers-manifest.tsv"}
case $MANIFEST in /*) ;; *) [ -f "$MANIFEST" ] || MANIFEST="$SCRIPT_DIR/$MANIFEST" ;; esac
DEST="$SCRIPT_DIR/papers"
ARXIV_DELAY=${ARXIV_DELAY:-3}
OTHER_DELAY=${OTHER_DELAY:-1}
UA="procrustes-analysis research archive (fetch-papers.sh; curl)"

command -v curl >/dev/null 2>&1 || { echo "fetch-papers.sh: curl not found" >&2; exit 1; }
[ -f "$MANIFEST" ] || { echo "fetch-papers.sh: manifest not found: $MANIFEST" >&2; exit 1; }
mkdir -p "$DEST" || exit 1

TAB=$(printf '\t')
NL='
'
total=$(grep -v -e '^#' -e '^key' "$MANIFEST" | grep -c .)
downloaded=0
skipped=0
failed=0
arxiv_done=0
failures=""
tmp=""

trap 'if [ -n "$tmp" ]; then rm -f "$tmp"; fi; echo; echo "interrupted" >&2; exit 130' INT TERM

echo "Fetching $total files listed in $(basename "$MANIFEST") into $DEST"

while IFS="$TAB" read -r key file url host _ _ status _; do
  case $key in
    '' | '#'* | key) continue ;;
  esac
  case $file in
    */* | .* | '') echo "skip  unsafe filename: $file" >&2; continue ;;
  esac
  out="$DEST/$file"
  if [ -s "$out" ]; then
    skipped=$((skipped + 1))
    continue
  fi

  # Be polite: arXiv asks automated clients to space their requests.
  if [ "$host" = "arxiv.org" ]; then
    if [ "$arxiv_done" -gt 0 ]; then sleep "$ARXIV_DELAY"; fi
    arxiv_done=$((arxiv_done + 1))
  else
    sleep "$OTHER_DELAY"
  fi

  tmp="$out.part"
  if curl -fsSL --retry 2 --retry-delay 5 --connect-timeout 20 --max-time 300 \
       -A "$UA" -o "$tmp" "$url" </dev/null; then
    magic=$(dd if="$tmp" bs=4 count=1 2>/dev/null)
    if [ "$magic" = "%PDF" ]; then
      mv "$tmp" "$out"
      downloaded=$((downloaded + 1))
      echo "ok    $file"
    else
      rm -f "$tmp"
      failed=$((failed + 1))
      failures="$failures$NL  $file: not a PDF (the server sent a web page); $url"
      echo "fail  $file (not a PDF)"
    fi
  else
    rc=$?
    rm -f "$tmp"
    failed=$((failed + 1))
    failures="$failures$NL  $file: curl exit $rc (HEAD status in manifest: $status); $url"
    echo "fail  $file (curl exit $rc)"
  fi
  tmp=""
done < "$MANIFEST"

kb=$(du -sk "$DEST" | cut -f1)
echo
echo "Summary"
echo "  listed:      $total"
echo "  downloaded:  $downloaded"
echo "  skipped:     $skipped (already in papers/)"
echo "  failed:      $failed"
echo "  papers/ now holds $((kb / 1024)) MB"
if [ "$failed" -gt 0 ]; then
  echo
  echo "Failed (try these in a browser):"
  printf '%s\n' "$failures" | sed '/^$/d'
fi
[ "$failed" -eq 0 ]
