#!/usr/bin/env bash
# Download the datasets used by the experiments.
#
# STUB. Nothing is implemented, deliberately.
#
# Several of these datasets carry redistribution terms, and some require
# registration. An auto-downloader that routes around either would be wrong
# regardless of convenience. Fill in one dataset at a time, after reading its
# licence, and record the retrieval date in docs/DATA.md.
#
# Usage (once implemented):
#   ./scripts/download_data.sh            # everything
#   ./scripts/download_data.sh sachs      # one dataset

set -euo pipefail

DATA_DIR="${DATA_DIR:-data}"
TARGET="${1:-all}"

cat <<'MSG'
TODO: dataset downloads are not implemented.

See docs/DATA.md for per-dataset provenance, licences and manual instructions.
Place downloaded files under data/<dataset>/ ; the directory is gitignored.

Datasets: sachs dream realcause ihdp acic twins jobs lalonde perturb_seq
MSG

echo
echo "Requested: ${TARGET}"
echo "Target directory: ${DATA_DIR}"
exit 1
