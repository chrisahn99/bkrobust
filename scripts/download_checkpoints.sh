#!/usr/bin/env bash
# Download the causal foundation model checkpoints for the audit.
#
# STUB. Nothing is implemented, deliberately.
#
# Checkpoint licences vary and must be read before download. Every checkpoint
# must also be pinned to a revision -- these models get updated, and an
# unpinned checkpoint makes every audit number unreproducible. Record source,
# revision and licence in docs/CHECKPOINTS.md as each is added.
#
# Usage (once implemented):
#   ./scripts/download_checkpoints.sh            # everything
#   ./scripts/download_checkpoints.sh causalpfn  # one checkpoint

set -euo pipefail

CHECKPOINT_DIR="${CHECKPOINT_DIR:-checkpoints}"
TARGET="${1:-all}"

cat <<'MSG'
TODO: checkpoint downloads are not implemented.

See docs/CHECKPOINTS.md for sources, revisions and licences.
Place downloaded checkpoints under checkpoints/<key>/ ; the directory is
gitignored.

Requires the [cfm] extra:  pip install -e ".[cfm]"

Checkpoints: causalpfn causalfm dopfn
MSG

echo
echo "Requested: ${TARGET}"
echo "Target directory: ${CHECKPOINT_DIR}"
exit 1
