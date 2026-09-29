#!/bin/bash
set -euo pipefail

# Run this any time for a running summary of AU-test failure rates across
# every family processed so far in the bulk run -- emails the result.
#
# Usage:
#   bash scripts/email_failure_report.sh [email_address]
#
# Defaults to the address in EMAIL_TO below if none is given on the command
# line.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$(dirname "$PROJECT_DIR")/data"

EMAIL_TO="${1:-chris.creevey@qub.ac.uk}"
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

python3 "$PROJECT_DIR/src/generate_failure_report.py" \
    --output-base "${DATA_DIR}/sim_outputs_bulk" \
    --manifest "${DATA_DIR}/manifest/manifest_full.json" \
    --out-tsv "${DATA_DIR}/reports/failure_report_${TIMESTAMP}.tsv" \
    --email "$EMAIL_TO"
