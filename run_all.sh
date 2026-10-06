#!/usr/bin/env bash
# Reproduce the analysis end to end. Run from the repository root.
set -euo pipefail
[ -f .env ] && set -a && . ./.env && set +a

python  scripts/prepare_cohort.py       # raw exports -> cohort
python  scripts/validate_panel.py       # quality framework -> validated cohort
Rscript scripts/estimate_intervals.R    # reference intervals
Rscript scripts/plot_intervals.R        # figures 1 and 2
Rscript scripts/plot_comparison.R       # figure 3
python  scripts/compare_published.py    # comparison with published intervals

echo "Done. Intervals: ${RESULT_DIR:-results}/reference_intervals.csv"
