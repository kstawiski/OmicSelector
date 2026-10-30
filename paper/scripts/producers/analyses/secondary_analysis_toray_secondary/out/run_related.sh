#!/usr/bin/env bash
set -euo pipefail
cd './analyses/secondary_analysis_toray_secondary/out'
bash run_r.sh scripts/prepare_fixture.R > logs/prepare.log 2>&1
bash run_r.sh scripts/run_analysis.R > logs/related.log 2>&1
bash run_r.sh scripts/retain_full_precision.R > logs/precision.log 2>&1
