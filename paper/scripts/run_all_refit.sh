#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python3 "$here/run.py" refit "$@"
python3 "$here/normalization/refit.py" --initialized "$@"
python3 "$here/calibration/refit.py" --initialized "$@"
python3 "$here/robustness/refit.py" all --initialized "$@"
