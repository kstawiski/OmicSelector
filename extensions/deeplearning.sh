#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate base
until Rscript deeplearning_settings.R; do
  echo "Restarting..."
  sleep 10
done