#!/bin/bash
until Rscript deeplearning_settings.R; do
  echo "Restarting..."
  sleep 10
done