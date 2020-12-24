#!/bin/bash
until Rscript deeplearning_settings.R; do
  echo "Failed.. Restarting..."
  sleep 10
done