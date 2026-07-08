#!/usr/bin/env bash
set -euo pipefail

APP_DIR="/OmicSelector/shiny/omicselector"
SHINY_HOST="${SHINY_HOST:-0.0.0.0}"
SHINY_PORT="${SHINY_PORT:-3838}"

if [ ! -d "$APP_DIR" ]; then
  echo "Shiny app not found at $APP_DIR"
  exit 1
fi

mkdir -p "${OMICSELECTOR_ANALYSES_DIR:-/OmicSelector/analyses}"

R -e "shiny::runApp('${APP_DIR}', host = '${SHINY_HOST}', port = ${SHINY_PORT})" &

exec /init
