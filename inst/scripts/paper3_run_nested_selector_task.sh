#!/usr/bin/env bash

# Run one task from paper3_build_nested_selector_tasks.R. Required environment:
# TASK_FILE, TASK_ID (or SGE_TASK_ID), CACHE_PATH, CACHE_SHA256, RESULT_ROOT,
# PACKAGE_ROOT, PACKAGE_COMMIT. Optional: SINGULARITY_IMAGE and R_LIBRARY.

set -euo pipefail

: "${TASK_FILE:?TASK_FILE is required}"
: "${CACHE_PATH:?CACHE_PATH is required}"
: "${CACHE_SHA256:?CACHE_SHA256 is required}"
: "${RESULT_ROOT:?RESULT_ROOT is required}"
: "${PACKAGE_ROOT:?PACKAGE_ROOT is required}"
: "${PACKAGE_COMMIT:?PACKAGE_COMMIT is required}"

task_id="${TASK_ID:-${SGE_TASK_ID:-}}"
case "$task_id" in
  ''|*[!0-9]*) echo "TASK_ID/SGE_TASK_ID must be a positive integer" >&2; exit 64 ;;
esac
if [ "$task_id" -lt 1 ]; then
  echo "TASK_ID/SGE_TASK_ID must be a positive integer" >&2
  exit 64
fi

line=$(awk -F '\t' -v target="$task_id" 'NR == target + 1 { print; exit }' "$TASK_FILE")
if [ -z "$line" ]; then
  echo "Task id is absent from task table: $task_id" >&2
  exit 65
fi
IFS=$'\t' read -r recorded_id unit seed label outer_k <<< "$line"
if [ "$recorded_id" != "$task_id" ]; then
  echo "Task table id/order mismatch" >&2
  exit 65
fi
case "$unit" in ''|*[!A-Za-z0-9_]*) echo "Unsafe unit id" >&2; exit 65 ;; esac
case "$label" in ''|*[!A-Za-z0-9_]*) echo "Unsafe task label" >&2; exit 65 ;; esac
case "$seed" in 101|202|303|404|505) ;; *) echo "Unexpected seed" >&2; exit 65 ;; esac
case "$outer_k" in 3|5) ;; *) echo "Unexpected outer_k" >&2; exit 65 ;; esac
case "$CACHE_SHA256" in
  *[!0-9a-f]*) echo "CACHE_SHA256 must be lowercase hexadecimal" >&2; exit 65 ;;
esac
case "$PACKAGE_COMMIT" in
  *[!0-9a-f]*) echo "PACKAGE_COMMIT must be lowercase hexadecimal" >&2; exit 65 ;;
esac
if [ "${#CACHE_SHA256}" -ne 64 ] || [ "${#PACKAGE_COMMIT}" -ne 40 ]; then
  echo "Cache/package pins must contain 64/40 characters" >&2
  exit 65
fi

package_root=$(realpath "$PACKAGE_ROOT")
cache_path=$(realpath "$CACHE_PATH")
result_root=$(realpath "$RESULT_ROOT")
script="$package_root/inst/scripts/paper3_nested_selector_unit.R"
test -f "$script"
test "$(git -C "$package_root" rev-parse HEAD)" = "$PACKAGE_COMMIT"
test -z "$(git -C "$package_root" status --porcelain=v1 --untracked-files=all)"
test "$(sha256sum "$cache_path" | awk '{print $1}')" = "$CACHE_SHA256"

output_dir="$result_root/$label"
if [ -e "$output_dir" ]; then
  echo "Task output already exists: $output_dir" >&2
  exit 66
fi

args=(
  "$script"
  "--cache=$cache_path"
  "--cache-sha256=$CACHE_SHA256"
  "--unit=$unit"
  "--seed=$seed"
  "--output-dir=$output_dir"
  "--package-commit=$PACKAGE_COMMIT"
  "--inner-folds=5"
  "--bootstrap-reps=1000"
  "--verify-base=true"
)

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ -n "${SINGULARITY_IMAGE:-}" ]; then
  singularity_image="$SINGULARITY_IMAGE"
  test -f "$singularity_image"
  # Singularity 2.x reserves SINGULARITY_IMAGE for its internal launcher. Do
  # not leak the caller-facing path variable into the container environment.
  unset SINGULARITY_IMAGE
  command=(singularity exec "$singularity_image" env)
  if [ -n "${R_LIBRARY:-}" ]; then
    command+=("R_LIBS_USER=$R_LIBRARY")
  fi
  command+=(Rscript)
else
  if [ -n "${R_LIBRARY:-}" ]; then export R_LIBS_USER="$R_LIBRARY"; fi
  command=(Rscript)
fi

"${command[@]}" "${args[@]}"
