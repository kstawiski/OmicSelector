#!/usr/bin/env bash

# Run one task from paper3_build_nested_selector_tasks.R. Required environment:
# TASK_FILE, TASK_ID (or SGE_TASK_ID), CACHE_PATH, CACHE_SHA256, RESULT_ROOT,
# PACKAGE_ROOT, PACKAGE_COMMIT, SINGULARITY_IMAGE, RUNTIME_IMAGE_SHA256,
# R_LIBRARY, R_LIBRARY_SNAPSHOT, and R_LIBRARY_SNAPSHOT_SHA256.

set -euo pipefail

: "${TASK_FILE:?TASK_FILE is required}"
: "${CACHE_PATH:?CACHE_PATH is required}"
: "${CACHE_SHA256:?CACHE_SHA256 is required}"
: "${RESULT_ROOT:?RESULT_ROOT is required}"
: "${PACKAGE_ROOT:?PACKAGE_ROOT is required}"
: "${PACKAGE_COMMIT:?PACKAGE_COMMIT is required}"
: "${SINGULARITY_IMAGE:?SINGULARITY_IMAGE is required}"
: "${RUNTIME_IMAGE_SHA256:?RUNTIME_IMAGE_SHA256 is required}"
: "${R_LIBRARY:?R_LIBRARY is required}"
: "${R_LIBRARY_SNAPSHOT:?R_LIBRARY_SNAPSHOT is required}"
: "${R_LIBRARY_SNAPSHOT_SHA256:?R_LIBRARY_SNAPSHOT_SHA256 is required}"

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
case "$RUNTIME_IMAGE_SHA256" in
  *[!0-9a-f]*) echo "RUNTIME_IMAGE_SHA256 must be lowercase hexadecimal" >&2; exit 65 ;;
esac
case "$R_LIBRARY_SNAPSHOT_SHA256" in
  *[!0-9a-f]*) echo "R_LIBRARY_SNAPSHOT_SHA256 must be lowercase hexadecimal" >&2; exit 65 ;;
esac
if [ "${#CACHE_SHA256}" -ne 64 ] || [ "${#PACKAGE_COMMIT}" -ne 40 ] ||
   [ "${#RUNTIME_IMAGE_SHA256}" -ne 64 ] ||
   [ "${#R_LIBRARY_SNAPSHOT_SHA256}" -ne 64 ]; then
  echo "Cache/runtime/library/package pins must contain 64/64/64/40 characters" >&2
  exit 65
fi

package_root=$(realpath "$PACKAGE_ROOT")
cache_path=$(realpath "$CACHE_PATH")
result_root=$(realpath "$RESULT_ROOT")
singularity_image=$(realpath "$SINGULARITY_IMAGE")
r_library=$(realpath "$R_LIBRARY")
r_library_snapshot=$(realpath "$R_LIBRARY_SNAPSHOT")
script="$package_root/inst/scripts/paper3_nested_selector_unit.R"
test -f "$script"
test -f "$singularity_image"
test -d "$r_library"
test -f "$r_library_snapshot"
test "$(git -C "$package_root" rev-parse HEAD)" = "$PACKAGE_COMMIT"
test -z "$(git -C "$package_root" status --porcelain=v1 --untracked-files=all)"
test "$(sha256sum "$cache_path" | awk '{print $1}')" = "$CACHE_SHA256"
test "$(sha256sum "$singularity_image" | awk '{print $1}')" = \
  "$RUNTIME_IMAGE_SHA256"
test "$(sha256sum "$r_library_snapshot" | awk '{print $1}')" = \
  "$R_LIBRARY_SNAPSHOT_SHA256"

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
  "--runtime-image-sha256=$RUNTIME_IMAGE_SHA256"
  "--r-library-snapshot-sha256=$R_LIBRARY_SNAPSHOT_SHA256"
  "--inner-folds=5"
  "--bootstrap-reps=1000"
  "--verify-base=true"
)

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Singularity 2.x reserves SINGULARITY_IMAGE for its internal launcher. Do not
# leak the caller-facing path variable into the container environment.
unset SINGULARITY_IMAGE
command=(
  singularity exec "$singularity_image" env "R_LIBS_USER=$r_library" Rscript
)

"${command[@]}" "${args[@]}"
