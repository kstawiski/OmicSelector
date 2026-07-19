#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
entrypoint="run"
attestation_target=""
runtime_image=""
expected_image_sha256=""
expected_launcher_sha256=""
forwarded=()

while (($#)); do
  case "$1" in
    --entrypoint)
      (($# >= 2)) || { echo "--entrypoint requires a value" >&2; exit 2; }
      entrypoint=$2
      shift 2
      ;;
    --capture-runtime-attestation)
      (($# >= 2)) || {
        echo "--capture-runtime-attestation requires a path" >&2
        exit 2
      }
      attestation_target=$2
      shift 2
      ;;
    --runtime-image)
      (($# >= 2)) || { echo "--runtime-image requires a path" >&2; exit 2; }
      runtime_image=$2
      forwarded+=("$1" "$2")
      shift 2
      ;;
    --expected-runtime-image-sha256)
      (($# >= 2)) || {
        echo "--expected-runtime-image-sha256 requires a value" >&2
        exit 2
      }
      expected_image_sha256=$2
      forwarded+=("$1" "$2")
      shift 2
      ;;
    --expected-container-launcher-sha256)
      (($# >= 2)) || {
        echo "--expected-container-launcher-sha256 requires a value" >&2
        exit 2
      }
      expected_launcher_sha256=$2
      forwarded+=("$1" "$2")
      shift 2
      ;;
    *)
      forwarded+=("$1")
      shift
      ;;
  esac
done

[[ "$entrypoint" == "run" || "$entrypoint" == "assemble" ||
   "$entrypoint" == "capture-attestation" ]] || {
  echo "--entrypoint must be run, assemble, or capture-attestation" >&2
  exit 2
}
[[ -f "$runtime_image" ]] || { echo "runtime image is absent" >&2; exit 2; }
[[ "$expected_image_sha256" =~ ^[0-9a-f]{64}$ ]] || {
  echo "runtime image SHA-256 pin is invalid" >&2
  exit 2
}
[[ "$expected_launcher_sha256" =~ ^[0-9a-f]{64}$ ]] || {
  echo "container launcher SHA-256 pin is invalid" >&2
  exit 2
}

observed_image_sha256=$(sha256sum -- "$runtime_image" | awk '{print $1}')
observed_launcher_sha256=$(sha256sum -- "${BASH_SOURCE[0]}" | awk '{print $1}')
[[ "$observed_image_sha256" == "$expected_image_sha256" ]] || {
  echo "runtime image differs from its exact pin" >&2
  exit 3
}
[[ "$observed_launcher_sha256" == "$expected_launcher_sha256" ]] || {
  echo "container launcher differs from its exact pin" >&2
  exit 3
}

container_runtime=$(command -v singularity || true)
if [[ -z "$container_runtime" ]]; then
  container_runtime=$(command -v apptainer || true)
fi
[[ -n "$container_runtime" ]] || {
  echo "Singularity or Apptainer is required" >&2
  exit 4
}

case "$entrypoint" in
  run)
    entry_script="$script_dir/paper3_run_ranktree_sensitivity_bundle.R"
    entry_args=("${forwarded[@]}")
    ;;
  assemble)
    entry_script="$script_dir/paper3_assemble_ranktree_sensitivity.R"
    entry_args=("${forwarded[@]}")
    ;;
  capture-attestation)
    [[ -n "$attestation_target" ]] || {
      echo "capture-attestation requires --capture-runtime-attestation" >&2
      exit 2
    }
    entry_script="$script_dir/paper3_run_ranktree_sensitivity_bundle.R"
    entry_args=("--write-runtime-attestation" "$attestation_target")
    ;;
esac
[[ -f "$entry_script" ]] || { echo "package entry script is absent" >&2; exit 4; }

export SINGULARITYENV_OMICSELECTOR_RUNTIME_IMAGE_SHA256="$expected_image_sha256"
export SINGULARITYENV_OMICSELECTOR_CONTAINER_LAUNCHER_SHA256="$expected_launcher_sha256"
export APPTAINERENV_OMICSELECTOR_RUNTIME_IMAGE_SHA256="$expected_image_sha256"
export APPTAINERENV_OMICSELECTOR_CONTAINER_LAUNCHER_SHA256="$expected_launcher_sha256"
exec "$container_runtime" exec --cleanenv "$runtime_image" \
  Rscript "$entry_script" "${entry_args[@]}"
