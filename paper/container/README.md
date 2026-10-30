# Reproduction environments

Import the checksum-verified environment assets, then build the CPU and CUDA images from the pinned base and OmicSelector 2.6.10 source archive.

- `python3 container/import_base.py --assets ../bundle/container`.
- `bash container/build.sh --bundle ../bundle`.
- `sudo apptainer build omicselector-cpu.sif container/Apptainer.cpu.def`.
- `sudo apptainer build omicselector-cuda.sif container/Apptainer.cuda.def`.
- `apptainer exec --cleanenv omicselector-cpu.sif Rscript -e 'stopifnot(packageVersion("OmicSelector") == "2.6.10")'`.
- `apptainer exec --nv --cleanenv omicselector-cuda.sif Rscript -e 'stopifnot(packageVersion("OmicSelector") == "2.6.10")'`.

Run the summary command inside the CPU environment with the paper and bundle directories mounted. For Docker, mount their common parent at `/study`, set the working directory to `/study/paper`, and run `bash scripts/run_all_summary.sh --bundle /study/bundle`. For Apptainer, bind the common parent at `/study` and use the same command. Set `OMICSELECTOR_PAPER_ROOT` to a writable work directory.

The base pins R 4.5.2, Python 3.12.3 and PyTorch 2.11 with CUDA 12.8. Statistical summary and secondary-analysis producers use the deposited R 4.6.1 runtime and libraries. Build and runtime require Linux x86-64. GPU execution requires a compatible NVIDIA driver. Apptainer 1.4.5 recipe build times were 1237 and 1196 seconds for the CPU and CUDA environments. The five-method replay comprises 15 fold comparisons. Its maximum foundation-model difference is 7.15255737304688e-7 against a tolerance of 1e-6, with unchanged model-state fingerprints.

Docker must permit the `faccessat2` system call or return `ENOSYS` so the Ubuntu 24.04 base can use its compatibility path. Use a Docker engine with a compatible syscall policy. `DOCKER_HOST` can select a dedicated compatible engine. The environment checker validates the pinned R, Python and package versions.
