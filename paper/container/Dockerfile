# Import the checksum-verified E1 base archive with container/import_base.py first.
ARG E1_BASE=omicselector-study-e1-lane:2fa387c6
FROM ${E1_BASE} AS common
LABEL org.opencontainers.image.title="OmicSelector paper reproduction" \
      org.opencontainers.image.revision="b27bec71373ef9cb63aaa270bc0d4f0f6f39dbf2" \
      org.omicselector.e1.source="sha256:2fa387c63144aed56acdd867e81de36aa192346f183dc36f19d6cddd2acfc030"
COPY container/E1-python-freeze.txt container/E1-r-packages.tsv container/check_environment.py /opt/paper-locks/
RUN /opt/venv/bin/python /opt/paper-locks/check_environment.py
COPY container/OmicSelector-source.tar.gz container/package-source.json /opt/paper-locks/
RUN echo "967009d63d9cb4891666f2a086f8515c19f3944a0abf4f7ae0e3390cc1d54273  /opt/paper-locks/OmicSelector-source.tar.gz" | sha256sum -c - \
 && mkdir -p /opt/paper-package \
 && tar -xzf /opt/paper-locks/OmicSelector-source.tar.gz -C /opt/paper-package \
 && R CMD INSTALL /opt/paper-package/OmicSelector \
 && Rscript --vanilla -e 'stopifnot(as.character(packageVersion("OmicSelector")) == "2.6.10")'
ENV OMICSELECTOR_PACKAGE_SOURCE=/opt/paper-package/OmicSelector
ENV RETICULATE_PYTHON=/opt/venv/bin/python RETICULATE_USE_MANAGED_VENV=no PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 NUMBA_NUM_THREADS=4
ENV ATEN_CPU_CAPABILITY=avx2 OPENBLAS_CORETYPE=Zen MKL_ENABLE_INSTRUCTIONS=AVX2 ONEDNN_MAX_CPU_ISA=AVX2 DNNL_MAX_CPU_ISA=AVX2
ENV HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1 TABPFN_DISABLE_TELEMETRY=1 HF_HUB_DISABLE_IMPLICIT_TOKEN=1
ENV R_ENABLE_JIT=0 R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
WORKDIR /workspace
FROM common AS cpu
ENV CUDA_VISIBLE_DEVICES=-1
RUN /opt/venv/bin/python -c "import torch; assert not torch.cuda.is_available(); print(torch.__version__)"
FROM common AS cuda
ENV NVIDIA_VISIBLE_DEVICES=all NVIDIA_DRIVER_CAPABILITIES=compute,utility
RUN /opt/venv/bin/python -c "import torch; assert torch.version.cuda == '12.8'; print(torch.__version__)"
