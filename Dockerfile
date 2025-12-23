# OmicSelector 2.0 - Unified WebUI + RStudio Server
#
# Runs the Shiny WebUI and RStudio Server in a single container.

FROM rocker/rstudio:4.4.0

LABEL maintainer="OmicSelector Team"
LABEL org.opencontainers.image.title="OmicSelector 2.0"
LABEL org.opencontainers.image.description="Shiny WebUI + RStudio Server for OmicSelector"
LABEL org.opencontainers.image.version="2.0.0"

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV TZ=UTC
ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV RENV_CONFIG_AUTOLOADER_ENABLED=false
ENV OMICSELECTOR_ANALYSES_DIR=/OmicSelector/analyses
ENV CRAN_REPO=https://packagemanager.posit.co/cran/__linux__/jammy/latest
ENV MLR_REPO=https://mlr-org.r-universe.dev
ENV TORCH_HOME=/usr/local/lib/torch

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    make \
    cmake \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libtiff-dev \
    libjpeg-dev \
    libpng-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libgit2-dev \
    libglpk-dev \
    libgmp-dev \
    libgsl-dev \
    libsodium-dev \
    libcairo2-dev \
    libxt-dev \
    libx11-dev \
    pandoc \
    git \
    wget \
    curl \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "options(repos = c(CRAN = Sys.getenv('CRAN_REPO'), MLR = Sys.getenv('MLR_REPO')), Ncpus = 2); install.packages(c('remotes','shiny','shinyjs','bslib','bsicons','shinyWidgets','sass','plotly','mice','readxl','ggplot2','ggrepel','digest'), dependencies = TRUE)"

RUN mkdir -p /OmicSelector /OmicSelector/analyses
RUN mkdir -p ${TORCH_HOME}

COPY . /OmicSelector

WORKDIR /OmicSelector

RUN rm -rf renv renv.lock .Rprofile 2>/dev/null || true

RUN R -e "options(repos = c(CRAN = Sys.getenv('CRAN_REPO'), MLR = Sys.getenv('MLR_REPO')), Ncpus = 2); remotes::install_deps(dependencies = NA)"
RUN R -e "options(repos = c(CRAN = Sys.getenv('CRAN_REPO'), MLR = Sys.getenv('MLR_REPO')), Ncpus = 2); install.packages(c('mlr3learners','ranger','glmnet','xgboost','lightgbm','e1071','kknn','nnet','DALEX','iml','vetiver','pins','plumber','jsonlite','cachem','memoise'))"
RUN R -e "options(repos = c(CRAN = Sys.getenv('CRAN_REPO'), MLR = Sys.getenv('MLR_REPO')), Ncpus = 2); install.packages(c('torch','mlr3torch'))"
RUN R -e 'if (requireNamespace("torch", quietly = TRUE)) { sys <- tolower(Sys.info()[["sysname"]]); arch <- R.version$arch; if (sys == "linux" && arch %in% c("x86_64","amd64")) { torch::install_torch(type = "cpu") } else { message("Skipping torch::install_torch() for ", sys, "/", arch) } }'
RUN R CMD INSTALL --no-multiarch --with-keep.source .

RUN chown -R rstudio:rstudio /OmicSelector ${TORCH_HOME}

COPY docker/start.sh /usr/local/bin/start-omicselector.sh
RUN chmod +x /usr/local/bin/start-omicselector.sh

EXPOSE 3838 8787

CMD ["/usr/local/bin/start-omicselector.sh"]
