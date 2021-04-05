![](vignettes/logo.png)

# OmicSelector

![Docker](https://github.com/kstawiski/OmicSelector/workflows/Docker/badge.svg) ![R package](https://github.com/kstawiski/OmicSelector/workflows/R%20package/badge.svg)

OmicSelector is the environment, docker-based application and R package for biomarker signiture selection (feature selection) & deep learning diagnostic tool development from high-throughput high-throughput omics experiments and other multidimensional datasets. It was initially developed for miRNA-seq (small RNA, smRNA-seq; hence the previous name was miRNAselector), RNA-seq and qPCR, but can be applied for every problem where numeric features should be selected to counteract overfitting of the models. Using our tool, you can choose features, like miRNAs, with the most significant diagnostic potential (based on the results of miRNA-seq, for validation in qPCR experiments). It can also develop the best deep learning model for your signature, as well as be an IDE for your more complex data mining project (contains R Studio, Jupyter notebooks and VS Code.. all integrated in one!).

The main purpose of OmicSelector is to provide you with the set of candidate features (biomarkers) for further validation of biomarker study from e.g. high-throughput experiments. The package performs feature selection first. In the next step, the sets of features are tested in the process called “benchmarking”. In benchmarking we test all of those sets of features (biomarkers) using various data-mining (machine learning) methods. Based on the average performance of sets in cross-validation or holdout-validation (testing on the test set and/or validation set) we can suggest which of the signatures (set of features) has the greatest potential in further validation. As the feautres are selected, OmicSelector can perform advanced modeling of deep feedforward neural networks with and without autoencoders. The best network is developed using comperhensive grid search of optimal hyperparameters. This section works with Tensorflow (via Keras), so the computations can be GPU-accelerated! The best network can be easily implemented in clinical practice using our interactive application.

Go to https://biostat.umed.pl/OmicSelector/ for more details.

## Try it out

[![Run on Ainize](https://ainize.ai/images/run_on_ainize_button.svg)](https://ainize.web.app/redirect?git_repo=https://github.com/kstawiski/OmicSelector)

Link: [https://master-omic-selector-kstawiski.endpoint.ainize.ai/](https://master-omic-selector-kstawiski.endpoint.ainize.ai/)

Please note that uploading real data to this instance is not safe. You're data and analysis files will be accessible by anyone. You can fork the repo if you want your own free tesing copy using [Ainize](https://ainize.web.app/redirect?git_repo=https://github.com/kstawiski/OmicSelector).

## Quick start

### Docker (with GUI):

1. GPU-based, using Nvidia CUDA: [kstawiski/omicselector-gpu](https://hub.docker.com/r/kstawiski/omicselector-gpu)

```
docker run --name OmicSelector --restart always -d -p 28888:80 --gpus all -v $(pwd)/:/OmicSelector/host/ kstawiski/omicselector-gpu
```

2. CPU-based: [kstawiski/omicselector](https://hub.docker.com/r/kstawiski/omicselector)

```
docker run --name OmicSelector --restart always -d -p 28888:80 -v $(pwd)/:/OmicSelector/host/ kstawiski/omicselector
```

As docker image updates itself, it may take few minutes for the app to be operational. You can check logs using `docker logs OmicSelector`. The GUI is accessable via `http://your-host-ip:28888/`. If you use command above, your working directory will be binded as `/OmicSelector/host/`.

Video tutorial is available here:

[![OmicSelector - feature selection and deep learnining with GUI - tutorial](https://yt-embed.herokuapp.com/embed?v=dKUdINEcOjk)](https://www.youtube.com/watch?v=dKUdINEcOjk "OmicSelector - feature selection and deep learnining with GUI - tutorial.")

This tutorial shows how OmicSelector' GUI works and how to perform (without programming knowledge):

- Feature selection
- Benchmarking (selecting best set of variables based on the performance of data-mining models)
- Deep learning model development (feedforward neural network up to 3 hidden layers and with/without autoencoders; grid search of hyperparameters)
- Exploratory analysis (differential expression using t-test, imputation of missing data using predictive mean matching, correcting the batch effect using ComBat, generating heatmaps and volcano plots).

### R package (without GUI):

#### Own enviorment:

```
library("devtools") # if not installed, install via install.packages('devtools')
source_url("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/vignettes/setup.R")
install_github("kstawiski/OmicSelector", force = T)
library(keras)
install_keras()
library(OmicSelector)
OmicSelector_setup()
```

#### Linux/macOS using conda

1. Installing the package in own Anaconda enviorment:

Use e.g. `conda create -n OmicSelector` and `conda activate OmicSelector` to set up your enviorment.  Please note the this will work only when running on **Linux (Ubuntu)** OS or macOS.

```
conda install --channel "conda-forge" --channel "anaconda" --channel "r" tensorflow keras jupyter jupytext numpy pandas r r-devtools r-rgl r-rjava r-mnormt r-purrrogress r-xml gxx_linux-64 libxml2 pandoc r-rjava r-magick opencv pkgconfig gfortran_linux-64
echo "options(repos=structure(c(CRAN='http://cran.r-project.org')))" >> ~/.Rprofile
Rscript -e 'update.packages(ask = F); install.packages(c("devtools","remotes")); remotes::install_cran("pkgdown");'
Rscript -e 'devtools::source_url("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/vignettes/setup.R")'
```

If you have compatible GPU you can consider changing `tensorflow` to `tensorflow-gpu` in `mamba install` command.

#### Windows OS

You can download our Windows-based R enviroment from here: https://studumedlodz-my.sharepoint.com/:u:/g/personal/btm_office365_umed_pl/EQUihquz915JoVhsQQShcnoBZaukMkwd3MnC1LER0iORNw?e=W6KEyu 

After unpacking, if you wish to use our enviorment please consider setting the our R version in your [R Studio](https://rstudio.com/products/rstudio/download/) installation:

![](vignettes/win1.png)

![](vignettes/win2.png)

![](vignettes/win3.png)


## Footnote

Citation:

`In press.`

Authors:

- [Konrad Stawiski, M.D. (konrad.stawiski@umed.lodz.pl)](https://konsta.com.pl)
- Marcin Kaszkowiak.
- Damian Mikulski, M.D.

Supervised by: prof. Wojciech Fendler, M.D., Ph.D. 

For any troubleshooting use [https://github.com/kstawiski/OmicSelector/issues](https://github.com/kstawiski/OmicSelector/issues).

Department of Biostatistics and Translational Medicine, Medical Univeristy of Lodz, Poland (https://biostat.umed.pl) 
