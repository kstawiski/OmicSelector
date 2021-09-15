#!/bin/bash
conda init
. ~/.bashrc
Rscript /update.R

# Check public
echo $PUBLIC > /PUBLIC

# GUI
cd /OmicSelector/
git clone https://github.com/kstawiski/OmicSelector.git
chmod -R 755 /OmicSelector/OmicSelector/static/
screen -dmS permissionsfix find /opt/conda ! -perm 777 -print0 | xargs -0 -I {} chmod 777 {}

# nignx+php
/usr/sbin/nginx -g "daemon off;" &
mkdir -p /run/
mkdir -p /run/php/
php-fpm7.3 -R -F &

# vscode
code-server --install-extension ms-python.python --force
code-server --install-extension ikuyadeu.r --force
code-server --install-extension formulahendry.code-runner --force
code-server --install-extension grapecity.gc-excelviewer --force
code-server --install-extension daghostman.vs-treeview --force
screen -dmS vscode code-server --auth none --port 8080

# Rstudio server
rstudio-server start


# Apps:
Rscript -e 'if(!dir.exists("/radiant-data")) { dir.create("/radiant-data") }'
cd /radiant-data/
screen -dmS radiant Rscript -e "radiant::radiant_url(port = 3839)"

# Shiny apps:
# screen -dmS shiny shiny-server
screen -dmS app-deeplearning_model R -e "shiny::runApp('/OmicSelector/OmicSelector/shiny/deeplearning_model', port = 20001)"
screen -dmS app-tool_batch R -e "shiny::runApp('/OmicSelector/OmicSelector/shiny/tool_batch', port = 20002)"
screen -dmS app-tool_heatmap R -e "shiny::runApp('/OmicSelector/OmicSelector/shiny/tool_heatmap', port = 20003)"
screen -dmS app-tool_impute R -e "shiny::runApp('/OmicSelector/OmicSelector/shiny/tool_impute', port = 20004)"
screen -dmS app-tool_de R -e "shiny::runApp('/OmicSelector/OmicSelector/shiny/tool_de', port = 20005)"
screen -dmS app-start R -e "shiny::runApp('/OmicSelector/OmicSelector/shiny/start', port = 20006)"

cd /OmicSelector/
# Jupyter
jupyter serverextension enable jupytext
jupyter notebook --no-browser