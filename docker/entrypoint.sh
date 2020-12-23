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
screen -dmS permissionsfix find /opt/conda ! -perm 775 -print0 | xargs -0 -I {} chmod 775 {}

# nignx+php
/usr/sbin/nginx -g "daemon off;" &
mkdir -p /run/
mkdir -p /run/php/
php-fpm7.3 -R -F &

# Rstudio server
/rstudio-server-conda/start_rstudio_server.sh 8787 &

# Jupyter
jupyter serverextension enable jupytext
jupyter notebook --no-browser