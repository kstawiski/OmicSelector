#!/bin/bash
conda init
. ~/.bashrc
Rscript /update.R
cd /OmicSelector/
git clone https://github.com/kstawiski/OmicSelector.git
chmod -R 755 /OmicSelector/OmicSelector/static/
/usr/sbin/nginx -g "daemon off;" &
mkdir -p /run/
mkdir -p /run/php/
php-fpm7.3 -R -F &
jupyter serverextension enable jupytext
jupyter notebook --no-browser