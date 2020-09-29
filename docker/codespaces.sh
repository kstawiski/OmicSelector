#!/bin/bash
cd /workspace/
wget https://www.php.net/distributions/php-7.4.9.tar.gz
tar xzf php-7.4.9.tar.gz
cd php-7.4.9
./configure --enable-fpm --with-mysqli --with-openssl --with-zlib --with-gd --with-jpeg-dir=/usr/lib/
make
sudo make install
curl -sS https://getcomposer.org/installer -o composer-setup.php && sudo php composer-setup.php --install-dir=/usr/local/bin --filename=composer

sudo curl -s https://packages.sury.org/php/apt.gpg | sudo apt-key add -