#!/bin/bash




# must start in ./docker/ dir
# docker login --username=kstawiski
docker builder prune
docker image prune -a

docker build --rm --force-rm -f ../Dockerfile.gpu -t OmicSelector-gpu ../
# if low memory machine: docker build --rm --force-rm -f ../Dockerfile.workflow -t OmicSelector ../
docker tag OmicSelector-gpu:latest kstawiski/OmicSelector-gpu:latest
docker push kstawiski/OmicSelector-gpu

docker build --rm --force-rm -t OmicSelector ../
# if low memory machine: docker build --rm --force-rm -f ../Dockerfile.workflow -t OmicSelector ../ 
docker tag OmicSelector:latest kstawiski/OmicSelector:latest
docker push kstawiski/OmicSelector

# for google cloud
# docker tag kstawiski/OmicSelector:latest gcr.io/konsta/OmicSelector:latest
# docker push gcr.io/konsta/OmicSelector
docker pull kstawiski/OmicSelector


docker run --name OmicSelector --rm -d -p 28888:80 -v /boot/temp/:/tmp/ -v /home/konrad/:/OmicSelector/host/ kstawiski/OmicSelector


# RAMDISK:
# sudo mkdir /mnt/ramdisk
# sudo mount -t tmpfs -o size=128g tmpfs /mnt/ramdisk
