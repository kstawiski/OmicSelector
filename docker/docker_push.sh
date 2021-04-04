#!/bin/bash


# start command: ./docker_push.sh | tee build.log
# must start in ./docker/ dir
# docker login --username=kstawiski
docker builder prune
docker image prune -a

docker build --rm --force-rm -f ../Dockerfile.gpu -t omicselector-gpu ../
# if low memory machine: docker build --rm --force-rm -f ../Dockerfile.workflow -t OmicSelector ../
docker tag omicselector-gpu:latest kstawiski/omicselector-gpu:latest
docker push kstawiski/omicselector-gpu

docker build --rm --force-rm -f ../Dockerfile.cpu -t omicselector ../
# if low memory machine: docker build --rm --force-rm -f ../Dockerfile.workflow -t OmicSelector ../ 
docker tag omicselector:latest kstawiski/omicselector:latest
docker push kstawiski/omicselector

docker build --rm --force-rm -f ../Dockerfile.public -t omicselector-public ../
# if low memory machine: docker build --rm --force-rm -f ../Dockerfile.workflow -t OmicSelector ../
docker tag omicselector-public:latest kstawiski/omicselector-public:latest
docker push kstawiski/omicselector-public

# for google cloud
# docker tag kstawiski/OmicSelector:latest gcr.io/konsta/OmicSelector:latest
# docker push gcr.io/konsta/OmicSelector
docker pull kstawiski/omicselector


docker run --gpus all --name OmicSelector --rm -d -p 28888:80 -v /OmicSelectorData/:/OmicSelector/ kstawiski/omicselector-gpu


# RAMDISK:
# sudo mkdir /mnt/ramdisk
# sudo mount -t tmpfs -o size=128g tmpfs /mnt/ramdisk
