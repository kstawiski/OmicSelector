#!/bin/bash
#docker run --gpus all --rm --name OmicSelector -d -p 28888:80 -v /home/konrad/:/OmicSelector/host/ kstawiski/omicselector-gpu
docker run --gpus all --name OmicSelector -d -p 28888:80 -v /home/konrad/:/OmicSelector/host/ --userns=host --privileged kstawiski/omicselector-gpu
