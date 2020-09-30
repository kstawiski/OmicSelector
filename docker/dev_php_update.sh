#!/bin/bash

# screen docker run --name OmicSelector -p 28888:80 -p 28889:8888 OmicSelector
docker cp /home/konrad/snorlax/OmicSelector/static/. OmicSelector:/OmicSelector/OmicSelector/static/
docker cp /home/konrad/snorlax/OmicSelector/docker/. OmicSelector:/OmicSelector/OmicSelector/docker/
docker cp /home/konrad/snorlax/OmicSelector/templetes/. OmicSelector:/OmicSelector/OmicSelector/templetes/
