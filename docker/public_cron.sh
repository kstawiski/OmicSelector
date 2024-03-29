#!/bin/bash
# crontab entry: 0 5 * * 1 /root/omicselector_cron.sh
docker stop omicselector-public
docker pull kstawiski/omicselector-public
docker run --name omicselector-public --rm --gpus all --cpus="12" --memory="32g" --memory-swap="32g" --env PUBLIC=1 -d -p 20019:80 kstawiski/omicselector-public