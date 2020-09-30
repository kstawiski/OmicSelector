#!/bin/bash
echo "Starting building..." 2>&1 | tee /build.log
cd /OmicSelector/OmicSelector
git pull 2>&1 | tee -a /build.log
cd /OmicSelector
date 2>&1 | tee -a /build.log
R CMD build OmicSelector 2>&1 | tee -a /build.log
R CMD check OmicSelector 2>&1 | tee -a /build.log