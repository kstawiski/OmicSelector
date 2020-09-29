#!/bin/bash
set -e
# Task:
echo "Preprocessing" > /task-name.txt

# Process:
echo "Starting task..." 2>&1 | tee /task.log
date 2>&1 | tee -a /task.log
cd /OmicSelector/
Rscript -e '
    library(rmarkdown);
    file.copy("/OmicSelector/OmicSelector/templetes/1_preprocessing.rmd", "/OmicSelector/1_preprocessing.Rmd", overwrite = TRUE);
    render("/OmicSelector/1_preprocessing.Rmd", output_file = "1_preprocessing.html", output_dir = "/OmicSelector");
    file.copy("/task.log", "/OmicSelector/1_preprocessing.Rmd", overwrite = TRUE);
    cat("[OmicSelector: TASK COMPLETED]"); writeLines("[2] PREPROCESSED", "var_status.txt", sep="");
    ' 2>&1 | tee -a /task.log