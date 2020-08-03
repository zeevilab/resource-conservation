#!/bin/bash
set -ex



Rscript "LMM_step_1.R" $2 $3
# chmod +x KEGG/LMM_step_2.sh
# KEGG/LMM_step_2.sh $1 
# Rscript "KEGG/LMM_step_3.R" $2 $3



