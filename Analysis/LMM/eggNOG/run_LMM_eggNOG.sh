#!/bin/bash
set -ex



Rscript "LMM_step_1.R" $1 $2 $3
chmod +x eggNOG/LMM_step_2.sh
eggNOG/LMM_step_2.sh $1 
Rscript "eggNOG/LMM_step_3.R" $1 $2 $3 $4



