#!/bin/bash
#
#
#PBS -l nodes=1:ppn=8,walltime=36:00:00
# Define working directory
export WORK_DIR=$HOME/main_project/ALSPAC_EWAS/EWAS/
# Change into working directory
cd $WORK_DIR
# Execute code

Trait="metabolite_ewas_phens.txt"
CellData="houseman"
CellAdj="Cells"
BorM="B"
TP="FOM"
Covariates="none"
WD="/panfs/panasas01/sscm/tb13101/main_project/ALSPAC_EWAS/EWAS"

time Rscript R/meffil_EWAS_script.r ${Trait} ${CellData} ${CellAdj} ${BorM} ${TP} ${Covariates} ${WD}