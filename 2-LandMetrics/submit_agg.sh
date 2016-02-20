#!/bin/bash
#PBS -N STM-agg
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -q qwork
#PBS -j oe

module load bioinformatics/R/3.2.3
cd $PBS_O_WORKDIR
Rscript ./2_agg_turnover.r
