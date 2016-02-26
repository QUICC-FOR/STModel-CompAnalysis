#!/bin/bash
#PBS -N STM-solve
#PBS -l walltime=24:00:00
#PBS -l nodes=12:ppn=1
#PBS -q qfbb
#PBS -j oe

module load bioinformatics/R/3.2.3
cd $PBS_O_WORKDIR
Rscript ./1_solve.r
