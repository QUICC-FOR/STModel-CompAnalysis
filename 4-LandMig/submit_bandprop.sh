#!/bin/bash
#PBS -N STM-bandprop
#PBS -l walltime=24:00:00
#PBS -l nodes=12:ppn=1
#PBS -q qfbb
#PBS -j oe

cd $PBS_O_WORKDIR
Rscript ./2_prop2095.r
