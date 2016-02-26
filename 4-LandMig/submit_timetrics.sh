#!/bin/bash
#PBS -N STM-timetrics
#PBS -l walltime=24:00:00
#PBS -l nodes=20:ppn=1
#PBS -q qfbb
#PBS -j oe

cd $PBS_O_WORKDIR
Rscript ./1_timetrics.r
