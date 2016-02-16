#!/bin/bash
#PBS -N STM-turnover
#PBS -l walltime=12:00:00
#PBS -l nodes=4:ppn=1
#PBS -q qwork@mp2
#PBS -j oe

module load bioinformatics/R/3.2.3

cd $PBS_O_WORKDIR
RScript ./Turnover.r
