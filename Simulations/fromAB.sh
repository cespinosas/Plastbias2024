#!/bin/bash

#PBS -l walltime=30:00:00
#PBS -l nice=19
#PBS -l nodes=1:ppn=1
#PBS -q linux-spool

cd $PBS_O_WORKDIR
./carreras ${PBS_ARRAYID} 12 36 500 0.002 40000 0.02 0.4 1
