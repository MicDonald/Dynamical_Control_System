#!/bin/sh

#PBS -l nodes=node14:ppn=1
#PBS -o output.txt
#PBS -j oe
#PBS -N ABC

nproc=`cat $PBS_NODEFILE | wc -l`
nnode=`cat $PBS_NODEFILE | sort | uniq |wc -l`
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
echo Time is `date`
echo Directory is `pwd`
echo $PBS_JOBID $PBS_JOBNAME on $PBS_O_HOST
echo Nodefile: $PBS_NODEFILE
echo This jobs runs on the following nodes:
echo `cat $PBS_NODEFILE`
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cores
echo =================================================
./console
echo =================================================
echo Ending time is `date`
