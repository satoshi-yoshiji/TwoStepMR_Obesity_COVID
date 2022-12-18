#!/bin/bash
#PBS -N MR
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=24G
#PBS -l vmem=24G
#PBS -t 1-491

mkdir -p output/harmonized
mkdir -p output/or
mkdir -p output/pleio
mkdir -p output/hetero

cd $PBS_O_WORKDIR
mkdir -p logs/

list=$(head -n ${PBS_ARRAYID} /scratch/richards/satoshi.yoshiji/11.pQTL/decode_batch/listbatch.txt | tail -n 1) # path to pQTLs list
for protein in `cat $list`;
do
protein_name=`basename $protein ".txt.gz"`

/admin/kleinman/custom_software/R-4.1.2-gcc-9.3.0/lib64/R/bin/Rscript --vanilla 01.MR.R $protein $protein_name
done
