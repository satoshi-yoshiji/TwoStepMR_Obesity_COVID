#!/bin/bash
#PBS -N MR_proxy
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=24G
#PBS -l vmem=24G
#PBS -t 1-84

mkdir -p output/harmonized
mkdir -p output/or
mkdir -p output/pleio
mkdir -p output/hetero
mkdir -p outout/steiger

cd $PBS_O_WORKDIR
mkdir -p logs/

list=$(head -n ${PBS_ARRAYID} /scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.2.proteins_for_step2_without_reverse/pQTL_batch/listbatch.txt | tail -n 1)
for protein in `cat $list`;
do
protein_name=`basename $protein ".tsv"`  #in the format of "A1BG.16561.9"

/admin/kleinman/custom_software/R-4.1.2-gcc-9.3.0/lib64/R/bin/Rscript --vanilla 01.MR.proxy.R $protein $protein_name
done
