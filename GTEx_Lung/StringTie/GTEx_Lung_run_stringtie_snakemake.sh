#!/bin/bash
#SBATCH --job-name=StringTie_snakemake # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=dlevy2@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 96:00:00
#SBATCH --qos=normal


source activate stringtie_environment

snakemake --snakefile GTEx_Lung_stringtie.snakefile -j 20 --keep-target-files --rerun-incomplete --cluster "sbatch -p private -n 8 -c 1 -t 96:00:00"