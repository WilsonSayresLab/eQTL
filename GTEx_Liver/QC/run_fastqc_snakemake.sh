#!/bin/bash
#SBATCH --job-name=FASTQC_snakemake # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=hnatri@asu.edu # send-to address
#SBATCH -t 96:00:00
#SBATCH -n 1

newgrp combinedlab

source activate fastqc_environment

snakemake --snakefile GTEx_Liver_fastqc.snakefile -j 20 --rerun-incomplete --cluster "sbatch -p private -n 1 --nodes 1 -c 8 -t 96:00:00"