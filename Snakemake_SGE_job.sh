#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q medium.q
#$ -N sm_final
#$ -pe smp 56

module load singularity/3.6.4
module load anaconda3/2020.02
#conda activate snakemake_test
conda activate TEnriAn

#one core will be used by snakemake to monitore the other processes
THREADS=$(expr ${NSLOTS} - 1)

#snakemake --cores ${THREADS} --until run_assembly --use-conda --use-singularity --verbose
#snakemake --cores ${THREADS} --until run_orthology --use-conda --use-singularity --verbose
snakemake --cores ${THREADS} --until run_alignment_filtering --use-conda --use-singularity --verbose

### when setting up the workflow for the first time, we need to create the environments:
### snakemake --use-conda --conda-create-envs-only --jobs 1
### OR
###  snakemake --use-singularity --use-conda --conda-create-envs-only --jobs 1

###  If job crashed or was stopped you need to first unlock the directories with 'snakemake --unlock'
