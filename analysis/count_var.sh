#!/bin/bash -e
#SBATCH --job-name=rg
#SBATCH --time=24:00:00
#SBATCH --mem=5GB
#SBATCH --output=rg-%j.out
#SBATCH --error=rg-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

# Load required modules
module purge
module load GATK/4.5.0.0-gimkl-2022a

 gatk CountVariants \
      -V FAW_geno-filt.vcf.gz
      -O all_Count
 