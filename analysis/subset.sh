#!/bin/bash -e
#SBATCH --job-name=bcf
#SBATCH --time=48:00:00
#SBATCH --mem=4GB
#SBATCH --output=bcf-%j.out
#SBATCH --error=bcf-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

# Load required modules
module purge
module load BCFtools

bcftools view -S indiv.txt FAW_geno-filt.vcf.gz -O z -o subset.filtered.vcf.gz