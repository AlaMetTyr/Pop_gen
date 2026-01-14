#!/bin/bash -e
#SBATCH --job-name=relatedness
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --output=relate.out
#SBATCH --error=relate.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##Module load
module load VCFtools 

vcftools --vcf ../../final_maf_filtered/faw_geno_LD_maf_v42.vcf \
    --relatedness2 \
    --out vcftools_relatedness
