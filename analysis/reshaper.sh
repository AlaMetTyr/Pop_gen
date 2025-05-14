#!/bin/bash -e
#SBATCH --job-name=reshaper
#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=12
#SBATCH --output=reshaper.out
#SBATCH --error=reshaper.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488


module load VCFtools
vcftools --vcf ausnz.vcf --keep popmap_reshape-ausnz.txt --out ausnz --recode

python ./reshaper_baypass.py ausnz.recode.vcf popmap_reshape-ausnz.txt ./baypass/ausnz.geno