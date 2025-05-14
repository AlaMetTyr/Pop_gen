#!/bin/bash -e
#SBATCH --job-name=vcftools
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --output=GATK-table.out
#SBATCH --error=GATKtable.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

# Load required modules
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

vcftools --gzvcf ./FAW_geno-filt.vcf.gz --counts

vcftools --gzvcf ./FAW_geno-filt.vcf.gz --freq --out output

