#!/bin/bash -e
#SBATCH --job-name=vcfcheck
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --output=vcfcheck-%j.out
#SBATCH --error=vcfcheck-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##modules
module purge
module load PLINK/1.09b6.16 
VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

plink --vcf faw_MAF_geno_LD.vcf --freq --missing --hardy --out snp_summary --verbose

vcftools --vcf faw_MAF_geno_LD.vcf --site-mean-depth --out snp_depth_summary --vcf
vcftools --vcf faw_MAF_geno_LD.vcf --site-quality --out snp_quality_summary --vcf
vcftools --vcf faw_MAF_geno_LD.vcf --site-pi --out snp_diversity_summary --vcf

