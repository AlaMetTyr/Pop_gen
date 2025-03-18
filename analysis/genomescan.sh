#!/bin/bash -e
#SBATCH --job-name=reshaper

module purge
module load BayPass/2.31-intel-2022a
module load VCFtools

##subset vcf to region for analysis
bcftools view -r "Scaffold HiC_3" faw_bialSNP_MAF_genov2.vcf.gz -o faw_HiC_3_subset.vcf

##reorder and then output in baypass format using reshaper script
vcftools --vcf faw_HiC_3_subset.vcf --keep popmap_reshape.txt --out reordered --recode
python ./reshaper_baypass.py reordered.recode.vcf popmap_reshape.txt FAW.geno

##baypass contrast analysis
i_baypass -gfile FAW.geno -contrastfile contrast.ecotype -outprefix FAW_baypass -nthreads 2
