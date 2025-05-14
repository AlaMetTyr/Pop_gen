#!/bin/bash -e
#SBATCH --job-name=baypass
#SBATCH --time=8:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=12
#SBATCH --output=baypass-nz.out
#SBATCH --error=baypass-nz.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488


module purge
module load BCFtools
module load VCFtools
module load BayPass/2.31-intel-2022a

##subset region wantto explore
#bcftools view -r "HiC_scaffold_8" faw_bialSNP_MAF_genov2.vcf.gz -o faw_HiC_8_subset.vcf
#bcftools view -r "HiC_scaffold_8" faw_bialSNP_MAF_genov2.vcf.gz -o faw_HiC_14_subset.vcf
#bcftools view -r "HiC_scaffold_8" faw_bialSNP_MAF_genov2.vcf.gz -o faw_HiC_29_subset.vcf

#vcftools --vcf faw_HiC_8_subset.vcf --keep popmap_reshape.txt --out HiC8 --recode


python ./reshaper_baypass.py ./vcf/NZ_bialSNP_MAF_geno.vcf nz_pops-Copy1.txt ./baypass/nz-filtered.geno


#run baypass on different regions
i_baypass -gfile ./baypass/nz-filtered.geno -contrastfile contrast-nz.ecotype -outprefix nz_NI-SI_baypass -nthreads 2