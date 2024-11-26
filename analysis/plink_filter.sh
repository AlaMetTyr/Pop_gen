#!/bin/bash -e
#SBATCH --job-name=plink
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --output=plink-%j.out
#SBATCH --error=plink-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##modules
module purge
module load PLINK/1.09b6.16 
module load BCFtools/1.19-GCC-11.3.0

#filter out multiallelic variants
#bcftools view -m2 -M2 -v snps 2_genotyped.vcf -o biallelic.vcf
#bcftools filter -i 'MAF > 0.05' -O v -o faw_MAF.vcf biallelic.vcf

#plink --vcf faw_MAF.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --geno 0.1 --recode vcf --out faw_MAF_geno
#plink --vcf faw_MAF_geno.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --indep-pairwise 50 5 0.2
#plink --vcf faw_MAF_geno.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --exclude plink.prune.out --recode vcf --out faw_MAF_geno_LD	

# prep plink PCA
plink --vcf biallelic.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out faw_MAF_geno_LD_pca

sed -e 's/ /\t/g' pca.eigenvec > pca.eigenvec

