#!/bin/bash -e
#SBATCH --job-name=fst_loop
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --output=fstindiv.out
#SBATCH --error=fstindiv.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##Module load
module load VCFtools
module load BCFtools
module load SAMtools
module load PLINK/2.00a5.14-GCC-12.3.0  

#vcftools --gzvcf FAW_geno-filt.vcf.gz --weir-fst-pop pop_faw-inv-aus.txt --weir-fst-pop pop_faw-inv-nz.txt --fst-window-size 50000 --out inv_vs_nat-corn-ausnz

#plink2 --bfile faw_bialSNP.vcf --fst CATPHENO --within popmap.txt --out fst_results


#awk '{print $2}' pop_faw-Copy2.txt | sort | uniq > pop_list.txt


#while read pop1; do
#  while read pop2; do
#    if [[ "$pop1" != "$pop2" ]]; then
#      awk -v p="$pop1" '$2 == p {print $1}' popmap.txt > pop1_nel.txt
 #     awk -v p="$pop2" '$2 == p {print $1}' popmap.txt > pop2_samples.txt
#
#      vcftools --vcf faw_bialSNP_v2.vcf --weir-fst-pop pop1_samples.txt --weir-fst-pop pop2_samples.txt -fst-window-size 50000 --out "Fst_${pop1}_${pop2}"
#    fi
#  done < pop_list.txt
#done < pop_list.txt

########## with 2 individuals using plinlk2#########


#bgzip faw_bialSNP.vcf > faw_bialSNP.vcf.gz
##then...###
#bcftools view -s FAW-wc-01,Nel-03 -o temporal.vcf.gz -Oz faw_bialSNP.vcf.gz
#bcftools index temporal.vcf.gz


#plink2 --vcf temporal.vcf.gz --make-bed --out faw_bialSNP --allow-extra-chr

plink2 --bfile faw_bialSNP \
        --within pop1_samples.txt \
       --fst CATPHENO method=hudson \
       report-variants \
       blocksize=20000 \
       --out fst_results_nel-WC \
       --allow-extra-chr \
  # --pheno-name POP 


#plink2 --bfile faw_bialSNP \
#       --pheno pop1_samples.txt \
 #      --pheno-name POP \
 #      --fst method=hudson \
 #      report-variants \
 #      blocksize=20000 \
 ##      --out fst_results \
 #      --allow-extra-chr
