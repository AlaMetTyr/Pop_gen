#!/bin/bash -e
#SBATCH --job-name=r
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --output=fst.out
#SBATCH --error=fst.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##Module load
module load VCFtools 

vcftools --gzvcf FAW_geno-filt.vcf.gz --weir-fst-pop pop_faw-inv.txt --weir-fst-pop pop_faw-nat.txt --fst-window-size 100000 --out --out inv_vs_nat

