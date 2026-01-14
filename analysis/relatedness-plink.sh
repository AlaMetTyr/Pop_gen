#!/bin/bash -e
#SBATCH --job-name=plink-maf
#SBATCH --account=ga03488    
#SBATCH --time=10:00:00        
#SBATCH --mem=10G             
#SBATCH -o logs/plink_maf%A_%a.out
#SBATCH -e logs/plink_maf%A_%a.err

module purge
module load PLINK
module load BCFtools


VCF="../../final_maf_filtered/faw_geno_LD_maf_v42.vcf"       
PREFIX="plink2_output"       


plink2 --vcf "$VCF" \
    --make-pgen \
    --out "$PREFIX"


plink2 --pfile "$PREFIX" \
    --make-king-table \
    --out "${PREFIX}_king_relatedness"
