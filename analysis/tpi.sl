#!/bin/bash -e
#SBATCH --job-name=tpi
#SBATCH --time=4:00:00
#SBATCH --mem=5gb
#SBATCH --account=ga03488
#SBATCH -o tpi_%A_%a.out
#SBATCH -e tpi_%A_%a.err

module load BCFtools/1.19-GCC-11.3.0
module load IQ-TREE
module load SAMtools
module load PLINK
module load VCFtools
module load tabix

#VCF="/nesi/nobackup/ga03488/Amy/FAW/2025/cohort/snps.pass.vcf.gz"   
#REF=/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa
#SCAFF=HiC_scaffold_29
#START=8352258
#END=8355138


awk '/^>/{if(n){print name,len};name=substr($0,2);len=0;n=1;next}{len+=length}END{if(n)print name,len}' TPI.consensus.fasta | column -t

iqtree2 -s TPI.consensus.fasta -m TEST -B 1000 --alrt 1000 -T AUTO

tabix -h ../../final_30M_basefilter/snps.pass.vcf.gz HiC_scaffold_29:8352257-8355137 > TPI.Invasive.vcf

vcftools --vcf TPI.Invasive.vcf --plink --out TPI.invasive

plink2 --pedmap TPI.invasive --make-bed --allow-extra-chr --out TPI.invasive

plink2 --bfile TPI.invasive \
       --mind 0.1 \
       --pca 10 \
       --allow-extra-chr \
       --out TPI.invasive



