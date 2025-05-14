#!/bin/bash -e
#SBATCH --job-name=treemix_stacks
#SBATCH --time=24:00:00
#SBATCH --mem=30GB
#SBATCH --output=tm.out
#SBATCH --error=tm.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##Module load
module load Stacks 
module load TreeMix/1.13-GCC-11.3.0

populations -V ../reordered.recode.vcf -O ./ -M pop_faw-reshape.txt --treemix
tail -n +2 reordered_faw.p.treemix > faw_4pop_rice2.treemix
gzip faw_4pop_rice2.treemix

module load Miniconda3
for i in {1..6};
do
treemix -i faw_4pop_rice2.treemix.gz -root rice -m ${i} -bootstrap -k 1000 -o 4pop_rice2_${i} > 4pop_rice2_${i}_log 
done