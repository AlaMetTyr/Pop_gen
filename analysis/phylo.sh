#!/bin/bash -e
#SBATCH --job-name=phylo
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --output=relate-phylo.out
#SBATCH --error=relate-phylo.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##Module load
module load PLINK 

plink --bfile faw_LD_maf \
       --distance square 1-ibs \
       --allow-extra-chr \
       --out faw_LD_maf_plink1
