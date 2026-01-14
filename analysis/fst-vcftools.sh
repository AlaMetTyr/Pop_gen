#!/bin/bash -e
#SBATCH --job-name=fst
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --output=fst.out
#SBATCH --error=fst.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##Module load
module load VCFtools 

#!/bin/bash

VCF="/home/a.vaughan/nobackup_ga03488/Amy/FAW/2025/cohort/final_maf_filtered/faw_geno_LD_maf_v42.vcf"     # <--- CHANGE THIS TO YOUR VCF
POP_DIR="./pop"               # folder containing pop_*.txt files

# Loop through all population files
pops=($POP_DIR/pop_*.txt)

for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        
        pop1=${pops[$i]}
        pop2=${pops[$j]}

        # clean names (pop_Cambodia.txt → Cambodia)
        n1=$(basename "$pop1" .txt | sed 's/^pop_//')
        n2=$(basename "$pop2" .txt | sed 's/^pop_//')

        out="FST_${n1}_vs_${n2}"
        
        echo "Running FST for: $n1 vs $n2"
        
        vcftools --vcf "$VCF" \
            --weir-fst-pop "$pop1" \
            --weir-fst-pop "$pop2" \
            --out "$out"
    done
done
