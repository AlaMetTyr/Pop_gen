#!/bin/bash -e
#SBATCH --job-name=faststructure
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --output=faststructure-%j.out
#SBATCH --error=faststructure-%j.err
#SBATCH --account=ga03488

##modules
module purge
module load PLINK/2.00a2.3
#module load BCFtools/1.19-GCC-11.3.0
module load fastStructure/1.0-gimkl-2020a-Python-2.7.18


# prep plink files
#bcftools query -f '%CHROM\t%POS\t%END\n' biallelic.vcf > faw.bed
#plink2 --vcf ../../maf_pca/faw_geno_LD_maf.vcf --allow-extra-chr --make-bed --out faw


##modules
module load NeSI
module load fastStructure/1.0-gimkl-2020a-Python-2.7.18

##run faststructure
#python /opt/nesi/CS400_centos7_bdw/fastStructure/1.0-gimkl-2020a-Python-2.7.18/bin/structure.py -K 10 --input=faw_MAF_geno_LD --output=faw

#$ python /opt/nesi/CS400_centos7_bdw/fastStructure/1.0-gimkl-2020a-Python-2.7.18/bin/chooseK.py --input=../../maf_pca/faw_geno_LD_maf

for K in {1..10}; do
  structure.py -K $K --input=faw --output=your_output_prefix_K$K --full --seed=100
done

