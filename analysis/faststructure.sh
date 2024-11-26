#!/bin/bash -e
#SBATCH --job-name=faststructure
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --output=faststructure-%j.out
#SBATCH --error=faststructure-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

##modules
module purge
##module load PLINK/2.00a2.3
#module load BCFtools/1.19-GCC-11.3.0
#module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1 

#filter out multiallelic variants
#bcftools view -m2 -M2 -v snps 2_filtered_genotyped.vcf -o biallelic.vcf

# prep plink files
#bcftools query -f '%CHROM\t%POS\t%END\n' biallelic.vcf > faw.bed


##modules
module load NeSI
module load fastStructure/1.0-gimkl-2020a-Python-2.7.18

##run faststructure
python /opt/nesi/CS400_centos7_bdw/fastStructure/1.0-gimkl-2020a-Python-2.7.18/bin/structure.py -K 10 --input=faw_MAF_geno_LD --output=faw