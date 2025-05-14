#!/bin/bash -e
#SBATCH --job-name=bayscan_reformat
#SBATCH --time=168:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=12
#SBATCH --output=reshaper.out
#SBATCH --error=reshaper.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

module purge
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load PLINK/1.09b6.16 

# Convert VCF to PLINK
#vcftools --vcf ./faw_bialSNP_MAF_genov2.vcf --plink --out faw_populations.plink

# Convert PLINK to BED
#plink --file faw_populations.plink --make-bed --allow-extra-chr --noweb --out faw_populations

module purge
module load Java/1.8.0_144

java -Xmx1024m -Xms512m -jar ./PGDSpider_3.0.0.0/PGDSpider3-cli.jar -inputfile ../faw_bialSNP_MAF_geno.vcf -inputformat VCF -outputfile faw_populations.pgd -outputformat PGD -spid VCF_PGD.spid

java -Xmx1024m -Xms512m -jar ./PGDSpider_3.0.0.0/PGDSpider3-cli.jar -inputfile faw_populations.pgd -inputformat PGD -outputfile faw_populations.bs -outputformat GESTE_BAYE_SCAN