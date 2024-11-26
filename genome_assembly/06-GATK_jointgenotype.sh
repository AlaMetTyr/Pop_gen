#!/bin/bash -e
#SBATCH --job-name=gatk_genotype
#SBATCH --cpus-per-task=12
#SBATCH --time=60:00:00
#SBATCH --mem=4GB
#SBATCH --output=genotype-%j.out
#SBATCH --error=genotype-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

########################################################################################
######## Script to make indivudal gvcfs per sample from bwa alignments in GATK##########
########################################################################################

### Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  
module load BCFtools/1.19-GCC-11.3.0

# Set the input folder containing SAM files and the output folder for BAM/GVCF files
combined_gvcf="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/merged_data/NZ_FAW.g.vcf.gz" 
output_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/vcf"
reference_genome="/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa"

# Define the output for joint genotyping
genotyped_vcf="${output_folder}/2_genotyped.vcf.gz"

# Perform joint genotyping with the combined GVCF
echo "Performing joint genotyping on the combined GVCF"
gatk GenotypeGVCFs \
  -R "${reference_genome}" \
  -V "${combined_gvcf}" \
  -O "${genotyped_vcf}"

# Define the output for the filtered VCF
filtered_vcf="${output_folder}/2_filtered_genotyped.vcf.gz"

# Apply filtering with specific criteria
echo "Applying filtering to the genotyped VCF"
gatk VariantFiltration \
  -R "${reference_genome}" \
  -V "${genotyped_vcf}" \
  -O "${filtered_vcf}" \
  --filter-expression "QD < 2.0" --filter-name "LowQD" \
  --filter-expression "FS > 60.0" --filter-name "HighFS" \
  --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
  --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum"

# Re-index the filtered VCF
echo "Re-indexing the filtered VCF"
bcftools index "${filtered_vcf}"

echo "Script completed."