#!/bin/bash -e
#SBATCH --job-name=gatk_test
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --output=GATK-%j.out
#SBATCH --error=GATK-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

### Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  

# Define the directory containing GVCF files
gvcf_dir="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/gvcf/"
combined_gvcf="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/gvcf/combined_gvcf.vcf.gz"

# Define reference genome
reference_genome="/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa"

# Start combining GVCFs
echo "Checking for missing index files and combining GVCFs..."

# Initialize a list to hold the valid GVCFs
gvcf_list=""

# Loop through all GVCFs in the directory
for gvcf_file in "${gvcf_dir}"*.g.vcf.gz; do
  # Get the expected index file
  index_file="${gvcf_file}.tbi"

  # Check if the index file exists
  if [[ ! -f "${index_file}" ]]; then
    echo "Skipping ${gvcf_file} because the corresponding index file ${index_file} is missing."
    continue
  fi

  # Add this GVCF to the list if it has an index
  gvcf_list+=" -V ${gvcf_file}"
done

# Check if there are any valid GVCFs to combine
if [[ -z "${gvcf_list}" ]]; then
  echo "No valid GVCFs with index files found."
  exit 1  # Exit if no valid GVCFs are found
fi

# Combine the valid GVCFs into a single GVCF
gatk CombineGVCFs \
  -R "${reference_genome}" \
  ${gvcf_list} \
  -O "${combined_gvcf}"

echo "GVCFs combined successfully."

# Define the output for joint genotyping
genotyped_vcf="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/gvcf/test_genotyped.vcf.gz"

# Perform joint genotyping with the combined GVCF
echo "Performing joint genotyping on the combined GVCF"
gatk GenotypeGVCFs \
  -R "${reference_genome}" \
  -V "${combined_gvcf}" \
  -O "${genotyped_vcf}"

# Define the output for the filtered VCF
filtered_vcf="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/gvcf/filtered_merged_test.vcf.gz"

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

