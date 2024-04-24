#!/bin/bash -e
#SBATCH --job-name=gatk_test
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=4GB
#SBATCH --output=GATK-%j.out
#SBATCH --error=GATK-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

########################################################################################
######## Script to make indivudal gvcfs per sample from bwa alignments in GATK##########
########################################################################################
### Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  

# Set the input folder containing SAM files and the output folder for BAM/GVCF files
input_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/corrected_bam/" ## This is the out put of 00-check-RG if didnt add the @RG in bwa
output_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/gvcf_samples/"

# Create the output folder if it doesn't exist
mkdir -p "${output_folder}"

# Loop through all SAM files in the input folder
for sam_file in "${input_folder}"*.sam; do
  # Extract the base filename without extension
  filename=$(basename "${sam_file}" .sam)

  # Extract the sample name and lane information
  sample_name=$(echo "${filename}" | cut -d'_' -f1)
  lane_info=$(echo "${filename}" | grep -o 'L[0-9]{3}')

  # Construct the sample name with lane info
  sample_with_lane="${sample_name}_${lane_info}"
  
  # Convert SAM to BAM to make file more manageable to handle
  bam_file="${output_folder}${filename}.bam"
  samtools view -Sb "${sam_file}" > "${bam_file}"

  # Sort the BAM file
  sorted_bam_file="${output_folder}${filename}_sorted.bam"
  samtools sort -o "${sorted_bam_file}" "${bam_file}"

  # Index the sorted BAM file
  samtools index "${sorted_bam_file}"

  # Call HaplotypeCaller to generate GVCF
  gvcf_file="${output_folder}${sample_with_lane}.g.vcf.gz"
  gatk HaplotypeCaller \
    -R /nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa \ ## the same reference genome you used to assemble in BWA
    -I "${sorted_bam_file}" \
    -O "${gvcf_file}" \
    -ERC GVCF

done
