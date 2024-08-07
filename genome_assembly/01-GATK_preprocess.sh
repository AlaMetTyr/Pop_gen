#!/bin/bash -e
#SBATCH --job-name=gatk_preprocessing
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH --output=GATK_preprocess-%j.out
#SBATCH --error=GATK_preprocess-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

### Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  

# Set the input folder containing SAM files and the output folder for BAM files
input_folder="/nesi/nobackup/ga03488/Amy/FAW/Raw_data/bwa/corrected_sam/"
output_folder="/nesi/nobackup/ga03488/Amy/FAW/Raw_data/bwa/corrected_sam/GATK_ready"

# Create the output folder if it doesn't exist
echo "Creating output folder at ${output_folder}"
mkdir -p "${output_folder}"

# Step 1: Convert SAM to BAM, then sort and index
for sam_file in "${input_folder}"/*.sam; do
  echo "Processing SAM file: ${sam_file}"
  
  # Extract the base filename without extension
  filename=$(basename "${sam_file}" .sam)
  
  # Extract the sample name
  sample_name=$(echo "${filename}" | cut -d'_' -f1)
  
  # Skip processing if it's the 01-A sample
  if [[ "${sample_name}" == "01-A" ]]; then
    echo "Skipping sample 01-A, already processed."
    continue
  fi
  
  echo "Processing SAM file: ${sam_file}"

  # Convert SAM to BAM
  bam_file="${output_folder}/${filename}.bam"
  echo "Converting SAM to BAM: ${bam_file}"
  samtools view -Sb "${sam_file}" > "${bam_file}"

  # Sort the BAM file
  sorted_bam_file="${output_folder}/${filename}_sorted.bam"
  echo "Sorting BAM: ${sorted_bam_file}"
  samtools sort -o "${sorted_bam_file}" "${bam_file}"

  # Index the sorted BAM file
  echo "Indexing BAM: ${sorted_bam_file}"
  samtools index "${sorted_bam_file}"
done

# Step 2: Merge BAM files with the same sample name
echo "Merging BAM files for samples..."
declare -A sample_files

# Collect BAM files to merge by sample name (excluding lane information)
for sorted_bam_file in "${output_folder}"/*_sorted.bam; do
  # Extract the sample name (excluding lane information)
  base_filename=$(basename "${sorted_bam_file}" _sorted.bam)
  sample_name=$(echo "${base_filename}" | cut -d'_' -f1)

  # Append to the list of BAM files for this sample name
  sample_files["${sample_name}"]="${sample_files["${sample_name}"]} ${sorted_bam_file}"
done

# Merge BAM files for each sample and then sort and index
for sample_name in "${!sample_files[@]}"; do
  merged_bam_file="${output_folder}/${sample_name}.bam"
  
  # Echo which files are being merged for this sample
  echo "Merging BAM files for sample ${sample_name}: ${sample_files["${sample_name}"]}"
  
  # Merge the BAM files for the sample
  samtools merge -f "${merged_bam_file}" ${sample_files["${sample_name}"]}
  echo "Merged BAM file: ${merged_bam_file}"
  
  # Sort the merged BAM file
  sorted_merged_bam_file="${output_folder}/merged_bam/${sample_name}_sorted.bam"
  echo "Sorting merged BAM: ${sorted_merged_bam_file}"
  samtools sort -o "${sorted_merged_bam_file}" "${merged_bam_file}"

  # Index the sorted BAM file
   echo "Indexing merged BAM: ${sorted_merged_bam_file}"
  samtools index "${sorted_merged_bam_file}"

  echo "Merged, sorted, and indexed BAM file for sample ${sample_name}. Output: ${sorted_merged_bam_file}"
done

echo "Script completed."

