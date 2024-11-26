#!/bin/bash -e
#SBATCH --job-name=Mark_duplicates
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=5GB
#SBATCH --output=GATK_MD-%j.out
#SBATCH --error=GATK_MD-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

## Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  

input_folder="/nesi/nobackup/ga03488/Amy/FAW/Raw_data/rg_corr/GATK_ready"
output_folder="/nesi/nobackup/ga03488/Amy/FAW/Raw_data/rg_corr/GATK_ready/marked_duplicates/"

mkdir -p "${output_folder}"

for sorted_bam in "${input_folder}"/*_sorted.bam; do
  base_filename=$(basename "${sorted_bam}" _sorted.bam)
  marked_bam_file="${output_folder}/${base_filename}_marked.bam"
  metrics_file="${output_folder}/${base_filename}_duplicate_metrics.txt"

  # Check if output file already exists
  if [ -f "$marked_bam_file" ]; then
    echo "Output file ${marked_bam_file} already exists. Skipping marking duplicates."
  else
    gatk MarkDuplicates \
      -I "${sorted_bam}" \
      -O "${marked_bam_file}" \
      -M "${metrics_file}"

    samtools index "${marked_bam_file}"

    echo "Marked duplicates for ${base_filename}. Output: ${marked_bam_file}, Metrics: ${metrics_file}."
  fi
done
