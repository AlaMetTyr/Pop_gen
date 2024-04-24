#!/bin/bash -e
#SBATCH --job-name=RG_data_add
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --output=RG-%j.out
#SBATCH --error=RG-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

# Load required modules
module load GATK/4.5.0.0-gimkl-2022a 
module load SAMtools/1.19-GCC-12.3.0

# Define your input and output folders
input_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/"  # Adjust to your SAM files location
output_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/corrected_bam"  # Location for output files

# Directory containing your SAM files
SAM_DIR="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/"  # Adjust to your path with SAM files

# Directory for output SAM files with read groups
OUTPUT_DIR="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/corrected_bam"  # Adjust to your output folder

# Make sure the output directory exists
mkdir -p $OUTPUT_DIR

# Loop through all SAM files in the directory
for FILE in $SAM_DIR/*.sam; do
    # Extract the base name without path and extension
    BASE_NAME=$(basename "$FILE" .sam)

    # Extract sample name before the first underscore
    SAMPLE_NAME=$(echo "$BASE_NAME" | cut -d'_' -f1)

    # Add read group with the sample name derived from the file name
    samtools addreplacerg \
      -r "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
      -o "${OUTPUT_DIR}/${BASE_NAME}_aln_rg.sam" \
      "$FILE"
    
    # Echo to SLURM output whether the RG line was added
    echo "Added Read Group to ${FILE} with RGID and RGSM: ${SAMPLE_NAME}. Output file: ${OUTPUT_DIR}/${BASE_NAME}_aln_rg.sam"
done

