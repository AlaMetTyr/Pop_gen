#!/bin/bash -e
#SBATCH --job-name=RG_data_add
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --output=RG-%j.out
#SBATCH --error=RG-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

####################################################################################################################################################
##### if ran BWA without the -R option to include the @RG group required by GATK4.0 for haplotypecaller we will add this manually with samtools#####
####################################################################################################################################################

# Load required modules
module load GATK/4.5.0.0-gimkl-2022a 
module load SAMtools/1.19-GCC-12.3.0

# Directory containing your SAM files and output location
SAM_DIR="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/" 
OUTPUT_DIR="/nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/corrected_bam" 

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
    
    # Echo to SLURM output whether the RG line was added and what the sample name has been saved as to check for any errors
    echo "Added Read Group to ${FILE} with RGID and RGSM: ${SAMPLE_NAME}. Output file: ${OUTPUT_DIR}/${BASE_NAME}_aln_rg.sam"
done

