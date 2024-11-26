#!/bin/bash -e
#SBATCH --job-name=extract_Ache
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --output=ache.out
#SBATCH --error=ache.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

########################################################################################
######## Script to make indivudal gvcfs per sample from bwa alignments in GATK##########
########################################################################################
### Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  

# Set the input folder containing SAM files and the output folder for BAM/GVCF files
input_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/gvcf/" ## This is the out put of 00-check-RG if didnt add the @RG in bwa
output_folder="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/merged_data/"


# Define the reference genome coordinates for XP_022819835
REGION="chr:start-end" # Replace this with actual output from faidx

# Loop through each BAM file and extract the specified region
for bamfile in /path/to/bamfiles/*.bam; do
    samtools view -b "$bamfile" "$REGION" > "${bamfile%.bam}_XP_022819835.bam"
done



for extracted_bam in /path/to/bamfiles/*_XP_022819835.bam; do
    samtools sort -o "${extracted_bam%.bam}_sorted.bam" "$extracted_bam"
    samtools index "${extracted_bam%.bam}_sorted.bam"
done
