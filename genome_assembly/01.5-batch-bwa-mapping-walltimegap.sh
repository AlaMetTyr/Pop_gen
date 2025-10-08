#!/bin/bash -e


###########################
# Setup SLURM Environment #
###########################
#SBATCH --job-name=bwa-FAW-batch
#SBATCH --cpus-per-task=12
#SBATCH --time=142:00:00
#SBATCH --mem=20GB
#SBATCH --output=bwa-%j.out
#SBATCH --error=bwa-%j.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488
#################################################################################################################################
#####When walltime expired to then go back through and finish bwa alignments for the ones that werent completed previously#######
#######################Script checks for output directory and skips over those where this exists#################################
################################Also now includes the @RG fline (th -R flag of bwa)##############################################

### Executables####
module load BWA/0.7.17-GCC-11.3.0 




### Variables ###
workingdir="/nesi/project/ga03488/Amy/FAW/2025"
scratchdir="/nesi/nobackup/ga03488/Amy/FAW/2025"
ref_dir="/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference"
# Define the directory containing sample folders
sample_dir="/nesi/nobackup/ga03488/Amy/FAW/2025/raw"
# Flowcell ID for all samples (edit if needed)
FLOWCELL="EXT051"

## Iterate over each sample folder##
find "$sample_dir" -mindepth 1 -maxdepth 1 -type d | while IFS= read -r sample_folder; do
    # Extract the sample name from the folder path
    sample_name=$(basename "$sample_folder")

    # Check if the output directory already exists for the current sample
    out_dir="$scratchdir/bwa/$sample_name"
    if [ -d "$out_dir" ]; then
        echo "Output directory already exists for sample $sample_name. Skipping..."
        continue
    fi

    # Find the input FASTQ files for the current sample recursively
    fastq_r1=$(find "$sample_folder" -type f -name "*R1_001.fastq.gz" | head -n 1)
    fastq_r2=$(find "$sample_folder" -type f -name "*R2_001.fastq.gz" | head -n 1)

    # Check if both R1 and R2 files exist
    if [ ! -f "$fastq_r1" ] || [ ! -f "$fastq_r2" ]; then
        echo "Error: Missing input files for sample $sample_name"
        continue
    fi

    # Create the output directory for the current sample
    mkdir -p "$out_dir"

    # Change to the output directory
    cd "$out_dir" || exit

    # Try to extract lane & barcode info from filename if present
    LANE=$(basename "$fastq_r1" | grep -o 'L[0-9]\{3\}' | head -n1)
    SAMPLE_BARCODE=$(basename "$fastq_r1" | grep -o 'S[0-9]\+' | head -n1)

    # Fallbacks if not found
    LANE=${LANE:-L001}
    SAMPLE_BARCODE=${SAMPLE_BARCODE:-S00}

    RGID="${FLOWCELL}_${LANE}"
    RGSM="${sample_name}"
    PU="${FLOWCELL}.${LANE}.${SAMPLE_BARCODE}"

    # Add read group directly during alignment
    bwa mem -t 8 \
        -R "@RG\tID:${RGID}\tSM:${RGSM}\tPL:ILLUMINA\tLB:lib1\tPU:${PU}" \
        "$ref_dir/sfC.ver7.fa" "$fastq_r1" "$fastq_r2" \
        > "${sample_name}_aln.sam"

    echo " Completed ${sample_name} with RGID=${RGID}, RGSM=${RGSM}, RGPU=${PU}"
done
