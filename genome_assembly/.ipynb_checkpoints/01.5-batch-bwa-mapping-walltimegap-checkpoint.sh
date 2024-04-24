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

### Executables
module load BWA/0.7.17-GCC-11.3.0 

###Executables
module load BWA/0.7.17-GCC-11.3.0 

#############
# Variables #
#############
workingdir="/nesi/project/ga03488/Amy/"
scratchdir="/nesi/nobackup/ga03488/Amy/FAW/assemblies"
scripts="$workingdir/Scripts/genome_assembly/bwa"
logdir="$scripts/logs"
ref_dir="/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference"

# Define the directory containing sample folders
sample_dir="/nesi/nobackup/ga03488/Amy/FAW/"

# Iterate over each sample folder
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

    # Perform alignment using BWA
    bwa mem "$ref_dir/sfC.ver7.fa" "$fastq_r1" "$fastq_r2" > "${sample_name}_aln.sam"
done
