#!/bin/bash -e
#SBATCH --job-name=cutadapt_popgen # job name (shows up in the queue)
#SBATCH --time=00:01:00      # Walltime (HH:MM:SS)
#SBATCH --mem=512MB          # Memory in MB
#SBATCH --account=ga03488

input_files=("SQ2095_HKKCWDMXY_s_1_fastq.txt" "SQ2095_HKKCWDMXY_s_2_fastq.txt")  # List of input file names
for input_file in "${input_files[@]}"; do
    output_file="$(basename "$input_file" | sed 's/\.[^.]*$//').trimmed.fastq"
cutadapt -g AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CGAGATCGGAAGAGCGGACTTTAAGC -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGCA -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -o "$output_file" "$input_file"
done
