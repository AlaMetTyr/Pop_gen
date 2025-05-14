#!/bin/bash -e
#SBATCH --job-name=gatk_test
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --output=GATK-table.out
#SBATCH --error=GATKtable.err
#SBATCH --nodes=1 #default
#SBATCH --account=ga03488

# Load required modules
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0 

########################################################################################
######## Script to make indivudal gvcfs per sample from bwa alignments in GATK##########
########################################################################################

### Executables
module load GATK/4.5.0.0-gimkl-2022a   
module load SAMtools/1.19-GCC-12.3.0  
module load BCFtools/1.19-GCC-11.3.0

# Set the input folder containing SAM files and the output folder for BAM/GVCF files
genotyped_vcf="/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/gvcf/FAW_geno-filt.vcf.gz" 
reference_genome="/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa"

gatk VariantsToTable \
  -R "${reference_genome}" \
  -V "${genotyped_vcf}" \
  -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
  -GF DP -GF AD -GF GQ \
  -O summary_table_output_varfilter.tsv

