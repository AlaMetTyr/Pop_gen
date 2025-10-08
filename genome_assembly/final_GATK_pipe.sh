#!/bin/bash -e
#SBATCH --job-name=GATK-pipeline-array
#SBATCH --output=gatk%x_%j.out     # log file
#SBATCH --error=gatk%x_%j.err      # error log file
#SBATCH --account=ga03488    # your NeSI project code
#SBATCH --time=24:00:00         # maximum run time hh:mm:ss
#SBATCH --mem=60G              # maximum memory available to GATK
#SBATCH --cpus-per-task=8
#SBATCH --array=0-$(($(wc -l < samples.txt)-1))%20

####################################################################
## create temporary directory for Java so it does not fill up /tmp##
TMPDIR=/nesi/nobackup/ga03488/GATK_tmp/
mkdir -p ${TMPDIR}

# remove other modules that may be loaded
# load specific GATK version
module purge
module load GATK/4.3.0.0-gimkl-2022a

# tell Java to use ${TMPDIR} as the temporary directory
export _JAVA_OPTIONS=-Djava.io.tmpdir=${TMPDIR} 
####################################################################check this is done!####
#have to have made a text file of the sam files so can run as an array######################
#ls /nesi/nobackup/ga03488/Amy/FAW/assemblies/test_data/alignments_sam/corrected_bam/*.sam \
#  > samples.txt
############################################################################################

# 
export TMPDIR
export TEMP=${TMPDIR}
export TMP=${TMPDIR}

############################################
EDIT THESE PATHS
############################################
REF=/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa ##double check correct
OUT_ROOT=/nesi/nobackup/ga03488/Amy/FAW/2025/assemblies

# Global known-sites VCFs (same reference; bgzip+tabix indexed). Leave blank to skip BQSR.
KNOWN_SITES_SNP=        # e.g. /nesi/nobackup/ga03488/Amy/FAW/2025/global/merged.Invasive.g.vcf.gz
KNOWN_SITES_INDEL=      # 

############################################
# Output dirs
############################################
BAM_DIR="${OUT_ROOT}/bam_sorted_markdup";   mkdir -p "$BAM_DIR"
BQSR_DIR="${OUT_ROOT}/bqsr";                mkdir -p "$BQSR_DIR"
GVCF_DIR="${OUT_ROOT}/gvcf_samples";        mkdir -p "$GVCF_DIR"
mkdir -p logs

############################################
# Reference indexes
############################################
[[ -s ${REF}.fai ]] || samtools faidx "$REF"
DICT=${REF%.fa}.dict; DICT=${DICT%.fasta}.dict
[[ -s $DICT ]] || gatk CreateSequenceDictionary -R "$REF"

############################################
# Pick sample from array
############################################
INPUT_SAM=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" samples.txt)
base=$(basename "$INPUT_SAM" .sam)
SAMPLE=$(echo "$base" | cut -d'_' -f1)
LANE=$(echo "$base" | grep -oE 'L[0-9]{3}' || true)
TAG=${SAMPLE}${LANE:+_"$LANE"}

echo ">>> Processing ${TAG}"
echo "SAM: $INPUT_SAM"

############################################
# Per-task fast scratch
############################################
WORK=${SLURM_TMPDIR:-/tmp}/${TAG}
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

############################################
# 1) SAM -> sorted BAM (streaming)
############################################
samtools view -bS "$INPUT_SAM" \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} --tmpdir "${TMPDIR}" -o "$WORK/${base}_sorted.bam" -
samtools index "$WORK/${base}_sorted.bam"

############################################
# 2) Mark duplicates
############################################
gatk --java-options "-Xmx12g" MarkDuplicatesSpark \
  -I "$WORK/${base}_sorted.bam" \
  -O "$WORK/${base}_marked.bam" \
  --conf "spark.executor.cores=${SLURM_CPUS_PER_TASK}" \
  --create-output-bam-index true

############################################
# 3) BQSR (skip cleanly if no known sites)
############################################
KNOWN_ARGS=()
[[ -n "${KNOWN_SITES_SNP}"   ]] && KNOWN_ARGS+=( --known-sites "${KNOWN_SITES_SNP}" )
[[ -n "${KNOWN_SITES_INDEL}" ]] && KNOWN_ARGS+=( --known-sites "${KNOWN_SITES_INDEL}" )

BAM_FOR_HC="$WORK/${base}_marked.bam"
if (( ${#KNOWN_ARGS[@]} > 0 )); then
  gatk --java-options "-Xmx12g" BaseRecalibrator \
    -R "$REF" -I "$WORK/${base}_marked.bam" \
    "${KNOWN_ARGS[@]}" \
    -O "$WORK/${base}.recal.table"

  gatk --java-options "-Xmx12g" ApplyBQSR \
    -R "$REF" -I "$WORK/${base}_marked.bam" \
    --bqsr-recal-file "$WORK/${base}.recal.table" \
    -O "$WORK/${base}.recal.bam"

  samtools index "$WORK/${base}.recal.bam"
  BAM_FOR_HC="$WORK/${base}.recal.bam"
else
  echo "NOTE: No known sites provided; skipping BQSR."
fi

############################################
# 4) HaplotypeCaller (multi-threaded) → gVCF
############################################
gatk --java-options "-Xmx12g" HaplotypeCaller \
  -R "$REF" -I "$BAM_FOR_HC" \
  -O "$WORK/${TAG}.g.vcf.gz" \
  -ERC GVCF \
  --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK}

############################################
# 5) Move results out of scratch
############################################
mv -f "$WORK/"*_sorted.bam "$BAM_DIR"/
mv -f "$WORK/"*_sorted.bam.bai "$BAM_DIR"/
mv -f "$WORK/"*_marked.bam "$BAM_DIR"/
mv -f "$WORK/"*_marked.bam.bai "$BAM_DIR"/

if [[ -s "$WORK/${base}.recal.bam" ]]; then
  mv -f "$WORK/${base}.recal.bam" "$BAM_DIR"/
  mv -f "$WORK/${base}.recal.bam.bai" "$BAM_DIR"/
  mv -f "$WORK/${base}.recal.table" "$BQSR_DIR"/
fi

mv -f "$WORK/${TAG}.g.vcf.gz" "$GVCF_DIR"/
mv -f "$WORK/${TAG}.g.vcf.gz.tbi" "$GVCF_DIR"/ || true

echo "<<< Done ${TAG}"
