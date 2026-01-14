#!/bin/bash -e
#SBATCH --job-name=GATK-pipeline-array
#SBATCH --output=gatk%x_%j.out     # log file
#SBATCH --error=gatk%x_%j.err      # error log file
#SBATCH --account=ga03488    # your NeSI project code
#SBATCH --time=24:00:00         # maximum run time hh:mm:ss
#SBATCH --mem=60G              # maximum memory available to GATK
#SBATCH --cpus-per-task=8
#SBATCH --array=0-$(($(wc -l < samples.txt)-1))%20


##this rtook a long time to run on full dataset, if runnning again need to remember. also joint genotyping (2weeks)
####################################################################
TMPDIR=/nesi/nobackup/ga03488/GATK_tmp/
mkdir -p ${TMPDIR}

module purge
module load GATK/4.3.0.0-gimkl-2022a

export _JAVA_OPTIONS=-Djava.io.tmpdir=${TMPDIR} 
export TMPDIR
export TEMP=${TMPDIR}
export TMP=${TMPDIR}

##root files##
REF=/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa 
OUT_ROOT=/nesi/nobackup/ga03488/Amy/FAW/2025/assemblies
KNOWN_SITES_SNP=/nesi/nobackup/ga03488/Amy/FAW/2025/global   

#txt file for array
printf "%s\n" /nesi/nobackup/ga03488/Amy/FAW/2025/GATK_ready/merged_bam/*.bam > samples_bam.txt

#outputs##
BAM_DIR="${OUT_ROOT}/bam_sorted_markdup";   mkdir -p "$BAM_DIR"
BQSR_DIR="${OUT_ROOT}/bqsr";                mkdir -p "$BQSR_DIR"
GVCF_DIR="${OUT_ROOT}/gvcf_samples";        mkdir -p "$GVCF_DIR"
mkdir -p logs

##indecing##
[[ -s ${REF}.fai ]] || samtools faidx "$REF"
DICT=${REF%.fa}.dict; DICT=${DICT%.fasta}.dict
[[ -s $DICT ]] || gatk CreateSequenceDictionary -R "$REF"


INPUT_SAM=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" samples.txt)
base=$(basename "$INPUT_SAM" .sam)
SAMPLE=$(echo "$base" | cut -d'_' -f1)
LANE=$(echo "$base" | grep -oE 'L[0-9]{3}' || true)
TAG=${SAMPLE}${LANE:+_"$LANE"}

echo ">>> Processing ${TAG}"
echo "SAM: $INPUT_SAM"


WORK=${SLURM_TMPDIR:-/tmp}/${TAG}
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT


samtools view -bS "$INPUT_SAM" \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} --tmpdir "${TMPDIR}" -o "$WORK/${base}_sorted.bam" -
samtools index "$WORK/${base}_sorted.bam"

##GATK mark duplicates
gatk --java-options "-Xmx12g" MarkDuplicatesSpark \
  -I "$WORK/${base}_sorted.bam" \
  -O "$WORK/${base}_marked.bam" \
  --conf "spark.executor.cores=${SLURM_CPUS_PER_TASK}" \
  --create-output-bam-index true

##BQSR
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

##haplotypecaller
gatk --java-options "-Xmx12g" HaplotypeCaller \
  -R "$REF" -I "$BAM_FOR_HC" \
  -O "$WORK/${TAG}.g.vcf.gz" \
  -ERC GVCF \
  --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK}

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
