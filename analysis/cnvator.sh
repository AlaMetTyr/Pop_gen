#!/bin/bash -e
#SBATCH --job-name=cnvnator
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH --output=cnvnator-%j.out
#SBATCH --error=cnvnator-%j.err
#SBATCH --account=ga03488

module load  CNVnator/0.4-GCC-7.4.0 


BAM="/home/a.vaughan/nobackup_ga03488/Amy/FAW/2025/bam_sorted_markdup/Tar-01_sorted_marked.bam"
REF="/nesi/project/ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa"
OUT="sample.cnvnator-Tar-01"
BIN=100   


cnvnator -root ${OUT}.root -tree ${BAM}
cnvnator -root ${OUT}.root -his ${BIN} -d $(dirname ${REF})
cnvnator -root ${OUT}.root -stat ${BIN}
cnvnator -root ${OUT}.root -partition ${BIN}
cnvnator -root ${OUT}.root -call ${BIN} > ${OUT}.calls.txt
