# Pop_gen
Step by step of data processing and analysis for pop gen
This repo is a step by step data analysis from raw fastq files for population genomic data from GBS sequencing.

## Adapter trimming from cutadapt
`cutadapt -a ADAPTER1 -g ADAPTER1 -a ADAPTER2 -g ADAPTER2 -o output.fastq input.fastq`  
Specifying 2x forward (5') adapters with the -g flag, and reverse (3') adapters with the -a flag.  
This is uploaded in `cutadapt.sl`
