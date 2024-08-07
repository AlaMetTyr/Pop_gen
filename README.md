# Population genomics - from raw data to analyses  
Various codes for different stages of some population genomic work.
Personal use scripts therefore slightly messy.

## Script functions
### Generic scripts
`00-recursive_rename` to rename the raw data folders to a more manageable sample format including sample name and lane data  
  
### Genome assembly
`01-batch-bwa` Batch BWA mapping to the reference genome of Illumina Novaseq data  
`01.5-*-walltimegap` A second format for BWA script that accounts for walltime error. Recursively searches for output directories and only progresses with alignment for those that have not been processed.  
`00-check_RG` BWA does not automatically add the @RG sample information required for haplotype calling. Adds this based on the file name retroactively (BWA requires addition of -R variable)  
`01-GATK_preprocess`Merge, sort and index the sam files for each run of each sample
`0.2-GATK_Haplotyper` Prepping files for GATK (sort and index), Running GATK haplotype caller to generate gvcf files for each sample  
