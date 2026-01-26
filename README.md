# Population genomics - from raw data to analyses  
Various codes for different stages of some population genomic work.
Personal use scripts therefore slightly messy.

## Script functions
### Generic scripts
`00-recursive_rename` to rename the raw data folders to a more manageable sample format including sample name and lane data  
  
### Genome assembly
`01-batch-bwa` Batch BWA mapping to the reference genome of Illumina Novaseq data  
`01.5-*-walltimegap` A second format for BWA script that accounts for walltime error. Recursively searches for output directories and only progresses with alignment for those that have not been processed.*  
`00-check_RG` BWA does not automatically add the @RG sample information required for haplotype calling. Adds this based on the file name retroactively (BWA requires addition of -R variable)  
`01-GATK_preprocess`Merge, sort and index the sam files for each run of each sample  
`0.2-GATK_Haplotyper` Prepping files for GATK (sort and index), Running GATK haplotype caller to generate gvcf files for each sample  






* _have updated `01.5-*-walltimegap` to run bwa with the -r variable to add the @RG sample name (Oct 2025)
* _ have also updated ther GATK pipeline to be one code that runs as an array, including mark duplicates, BQSR and haplotypecaller


# Updated analyses (Jan 2026)  
### Genome assembly  
Updated to include BQSR in GATK pipeline  

### Analysis  
`phylo.sh` SNV phylogeny  
`relatedness-plink.sh` relatedness function for populations in plink2  
`vcftools-relatedness.sh` relatedness function for popultions in vcftools- compared with the plink2 function and not used further  
`tpi.sl` tpi extraction and analysis for constructing PCA for strain determination  
`fst-vcftools.sh` updated loop for weir-fst-pop between all populations in vcftools  
`faststructure.sh` updated loop for fastsructure  
`cnvator.sh` script to use CNVnator for copy number variant discovery. Used for detecting duplication in the AChE locus
