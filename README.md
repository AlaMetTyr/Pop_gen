# Pop_gen
Step by step of data processing and analysis for pop gen
This repo is a step by step data analysis from raw fastq files for population genomic data from GBS sequencing.

## Adapter trimming from cutadapt
    cutadapt -a ADAPTER1 -g ADAPTER1 -a ADAPTER2 -g ADAPTER2 -o output.fastq input.fastq  
Specifying 2x forward common adapters.  
This is uploaded in `cutadapt.sl`  

## stacks process radtags  
    process_radtags -p <input_file_1> -o <output_directory_1> -b <barcodes_file_1> -e <enzyme_name> -q <quality_threshold> -r -c -t <output_format>
    process_radtags -p <input_file_2> -o <output_directory_2> -b <barcodes_file_1> -e <enzyme_name> -q <quality_threshold> -r -c -t <output_format>

## Ustacks
    ustacks -t fastq -f <input_file> -o <output_directory> -i <individual_id> -m <max_num_diff> -M <min_stack_depth>
process each seperately in ustacks.  

### Merging ustacks output
    cp -R <output_directory_1>/stacks <merged_directory>/stacks
    cp -R <output_directory_2>/stacks/* <merged_directory>/stacks/

## Run cstacks  
    cstacks -n <num_threads> -s <merged_directory>/stacks -o <output_directory> -M <popmap_file>

## Run sstacks  
    sstacks -b <batch_id> -s <stacks_directory> -o <output_directory> -M <popmap_file>

## Run populations  
    populations -b <batch_id> -s <stacks_directory> -P <populations_options>
