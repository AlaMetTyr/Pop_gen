library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(Biostrings)


# Read the GFF file
gff <- import("./OGS7.0_20190530.gff", format = "gff3")

# Extract gene names and descriptions from the metadata
gff$Gene <- as.character(mcols(gff)$Name)
gff$Description <- sapply(mcols(gff)$Note, function(x) ifelse(length(x) > 0, as.character(x), NA))
head(gff)

# Load the genome sequence (FASTA file)
genome_seq <- readDNAStringSet("~/project_ga03488/Amy/Scripts/genome_assembly/reference/sfC.ver7.fa", format = "fasta")

# Load SNP data
selected_SNPs <- read.table("selected_SNPs.txt", header = TRUE, stringsAsFactors = FALSE)

# Extract scaffold and position from LocusName
selected_SNPs$scaffold <- sub(":.*", "", selected_SNPs$V1)
selected_SNPs$position <- as.numeric(sub(".*:", "", selected_SNPs$V1))

# Check SNP data structure
head(selected_SNPs)

# Convert the GFF data to a GRanges object
gff_granges <- GRanges(
  seqnames = seqnames(gff),  # Scaffold names
  ranges = IRanges(start = start(gff), end = end(gff)),  # Start and end positions
  gene = gff$Gene,  # Gene name
  description = gff$Description  # Gene description
)

# Convert the SNPs to a GRanges object
snp_granges <- GRanges(
  seqnames = selected_SNPs$scaffold,  # Scaffold names
  ranges = IRanges(start = selected_SNPs$position, end = selected_SNPs$position)  # SNP positions
)

# Check the structure of the GRanges objects
head(gff_granges)
head(snp_granges)

# Find overlaps between SNPs and genes in the GFF file
overlaps <- findOverlaps(snp_granges, gff_granges)

# Extract relevant information from the GFF file for overlapping genes
gene_annotations <- data.frame(
  SNP = selected_SNPs$V1[queryHits(overlaps)],  # SNP identifiers
  Gene = mcols(gff_granges)$gene[subjectHits(overlaps)],  # Gene names
  Description = mcols(gff_granges)$description[subjectHits(overlaps)]  # Gene descriptions
)

# View the gene annotations
head(gene_annotations)

# Merge the gene annotations with the selected SNPs
annotated_SNPs <- merge(selected_SNPs, gene_annotations, by.x = "V1", by.y = "SNP", all.x = TRUE)

# View the annotated SNPs
head(annotated_SNPs)

# Save the annotated SNPs to a new file
write.table(annotated_SNPs, "annotated_SNPs-0.01qval-HiC8.txt", sep = "\t", quote = FALSE, row.names = FALSE)

)
