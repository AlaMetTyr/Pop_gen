library(qvalue)
library(dplyr)
library(ggplot2)

setwd("/nesi/nobackup/ga03488/Amy/FAW/assemblies/bwa/corrected_sam/GATK_ready/all_faw/")
FAW.C2=read.table("FAW_HiC8_summary_contrast.out",h=T)
hist(10**(-1*FAW.C2$log10.1.pval.),freq=F,breaks=40)
abline(h=1)

scaffolds <- read.table("scaffold_locations8.txt")  ## Extracted scaffold names from VCF
FAW.C2 <- as.data.frame(cbind(FAW.C2, scaffolds))

## Make a simple Manhattan plot to check the distribution of outlier SNPs
plot(FAW.C2$log10.1.pval.)
abline(h = 3, col="red")  #0.001 p-value threshold_

# Convert log p-values to regular p-values
pvalues <- as.data.frame(10^-FAW.C2[,6])
# Check p-values summary
summary(pvalues)
FAW.C2 <- cbind(pvalues,FAW.C2)  ##append it to the df
pvalues <- as.vector(FAW.C2$`10^-FAW.C2[, 6]`)
# Compute q-values
qval <- qvalue(p = pvalues)
# Check q-values summary
summary(qval$qvalues)
plot(qval$qvalues)
# Convert q-values to a dataframe
qvalues = as.data.frame(qval$qvalues)

# Append q-values to FAW.C2
FAW.C2 <- cbind(FAW.C2, qvalues)
FAW.C2 = rename(FAW.C2, qval=`qval$qvalues`, LocusName = V1)

# Extract SNPs with q-value < 0.01
selected_SNPs <- FAW.C2[FAW.C2$qval < 0.01, ]

# Save selected SNPs to file
write.table(selected_SNPs, "FAW_C2SNPs_qval01_01_scaf8.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################
##fixed code?###
################
# Read the scaffold_list.txt file into R
scaffold_list <- readLines("scaffold_locations8.txt")

# Split each line into parts and keep only the first three
scaffold_list_reformatted <- sapply(scaffold_list, function(x) {
  parts <- strsplit(x, "\\s+")[[1]]  # Splitting by whitespace (space, tab, etc.)
  paste(parts[1:3], collapse = " ")  # Join the first three parts with space
})

scaffold_list_reformatted <- sapply(scaffold_list, function(x) {
  parts <- strsplit(x, "\\s+")[[1]]  # Splitting by whitespace (space, tab, etc.)
  paste(parts[1:3], collapse = " ")  # Join the first three parts with space
})

FAW.C2_2 <- cbind(FAW.C2, scaffold_list_reformatted)

ggplot(FAW.C2, aes(x=MRK, y = log10.1.pval., color = V1)) + 
  geom_point(show.legend = FALSE, alpha = 1, size = 2) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  geom_point(selected_SNPs, mapping = aes(x=MRK, y = log10.1.pval.), color = "#EB4511", size = 2) +
  theme_classic()

