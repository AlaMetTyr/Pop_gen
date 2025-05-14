##calculated weighted mean fst from sliding window between populations

# Load necessary library
library(dplyr)

# Read the Fst data
df <- read.csv("subset-stats-popsize.Fst.Dxy.pi.csv", header=T)  # Replace with actual file name
colnames(df)
# Remove rows with missing Fst values for any of the specified columns
fst_columns <- c("Fst_Nelson.hetero_Nelson.homo")

df <- df %>% filter(complete.cases(df[, fst_columns]))

# Loop through each Fst column to compute weighted mean Fst
for (fst_col in fst_columns) {
  overall_fst <- sum(df$sites * df[[fst_col]], na.rm = TRUE) / sum(df$sites, na.rm = TRUE)
  cat(fst_col, "Overall Fst:", round(overall_fst, 6), "\n")
}


# Load necessary libraries
library(ggplot2)
library(reshape2)

fst_data <- read.table("clipboard", header = TRUE, sep = "\t")


# Create a matrix for the heatmap
fst_matrix <- dcast(fst_data, Pop1 ~ Pop2, value.var = "Fst")

# Convert to long format for ggplot
fst_long <- melt(data_matrix, id.vars = "pop1")

# Plot the heatmap
ggplot(fst_long, aes(x = Pop1, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10)) +
  labs(x = "Population", y = "Population", fill = "Fst", title = "Pairwise Fst Heatmap")

