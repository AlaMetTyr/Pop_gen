library(dplyr)
library(tidyr)
library(ggplot2)

samplelist <- read.table("./pop_faw.txt",header=TRUE)
all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 2:5){
  data <- read_delim(paste0("faw.",k,".meanQ"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$sample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data
# Remove rows with any NA values in all columns
all_data <- all_data %>%
  drop_na()

# Alternatively, to remove NAs in specific columns (e.g., "value" or "Q"):
# all_data <- all_data %>%
#   drop_na(value, Q)

# Merge the 'state' information from samplelist to all_data
all_data <- left_join(all_data, samplelist, by = "sample")

# Order samples by state
all_data$sample <- factor(all_data$sample, levels = unique(all_data$sample[order(all_data$state)]))

# Now plot the data without NAs
all_data %>%
  ggplot(aes(x = sample, y = value, fill = factor(Q))) + 
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette = "Set1", name = "K", labels = seq(1:5)) +
  facet_wrap(~k, ncol = 1)
