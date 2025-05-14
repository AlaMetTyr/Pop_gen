ibrary(ggplot2)
library(ggrepel)

# Filter samples from New Zealand and Australia
vs$label <- ifelse(vs$pop %in% c("NewZealand", "Australia"), vs$sample_ID, NA)

p <- ggplot(vs, aes(x=V2, y=V3, col=pop, shape=status)) +
  geom_point(size=4) +
  geom_text_repel(aes(label=label), size=5, 
                  box.padding = 0.8,  # Increase padding around labels
                  point.padding = 0.5, # Space labels from points
                  max.overlaps = Inf,  # Allow all labels to be shown
                  force = 2,           # Increase repulsion strength
                  min.segment.length = 0) +  # Always draw segment lines
  xlab(paste("PC1, ", format(p$V1[1] * 100 / sum(p$V1), digits=4), "%", sep='')) +
  ylab(paste("PC2, ", format(p$V1[2] * 100 / sum(p$V1), digits=4), "%", sep='')) + 
  theme_bw() +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=20),
        legend.title = element_blank())

print(p)

# Save as SVG
svg("my_plot.svg", height=5, width=7.2)
print(p)
dev.off()

#######################plot 3d pca ###################
#####################################################

library(ggplot2)
library(ggrepel)
library(plotly)

# Filter samples from New Zealand and Australia for labeling
vs$label <- ifelse(vs$pop %in% c("NewZealand", "Australia"), vs$sample_ID, NA)

# Create an interactive 3D PCA plot
p3d <- plot_ly(vs, 
               x = ~V2, y = ~V3, z = ~V4,  # PC1, PC2, PC3
               color = ~pop, symbol = ~status, 
               text = ~sample_ID,  # Hover labels show sample ID
               type = "scatter3d", mode = "markers") %>%
  layout(title = "3D PCA Plot",
         scene = list(xaxis = list(title = paste("PC1, ", format(p$V1[1] * 100 / sum(p$V1), digits=4), "%")),
                      yaxis = list(title = paste("PC2, ", format(p$V1[2] * 100 / sum(p$V1), digits=4), "%")),
                      zaxis = list(title = paste("PC3, ", format(p$V1[3] * 100 / sum(p$V1), digits=4), "%"))))

# Show the plot
p3d
