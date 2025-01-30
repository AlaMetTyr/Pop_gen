library(ggplot2)
library(maps)

##insert tent data
coord <- read.table(file = "clipboard", 
                    sep = "\t", header=TRUE)

library(ggplot2)
library(rnaturalearth)
library(ggspatial)
# Load the world map from rnaturalearth
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot with a nicer, modern look
ggplot(data = world) +
  geom_sf(fill = "lightgray", color = "white") +  # Map base
  geom_point(data = coord, aes(x = long, y = lat), 
             color = "black", size = 3) +           # Your points
  annotation_scale(location = "bl", width_hint = 0.5) +  # Add scale bar
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +  # Add compass
  theme_minimal() + 
  labs(title = "Locations on a Modern World Map",
       x = "Longitude", y = "Latitude")


##nz
# Plot with New Zealand centered
ggplot(data = world) +
  geom_sf(fill = "lightgray", color = "white") +  # Map base
  geom_point(data = coord, aes(x = long, y = lat), 
             color = "black", size = 3) +           # Your points
  #annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +  # Add compass
  coord_sf(crs = "+proj=robin +lon_0=180") +  # Center the map around NZ
  theme_minimal() + 
  labs(title = "Locations on a Modern World Map (Centered on New Zealand)",
       x = "Longitude", y = "Latitude")
