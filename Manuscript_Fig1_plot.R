### Script to create sampling map for Masters ###
### James Reeve - Götebogs universitet
### 2021-03-19

### Prepare the envrionment ####
rm(list = ls())
dev.off() 
options(stringsAsFactors = FALSE) 
setwd("/Users/james/Documents/Courses/RGIS_course_2021-03-12/practicals")

### Packages ####
# Spatial analysis tools
library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(cleangeo)
library(sp)    ### I can probably replace this with an sf function !!!
library(ggnewscale)
library(sdmpredictors)
library(leaflet)
library(ggplot2)

### Upload data for base map ####
### Vector layers - these datasets
# Map of States & Provinces
tdwg <- readOGR("https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level4.geojson")
# Map of countries
countries <- readOGR("data/World_Countries/TM_WORLD_BORDERS-0.3.shp")
# Lakes polygons
lakes <- readOGR("data/World_Lakes/World_Lakes.shp")

### Raster layers
# Terrestrial temperature range
wclim <- raster::getData("worldclim", var = "bio", res = 5)
# Sea surface temperature range
Sea_temp_range <- load_layers("BO21_temprange_ss") # load_layers comes from the sdmpredictors package!

### Create dataset of sampling points ####
samples <- data.frame("lat" = c(66.475757, 63.747072, 
                                54.051279, 
                                61.596173, 43.342742, 
                                58.350077, 48.654292),
                      "long" = c(-70.498218, -68.489321, 
                                 -110.896082, 
                                 -149.345746, -124.350450, 
                                 -134.603122, -123.393027),
                      "pop_ID" = c("NsNUn", "NsNUd", 
                                   "NsABm &\n NsABk", 
                                   "TsAK", "TsOR", 
                                   "TuAK", "TuBC"),
                      "species" = c("Pungitius pungitius", "Pungitius pungitius", 
                                    "Pungitius pungitius",
                                    "Gasterosteus aculeatus", "Gasterosteus aculeatus",
                                    "Aulorhynchus flavidus", "Aulorhynchus flavidus"))

### 1. Cut base map to just North America ####
### First check that all datasets have this CRS:
### "+proj=longlat +datum=WGS84 +no_defs" using `crs(data)`
CoOrd <- c(-170,-12,30,90) # Note order: c(long_min, long_max, lat_min, lat_max)!
# Crop the provinces map
provinces_NA <- crop(tdwg, CoOrd)
# Crop lakes
lakes_NA <- crop(lakes, CoOrd)
# Crop countries
# But first I need to fix a topology error (i.e. part of a borderline loops back on itself)
# These issue can be detected and fixed automaticall with cleangeo
tmp <- clgeo_Clean(countries, strategy = "BUFFER")
countries_NA <- crop(tmp, CoOrd)
rm(tmp)
# Crop the rasters
wclim_NA <- crop(wclim, countries_NA)
Sea_temp_range_NA <- crop(Sea_temp_range, countries_NA)

### 2. Convert map projection to Azimuthal equal distance ####
### The standard map projection doesn't really give a good view of area and distance
### in Canada. To get a better map, I am converting all shape filed to Azimuthal equal
### distance. The nice thing about spatial packages in R is that I only have to specify
### this for one map, the others can use the CRS [Coordinate Reference System] of the 
### first map!

### Convert the CRS of the provinces map
# Find longitude of centre of Canada
cnt_NAm <- gCentroid(provinces_NA[provinces_NA$ISO_Code == "CA",])
cnt_NAm@coords[1] # Centorid Longitude = -98.35648
# Transform reference system - "+proj=aeqd" is proj4 code for Azimuthal equal distance projection
provinces_NA_aeqd <- spTransform(provinces_NA, 
                                 CRS = paste0("+proj=aeqd +lat_0=90 +lon_0=",cnt_NAm@coords[1]))

### Convert the rest of the data
countries_NA_aeqd <- spTransform(countries_NA, CRS = crs(provinces_NA_aeqd))
lakes_NA_aeqd <- spTransform(lakes_NA, CRS = crs(provinces_NA_aeqd))
# Rasters have a different function; projectRaster()
wclim_NA_aeqd <- projectRaster(wclim_NA, crs = crs(provinces_NA_aeqd))
Sea_temp_range_NA_aeqd <- projectRaster(Sea_temp_range_NA, crs = crs(provinces_NA_aeqd))

### Convert the sampling points
# First I need to change it to a spatial data type
coordinates(samples) <- ~ long + lat # Note: longitude is always 1st
# Assing a CRS value
crs(samples) <- crs(countries_NA)
# Now we can transform this CRS to azimuthal equal distance
samples_aeqd <- spTransform(samples, CRS = crs(provinces_NA_aeqd))

### 3. Plotting ####
### Note ggplot doesn't work with raster object, so as a quick workaround I converted
### them to a classical dataframe. I need to explore some packages that allow ggplot to
### plot rasters.

# Convert land temperature range
land_tmp <- as(wclim_NA_aeqd, "SpatialPixelsDataFrame")
land_tmp <- data.frame("long" = land_tmp@coords[,"x"],
                       "lat" = land_tmp@coords[,"y"],
                       "tmp_var" = land_tmp$bio2 / 10) 
# Note: it's common for GIS data to multiple temperature by 10 to convert floats to 
# integers as a way of saving space

# Convert sea surface temperature range
sea_tmp <- as(Sea_temp_range_NA_aeqd, "SpatialPixelsDataFrame")
sea_tmp <- data.frame("long" = sea_tmp@coords[,"x"],
                      "lat" = sea_tmp@coords[,"y"],
                      "tmp_var" = sea_tmp$BO21_temprange_ss)
# Make the ggplot
p1 <- ggplot()+
  #Sea surface temperature range(°C)
  geom_tile(data = sea_tmp, aes(x = long, y = lat, fill = tmp_var), alpha = 0.6)+ 
  scale_fill_gradientn(colours = topo.colors(10), breaks = c(0, 5, 10, 15, 20, 25))+
  labs(x = "Longitude", y = "Latitude", fill = "Temperature\nRange °C")+
  new_scale_fill()+ # This resets the gradient scale, to allow me to plot two different colour gradients!
  # Land temperature range (°C)
  geom_tile(data = land_tmp, aes(x = long, y = lat, fill = tmp_var), show.legend = FALSE)+ 
  scale_fill_gradientn(colours = grey.colors(10))+ 
  # Map of North American states/provinces/territories
  geom_sf(data = st_as_sf(provinces_NA_aeqd), colour = "grey", fill = NA)+ 
  # Map of North American countries
  geom_sf(data = st_as_sf(countries_NA_aeqd), colour = "black", fill = NA)+ 
  # Lakes
  geom_sf(data = st_as_sf(lakes_NA_aeqd), colour = "black", fill = "skyblue")+ 
  # Sampling points
  geom_sf(data = st_as_sf(samples_aeqd), aes(colour = species, shape = species),
          size = 3, show.legend = FALSE)+ 
  # Labels for sampling points
  geom_sf_label(data = st_as_sf(samples_aeqd), aes(label = pop_ID, colour = species), 
                size = 2.5, position = position_nudge(x = c(rep(6.8e5, 6), 7.8e5), 
                                                      y = c(rep(c(1.9e5, -1.9e5), 3), 0)),
                show.legend = FALSE)+ 
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(vjust = 1),
        legend.key.width = unit(0.05, "npc"), legend.key.height = unit(0.02, "npc"),)

### 4. Add pictures of fishes on right ####
library(jpeg)
library(grid)
library(gridExtra)

# Threespine picture
Ts_img <- readJPEG("/Users/james/Dropbox (Personal)/Master's/Species_photos/Ts_stickleback_Hazel_Cameron_Inglis_2018.jpg", TRUE)
Ts_grob <- rasterGrob(Ts_img, interpolate = TRUE)
Ts_grob$width <- unit(0.98, "npc")
Ts_grob$height <- unit(0.9, "npc")

# Ninespine picture
Ns_img <- readJPEG("/Users/james/Dropbox (Personal)/Master's/Species_photos/Ns_stickleback_Piet_Spanns_2006.jpg", TRUE)
Ns_grob <- rasterGrob(Ns_img, interpolate = TRUE)
Ns_grob$width <- unit(0.98, "npc")
Ns_grob$height <- unit(0.9, "npc")

# Tubesnout picture
Tu_img <- readJPEG("/Users/james/Dropbox (Personal)/Master's/Species_photos/Tubesnout_Hazel_Cameron_Inglis_2018.jpg", TRUE)
Tu_grob <- rasterGrob(Tu_img, interpolate = TRUE)
Tu_grob$width <- unit(0.98, "npc")
Tu_grob$height <- unit(0.9, "npc")

# Plot all three images
jpeg("/Users/james/Dropbox (Personal)/Master's/Species_photos/Multi-species_panel.jpg",
     height = 560, width = 300)

grid.arrange(grobs = list(Ts_grob, Ns_grob, Tu_grob), ncol = 1)

# Add coloured borders
grid.rect(y = 0.83, width = 0.98, height = 0.3, gp = gpar(lwd = 12, col = "green", fill = NA))
grid.rect(y = 0.5, width = 0.98, height = 0.3, gp = gpar(lwd = 12, col = "blue", fill = NA))
grid.rect(y = 0.17, width = 0.98, height = 0.3, gp = gpar(lwd = 12, col = "orange", fill = NA))

# Add plot symbols
grid.points(x = rep(0.9, 3), y = c(0.93, 0.58, 0.25), pch = c(17,15,16), 
            size = unit(c(0.1, 0.13, 0.11), "npc"), gp = gpar(col = c("green", "blue", "orange")),
            default.units = "npc")
dev.off()

# Load jpeg back into R
img <- readJPEG("/Users/james/Dropbox (Personal)/Master's/Species_photos/Multi-species_panel.jpg", TRUE)
p2 <- rasterGrob(img, interpolate = TRUE)

# Combine with map (p1)
Fig1 <- arrangeGrob(p1, p2, widths = c(4,1), layout_matrix = cbind(1,c(2,NA)), heights = c(1, 0.25))
grid.draw(Fig1)

png("/Users/james/Dropbox (Personal)/Comp_geno_2020/Figure1.png", width = 21, height = 14,
    units = "cm", res = 300)
grid.draw(Fig1)
dev.off()

