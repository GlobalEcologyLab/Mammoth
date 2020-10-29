##### Setup #####
library(raster)
library(sf)
library(tidyverse)
library(cowplot)
library(paletteer)
library(extrafont)
library(zoo)

##### Create sampling effort map #####
megafauna <- read_csv("Datasets/sampling sites.csv")
grid <- raster("Datasets/mammoth_mask.grd")
projCRS <- "+proj=aea +lat_1=30 +lat_2=60 +lon_0=84.5 +ellps=WGS84"
newExt <- extent(-7078550, 6493945, 4751297, 10309582)
effort <- rasterize(megafauna[,1:2], grid, fun = "count") %>% projectRaster(crs = projCRS) %>% 
  crop(y = newExt)
holocene_cells <- raster("Datasets/mammoth_holocene_cells.grd") %>% projectRaster(crs = projCRS) %>% 
  crop(y = newExt)
effort_df <- as.data.frame(effort, xy = T)
holocene_df <- as.data.frame(holocene_cells, xy = T) %>% filter(!is.na(ensembleWeighted_extirpPattern_tol0.01_wICE))

poly_clip <- st_read("Datasets/polygon_extent_clip.shp")

cont_lines <- st_read("Datasets/ne_110m_coastline.shp") %>%
  st_transform(crs = projCRS) %>%
  st_crop(., newExt)

ocean <- st_read("Datasets/ne_110m_ocean.shp") %>%
  st_difference(., poly_clip) %>%
  st_transform(crs = projCRS) %>%
  st_buffer(., dist = 0) %>%
  st_crop(., newExt)

ggplot() +
  geom_sf(data = ocean, colour = NA, fill = "gray80") +
  geom_tile(data = holocene_df, mapping = aes(x,y), fill = "#C994C7FF") +
  geom_tile(data = effort_df, mapping = aes(x, y, fill = layer)) +
  geom_sf(data = cont_lines, fill = NA, colour = "black") +
  guides(fill = guide_colorbar(title = "# of Sites")) + 
  scale_fill_paletteer_c("scico::imola", na.value = "transparent", trans = "log10", guide = "colourbar") +
  coord_sf(expand = FALSE) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, colour = "black", face = "bold") ,
    axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5, colour = "black"),
    axis.text.x = element_text(size = 14, hjust = 0.5, colour = "black"),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.ticks = element_line(colour = "black", size = 0.5),
    panel.ontop = TRUE,
    panel.grid.major = element_line(color = "gray40", linetype = 3),
    legend.background = element_blank(),
    legend.position = c(0.93, 0.38)) +
  labs(x = NULL, y = NULL)

ggsave("Visualizations/sampling effort map.png", width = 10, height = 4, unit = "in")
projectRaster(effort, crs = "+proj=longlat +datum=WGS84 +no_defs") %>% as.data.frame(xy = T) %>%
group_by(y) %>% summarize(Sites = sum(layer, na.rm = T)) %>%
ggplot(aes(y, Sites)) + geom_bar(stat="identity") +
scale_x_binned(position = "top", limits = c(15,90)) +
xlab("Latitude") + coord_flip() + theme_minimal()
ggsave("Visualizations/sampling effort histogram.png", width = 2, height = 4, unit = "in")
