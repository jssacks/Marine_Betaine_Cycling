

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
#library(marmap)
library(ggrepel)
library(maps)
library(mapdata)
library(ggpubr)


#define inputs:
exp.loc.file <- "Meta_Data/Ingalls_Lab_Data/Kin_Exp_MetaData.csv"
region.pal <- c("#014d64", "#01a2d9", "#00887d", "#6794a7")

  
#load in experiment data:
exp.loc.dat <- read_csv(exp.loc.file) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("Gyre", "Equatorial", "Coastal", "NPTZ"))) %>%
  mutate(Compound = str_replace(Compound, "homarine", "Homarine")) %>%
  mutate(Long = -1*abs(Long))

big.map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-170, -110), ylim = c(-6, 53)) +
  geom_point(data = exp.loc.dat, 
  aes(x = Long, y = Lat, fill = Region, shape = Compound), 
  alpha = 0.8, 
  size = 5) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = region.pal) +
#  geom_jitter(data = dat.3, aes(Long, Lat, color = Local_Date), size = 4, alpha = 0.3) + 
#  scale_color_viridis() +
  theme_bw() +
  geom_label_repel(data = exp.loc.dat, 
                aes(x = Long, y = Lat, label = KinExp_ID), 
                box.padding = 1.2, min.segment.length = 0, 
                ylim = c(-2, 48),
                alpha = 0.9,
                size = 3, 
                segment.size = 0.2, 
                segment.curvature = -1e-20) +
  geom_rect(xmin = -120, xmax = -125,
            ymin = 46, ymax = 51, alpha = 0, color = "black") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  guides(shape = "none")
big.map



####Make map of Puget Sound Stations:

#get map data:
d <- map_data("worldHires", c("Canada", "usa", "Mexico"))
  
#get just PS experiment locatin data:
exp.loc.dat.PS <- exp.loc.dat %>%
  filter(Region == "Coastal")


#Plot
PS.map <- ggplot() + 
  geom_polygon(data = d, aes(x=long, y = lat, group = group),
                        fill = "gray90", color = "black", linewidth = 0.2) +
  theme_bw()+
  coord_map("conic", lat0 = 18, xlim=c(-123.5, -122), ylim=c(47.5,49)) +
  geom_jitter(data = exp.loc.dat.PS, 
              aes(x = Long, y = Lat, shape = Compound), 
              alpha = 0.9, 
              size = 5,
              fill = "#00887d", width = 0.03, height = 0.03) +
  scale_shape_manual(values = c(21, 22)) +
  geom_label_repel(data = exp.loc.dat.PS, 
                  aes(x = Long, y = Lat, label = KinExp_ID), 
                  box.padding = 1.5, min.segment.length = 0.2, 
                  alpha = 0.8,
                  color = "black",
                 # point.padding = 0.5,
                 # ylim = c(-2, 48),
                  size = 3, 
                  segment.size = 0.2, 
                  segment.curvature = -1e-20) +
  theme_test() +
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(override.aes = list(fill = "gray20")))
  
PS.map


map.plot <- ggarrange(NA, NA, NA, NA, NA,
                      NA, big.map, NA, PS.map, NA,
                      NA, NA, NA, NA, NA, 
                      nrow = 3, ncol = 5,
                      heights = c(0.03, 0.94, 0.03),
                      widths = c(0.02, 0.56, 0.02, 0.38, 0.02),
                      labels = c(NA, NA, NA, NA, NA, NA, "A", NA, "B", NA, NA, NA, NA, NA, NA))
map.plot

ggsave(map.plot, filename = "Figures/Flux_Paper_Figures/Maps_Plot.png",
       scale = 1.2,
       height = 6, width = 8, units = "in", bg = "white")


















































































































