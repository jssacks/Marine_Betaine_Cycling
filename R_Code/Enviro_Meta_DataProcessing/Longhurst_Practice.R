#install.packages("sf")
#install.packages("rnaturalearth")


library(rnaturalearth)

library(sf)


longhurst <- sf::read_sf("Meta_Data/Longhurst_world_v4_2010.shp")
names(longhurst)

x <- st_drivers()

y <- read_sf("Meta_Data/Longhurst_world_v4_2010.shp")

z <- st_read("Meta_Data/Longhurst_Provinces/Longhurst_world_v4_2010.shp")



# World map
#library(rnaturalearth)
world_map <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))


# Base map
kk <- ggplot() +
  geom_sf(data = world_map, size = .2, fill = "gray80", col = "gray90") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))
kk

longhurst <- sf::read_sf("Meta_Data/Longhurst_Provinces/Longhurst_world_v4_2010.shp")
names(longhurst)

sf_use_s2(FALSE)

longhurst <- longhurst %>% 
  sf::st_simplify(dTolerance = 0.01) %>% 
  dplyr::group_by(ProvCode,ProvDescr) %>% 
  dplyr::summarise()
# plot(longhurst)

y <- longhurst %>%
  filter(ProvCode == "PEQD")

x <- as.data.table(y$geometry)

u <- x[[1]][[1]][[1]]

# plot
kk+  
  geom_sf(data = longhurst, aes(fill = ProvCode), size = .2, col = "grey50", alpha=.4)+
  ggtitle(paste("Longhurst Biogeochemical Provinces -", length(unique(longhurst$ProvCode)),"provinces"))+
  theme(legend.position="none")+
  geom_sf_text(data = longhurst %>% group_by(ProvDescr) %>% summarize(n()), aes(label = ProvDescr), colour = "grey20", check_overlap=TRUE)+
  coord_sf(expand = FALSE)

