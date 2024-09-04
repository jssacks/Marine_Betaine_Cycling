####

library(tidyverse)
library(lubridate)


##Enviro_Data_Curration: 

#Define inputs:

#TN397 (Gradeints 4)
g4.nut.uw.file <- "Meta_Data/Enviro_Meta_Data/TN397/Gradients4_TN397_Nutrients_UW.csv"
g4.pc.uw.file <- "Meta_Data/Enviro_Meta_Data/TN397/Gradients4_TN397_PPPCPN_UW.csv"
g4.nut.pc.discrete.file <- "Meta_Data/Enviro_Meta_Data/TN397/TN397_Gradients4_NutrientsAndParticulates.csv"
g4.ts.file <- "Meta_Data/Enviro_Meta_Data/TN397/TN397_Gradients4_uw_tsg.csv"
g4.bact.file <- "Meta_Data/Enviro_Meta_Data/TN397/Influx_Stations_Gradients_2021v1_1.csv"
g4.chla.file <- "Meta_Data/Enviro_Meta_Data/TN397/Gradients4_TN397_Optics_LISST_ACS_ECO.csv"

#Read in datasets:
g4.nut.uw.dat <- read_csv(g4.nut.uw.file)
g4.pc.uw.dat <- read_csv(g4.pc.uw.file)
g4.nut.disc.dat <- read_csv(g4.nut.pc.discrete.file)
g4.ts.dat <- read_csv(g4.ts.file)
g4.bact.dat <- read_csv(g4.bact.file)
g4.chla.dat <- read_csv(g4.chla.file)

g4.pc.1 <- g4.pc.uw.dat %>%
  mutate(time = as_datetime(time),
         time = round_date(time, unit = "hour"))


ggplot(g4.pc.1, aes(x = dt, y = pc, color = lat)) + 
  scale_color_viridis() +
  geom_point(size = 3)

ggplot(g4.pc.1, aes(x = dt, y = pp)) + 
  geom_point()



##surface.discrete.samps 
surf.disc <- g4.nut.disc.dat %>%
  filter(depth < 10)

ggplot(surf.disc, aes(x = time, y = NH4, color = Station)) +
  geom_point()


##g4 ts plots

g4.ts.dat.sml <- g4.ts.dat %>%
  mutate(dt = as_datetime(time),
         dt = round_date(dt, unit = "hour")) %>%
  group_by(dt) %>%
  mutate(mean.sst = mean(SST),
         mean.sal = mean(salinity),
         mean.lat = mean(lat),
         mean.lon = mean(lon)) %>%
  select(dt, mean.sst, mean.sal, mean.lat, mean.lon) %>%
  unique()

ggplot(g4.ts.dat.sml, aes(x = mean.lat, y = mean.sal, color = mean.lat)) +
  geom_point() +
  scale_color_viridis()

ggplot(g4.ts.dat.sml, aes(x = mean.lat, y = mean.sst, color = mean.lat)) +
  geom_point() +
  scale_color_viridis()


library(viridis)
ggplot(g4.ts.dat.sml, aes(x = mean.sal, y = mean.sst, color = mean.lat)) +
  geom_point() +
  scale_color_viridis()

###
ts.sml <- g4.ts.dat.sml %>%
  mutate(time = round_date(dt, unit = "hour")) %>%
  group_by(time) %>%
  mutate(mean.sst = mean(mean.sst),
         mean.sal = mean(mean.sal),
         mean.lat = mean(mean.lat),
         mean.lon = mean(mean.lon)) %>%
  select(time, mean.sst, mean.sal, mean.lat, mean.lon) %>%
  unique()
  



###G4 heterotrophic bacteria plots:
g4.bact.dat.2 <- g4.bact.dat %>%
  filter(depth < 10) %>%
  select(time, lat, lon, station, cell_abundance_bacteria, carbon_biomass_bacteria)

ggplot(g4.bact.dat.2, aes(x = lat, y = cell_abundance_bacteria, color = station)) +
  geom_point(size = 3) 

ggplot(g4.bact.dat.2, aes(x = time, y = cell_abundance_bacteria, color = station)) +
  geom_point(size = 3) 

ggplot() +
  geom_point(data = g4.bact.dat.2, aes(y = lat, size = cell_abundance_bacteria, x = lon)) +
  geom_point(data = G4.kin.exps, aes(x = Long, y = Lat), color = "red", size = 2)

G4.kin.exps <- read_csv("Meta_Data/Ingalls_Lab_Data/Hom_Kin_Locations.csv") %>%
  filter(Compound == "Homarine")



####G4 Chla
g4.chla.dat.2 <- g4.chla.dat %>%
  filter(chla_acs < 0.4) %>%
  mutate(time = as_datetime(time),
         time = round_date(time, unit = "hour")) %>%
  group_by(time) %>%
  reframe(mean_chla_acs = mean(chla_acs),
          mean_lat = mean(lat)) %>%
  unique()

ggplot(g4.chla.dat.2, aes(x = lat, y = chla_acs)) +
  geom_point(alpha = 0.1) +
  geom_smooth()

ggplot(g4.chla.dat.2, aes(x = time, y = chla_acs, color = lat)) +
  geom_point() +
  scale_color_viridis()

chla.sml <- g4.chla.dat.2 %>%
  select(time, mean_chla_acs)
##How to merge all of this data together:

g4.loc.dat <- read_csv("Meta_Data/Ingalls_Lab_Data/G4_Station_Locations.csv") %>%
  unite(c(Local_Date, Local_Time), col = dt1, sep = " ") %>%
  mutate(time = dmy_hms(dt1),
         time = round_date(time, unit = "hour")) %>%
  filter(!is.na(Samp_ID))

g4.loc.2 <- g4.loc.dat %>%
  select(Samp_ID, time, Lat_N, Long_W) %>%
  left_join(., chla.sml) %>%
  left_join(ts.sml)

ggplot(g4.loc.2, aes(color = mean_chla_acs, y = mean.sst, x = mean.sal)) +
  geom_point(size = 4) +
  scale_color_viridis()

####
dat.to.impute.to <- ts.sml %>%
  left_join(., g4.pc.1) %>%
  arrange(time)

x <- na.stinterp(dat.to.impute.to$pc, along = time(dat.to.impute.to$time), na.rm = F)

interp.out <- data_frame(x) %>%
  rename("pc_interp" = x)

dat.test.1 <- dat.to.impute.to %>%
  cbind(., interp.out)

dat.test.1.sml <- dat.test.1 %>%
  filter(mean.lon > -150) %>%
  filter(mean.lon < -138)


ggplot(dat.test.1.sml, aes(x = mean.lat, y = pc_interp)) +
  geom_point(aes(color = time)) +
  #geom_line(aes(color = time, group = time)) +
  scale_color_viridis() +
  geom_path() +
  geom_point(aes(x = mean.lat, y = pc), color = "red")

g4.pc.uw.dat.2 <- g4.pc.uw.dat %>%
  filter(lon > -150) %>%
  filter(lon < -138) 

ggplot(g4.pc.uw.dat.2, aes(x = lat, y = pc, color = time)) +
  geom_point() + 
  scale_color_viridis() +
  geom_path()



#install.packages("imputeTS")
#library(imputeTS)
#ggplot()
#install.packages("forecast")
#install.packages('forecast', dependencies = TRUE)

#install.packages("stinepack")
library(stinepack)
#library(impute)

#library(imputeTS)
#dat.x <- 
