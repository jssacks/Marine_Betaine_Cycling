

#load required packages
library(tidyverse)
library(stinepack)
library(lubridate)


#Define inputs:

#TN412 (Gradeints 5)
nut.uw.file <- "Meta_Data/Enviro_Meta_Data/TN412/Gradients5_TN412_NutrientsUW.csv"
#pc.uw.file <- "Meta_Data/Enviro_Meta_Data/TN397/Gradients4_TN397_PPPCPN_UW.csv"
don.file <- "Meta_Data/Enviro_Meta_Data/TN412/TN412_Gradients5_NutrientsWaterColumn.csv"
ts.file <- "Meta_Data/Enviro_Meta_Data/TN412/TN412_Gradients5_uw_tsg.csv"
#bact.file <- "Meta_Data/Enviro_Meta_Data/TN397/Influx_Stations_Gradients_2021v1_1.csv"
chla.file <- "Meta_Data/Enviro_Meta_Data/TN412/Gradients5_TN412_Optics_LISST_ACS_ECO.csv"
pp.13c.file <- "Meta_Data/Enviro_Meta_Data/TN412/Gradients5_TN412_15N13C.csv"

#Read in datasets:
nut.uw.dat <- read_csv(nut.uw.file)
#pc.uw.dat <- read_csv(pc.uw.file)
don.dat <- read_csv(don.file)
ts.dat <- read_csv(ts.file)
#bact.dat <- read_csv(bact.file)
chla.dat <- read_csv(chla.file)
pp.13c.dat <- read_csv(pp.13c.file)


#select desired variables and summarize datasets at the level of hour

##temp and sal
ts.hr <- ts.dat %>%
  mutate(time = round_date(time, unit = "hour")) %>%
  group_by(time) %>%
  mutate(sst = mean(TS),
         sss = mean(SSPS),
         lat = mean(lat),
         lon = mean(lon)) %>%
  select(time, sst, sss, lat, lon) %>%
  unique()

##Chl-a
chla.hr <- chla.dat %>%
  mutate(time = round_date(time, unit = "hour")) %>%
  group_by(time) %>%
  mutate(chla = mean(chla_acs, na.rm = TRUE)) %>%
  select(time, chla) %>%
  ungroup() %>%
  unique() %>%
  filter(!chla == "NaN")

##PC and PN
#pcpn.hr <- pc.uw.dat %>%
#  mutate(time = round_date(time, unit = "hour")) %>%
#  group_by(time) %>%
#  mutate(pc = mean(pc),
#         pn = mean(pn)) %>%
#  select(time, pc, pn) %>%
#  unique()

#PC and PN:


##Nutrients
#get surface bottle nutrient data
nut.hr.bottle <- don.dat %>%
  filter(depth < 10) %>%
  mutate(time = round_date(time, unit = "hour")) %>%
  group_by(time) %>%
  mutate(N_N = mean(NO3_NO2)) %>%
  select(time, N_N) %>%
  unique()

nut.hr.uw <- nut.uw.dat %>%
  filter(depth < 10) %>%
  mutate(time = round_date(time, unit = "hour")) %>%
  group_by(time) %>%
  mutate(N_N = mean(NplusN)) %>%
  select(time, N_N) %>%
  unique()

nut.hr <- rbind(nut.hr.bottle, nut.hr.uw) %>%
  filter(!time == as_datetime("2023-01-31 15:00:00")) %>% #remove outlier
  filter(!time == as_datetime("2023-02-03 00:00:00"))     #remove outlier

#combine bottle and underway nutrient data
#nut.hr <- nut.uw.dat %>%
#  mutate(time = round_date(time, unit = "hour")) %>%
#  group_by(time) %>%
#  mutate(N_N = mean(Nitrate_Nitrite),
#         SRP = mean(SRP)) %>%
#  select(time, N_N, SRP) %>%
# unique()

##DON
#don.hr <- don.dat %>%
#  filter(depth < 10) %>%
#  mutate(time = round_date(time, unit = "hour")) %>%
#  group_by(time) %>%
#  mutate(DON = mean(DON)) %>%
#  select(time, DON) %>%
#  unique()

##Bacteria
#bact.hr <- bact.dat %>%
#  filter(depth < 10) %>%
#  mutate(time = round_date(time, unit = "hour")) %>%
#  group_by(time) %>%
#  mutate(bact_abu = cell_abundance_bacteria,
#         bact_c = carbon_biomass_bacteria) %>%
#  select(time, bact_abu, bact_c) %>%
#  unique()

###PP from 13C measurements
pp.13c.hr <- pp.13c.dat%>%
  filter(depth < 10) %>%
  mutate(time = round_date(time, unit = "hour")) %>%
  group_by(time) %>%
  mutate(pp_13c = mean(avg_fix_rate_13C)) %>%
  select(time, pp_13c)


###_________Put it all together and attempt interpolation
full.dat <- ts.hr %>%
  left_join(., chla.hr) %>%
#  left_join(., pcpn.hr) %>%
  left_join(., nut.hr) %>%
  left_join(., pp.13c.hr) %>%
#  left_join(., don.hr) %>%
#  left_join(., bact.hr) %>%
  arrange(time)


##Use stineman interpolation to interpolate values for all variables 

#chla
chla.interp <- data_frame(chla_interp = na.stinterp(full.dat$chla, along = time(full.dat$time), na.rm = F))

#pc
#pc.interp <- data_frame(pc_interp = na.stinterp(full.dat$pc, along = time(full.dat$time), na.rm = F))

#pn
#pn.interp <- data_frame(pn_interp = na.stinterp(full.dat$pn, along = time(full.dat$time), na.rm = F))

#N_N
N_N.interp <- data_frame(N_N_interp = na.stinterp(full.dat$N_N, along = time(full.dat$time), na.rm = F))

#SRP
#SRP.interp <- data_frame(SRP_interp = na.stinterp(full.dat$SRP, along = time(full.dat$time), na.rm = F))

#DON
#DON.interp <- data_frame(DON_interp = na.stinterp(full.dat$DON, along = time(full.dat$time), na.rm = F))

#bact_abu
#bact_abu.interp <- data_frame(bact_abu_interp = na.stinterp(full.dat$bact_abu, along = time(full.dat$time), na.rm = F))

#bact_c
#bact_c.interp <- data_frame(bact_c_interp = na.stinterp(full.dat$bact_c, along = time(full.dat$time), na.rm = F))

#13C PP
pp_13c.interp <- data_frame(pp_13c_interp = na.stinterp(full.dat$pp_13c, along = time(full.dat$time), na.rm = F))

###Put all datasets together:
full.interp <- full.dat %>%
  cbind(., chla.interp, 
        #pc.interp, pn.interp, 
        N_N.interp, pp_13c.interp
        #SRP.interp, DON.interp, bact_abu.interp, bact_c.interp
        )

#write to csv
write_csv(full.interp, file = "Intermediates/G5_MetaData_Interpolated.csv")






#Prelim_Viz
ggplot(full.interp, aes(x = lat, y = sss)) +
  geom_point(aes(color = time)) # +
 # geom_point(aes(x = lat, y = N_N), color = "red")

