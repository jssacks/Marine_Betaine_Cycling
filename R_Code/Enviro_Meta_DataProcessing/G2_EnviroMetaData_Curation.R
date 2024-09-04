
#load packages
library(tidyverse)


#Define inputs:

#Metabolite sample information:
metab.file <- "Meta_Data/Enviro_Meta_Data/MGL1704/CTD_Info_for_G2_Samps.csv"


#MGL1704 (Gradeints 2)
pc.uw.file <- "Meta_Data/Enviro_Meta_Data/MGL1704/Gradients2_MGL1704_PPPCPN_UW.csv"
nut.bottle.file <- "Meta_Data/Enviro_Meta_Data/MGL1704/MGL1704_Gradients2_Nutrients.csv"
ctd.file <- "Meta_Data/Enviro_Meta_Data/MGL1704/MGL1704_Gradients2.csv"


#________
#Read in sample info:
metab.info <- read_csv(metab.file) %>%
  rename("Cast_Number" = Cast_Number_For_MetaData) %>%
  select(-Event_Number)
  

#Read in datasets:
pc.uw.dat <- read_csv(pc.uw.file)
nut.bottle.dat <- read_csv(nut.bottle.file)
ctd.dat <- read_csv(ctd.file)

##Compile and average information at the depth level for environmental parameters

#CTD Dat:
ctd.sml <- ctd.dat %>%
  select(time, stnnbr, castno, depth, ctdtmp, ctdsal, chla, lat, lon) %>%
  mutate(depth = round(depth)) %>% #round depth to nearest whole meter
  mutate(time = round_date(time, unit = "day")) %>% #round time to day to match with pcpn underway measurments
  group_by(stnnbr, castno) %>%
  mutate(chla = mean(chla, na.rm = T)) %>%
  group_by(stnnbr, castno, depth) %>%
  mutate(ctdtmp = mean(ctdtmp, na.rm = T),
         ctdsal = mean(ctdsal, ra.rm = T)) %>%       #chla is Fluorometric Chlorophyll a measured by A. White
  ungroup() %>%
  unique() %>%
  rename("Station_Number" = stnnbr,
         "Cast_Number" = castno,
         "Depth" = depth)

#Nutrient Dat:
nut.sml <- nut.bottle.dat %>%
  select(Station, NO3_NO2, PO4) %>%
  filter(!str_detect(Station, "Underway")) %>%
  mutate("Station_Number" = as.numeric(Station)) %>%
  select(-Station) %>%
  mutate(Depth = 15)


#PCPN Dat:
pcpn.sml <- pc.uw.dat %>%
  mutate(time = round_date(time, unit = "day")) %>%
  group_by(time) %>%
  mutate(pc = mean(pc),
         pn = mean(pn)) %>%
  select(time, pc, pn) %>%
  unique()

###Comibine PCPN data with CTD based on date and nutrient data based on station number and depth:
ctd.pcpn.nut <- ctd.sml %>%
  mutate(Depth = case_when(Depth == 17 ~ 15,           #adjust slight difference in depth to match with other metadata
                           TRUE ~ Depth)) %>%
  left_join(pcpn.sml) %>%
  left_join(., nut.sml) 


###Combine enviro metadata with Sample Information:
g2.enviro.full <- left_join(metab.info, ctd.pcpn.nut)














