


library(tidyverse)


####Organize D1 metadata:


#define inputs:

#sample name information:
samp.file <- "Meta_Data/Ingalls_Lab_Data/RC078_metadata.csv"

#lat and long information:
loc.file <- "Meta_Data/Ingalls_Lab_Data/RC078_ctd_downcast_tidy.csv"

#PCPN
pcpn.file <- "Meta_Data/Enviro_Meta_Data/RC078/pocpn_tidy.csv"

#nutrients.file
nut.file <- "Meta_Data/Enviro_Meta_Data/RC078/nutrients_tidy.csv"

#ctd depth info
#ctd.file <- "Meta_Data/Enviro_Meta_Data/RC078/ctd_tidy.csv"


##Read in and tidy up files files

#Sample information
samp.dat <- read_csv(samp.file) %>%
  select(parent_id, sample_id, station, depth_m) %>%
  unique() %>%
  filter(!is.na(depth_m))

#Location and CTD information (Average over top 10 m of water column to get station information)
ctd.dat <- read_csv(loc.file) %>%
  select(Station, Cast, latitude, longitude, depth_m, practical_salinity_psu, temperature_degC, oxygen_V, beam_attenuation_pm, fluorescence_mgpm3) %>%
  filter(depth_m >= 10) %>%
  group_by(Station, Cast) %>%
  reframe(lat = mean(latitude),
          lon = mean(longitude),
          sal = mean(practical_salinity_psu),
          temp = mean(temperature_degC),
          oxygen = mean(oxygen_V),
          beam_atten = mean(beam_attenuation_pm),
          Chl_fluor = mean(fluorescence_mgpm3)) %>%
  mutate(cast = as.numeric(str_remove(Cast, "C")),
         station = as.numeric(str_remove(Station, "S"))) 

#Duplicate S7_C2 data for S7_C3 to provide some contextualizing information for that cast and add to ctd.dat
ctd.s7.c3 <- ctd.dat %>%
  filter(station == 7,
         cast == 2) %>%
  mutate(cast = 3)

ctd.dat <- rbind(ctd.dat, ctd.s7.c3)


#POC and PN data (convert to uM to standardize with other measurments)
pcpn.dat <- read_csv(pcpn.file) %>%
  select(parent_id, sample_id, depth_m, station, cast, c_mol_l, n_mol_l) %>%
  mutate(POC_uM = c_mol_l*1E6,
         PN_uM = n_mol_l*1E6) 

#Nutrient data (pivot to wide format)
nut.dat <- read_csv(nut.file) %>%
  select(parent_id, sample_id, analyte_short, concentration) %>%
  pivot_wider(id_cols = c(parent_id, sample_id), names_from = analyte_short, values_from = concentration)


###Combine together:
all.dat <- pcpn.dat %>%
  left_join(., nut.dat) %>%
  left_join(., ctd.dat)
  
#write to csv
write_csv(all.dat, file = "Intermediates/RC078_MetaData_Compiled.csv")


