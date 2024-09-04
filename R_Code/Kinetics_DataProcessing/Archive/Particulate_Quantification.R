

#
#This script....


#####
library(readr)
library(tidyverse)

#Inputs

#Stds and RFs and RFratios
stds.file <- "Meta_Data/Ingalls_Lab_Data/Ingalls_Lab_Standards_03172023.csv" 
rf.file <- "Intermediates/Particulate_Stds_RFs_RFratios.csv"

#IS
hilic.is.file <- "Intermediates/particulate_IS_data_raw.csv"

### area dat
hilic.file <- "Intermediates/Particulate_HILIC_Pos_BMISed_dat.csv"
hilic.file.notnorm <- "Intermediates/particulate_betaine_data_raw.csv"


###Sample volume dat:
g4.vol.file  <- "Meta_Data/Ingalls_Lab_Data/G4_vol_filt.csv"
g3.vol.file <- "Meta_Data/Ingalls_Lab_Data/G3_vol_filt.csv"
d1.vol.file <- "Meta_Data/Ingalls_Lab_Data/RC078_metadata.csv"



#Quantify values using standards 
rf.dat <- read_csv(rf.file)


##Pull in and combine volume filtered data:

#G4
g4.vol.filt <- read_csv(g4.vol.file) %>%
  mutate(SampID = paste("220902_Smp_", Samp_ID, sep = "")) %>%
  rename("Vol_L" = Vol_Filt_L) %>%
  select(SampID, Vol_L) %>%
  mutate(Cruise = "TN397")

#G3
g3.vol.filt <- read_csv(g3.vol.file) %>%
  rename("SampID" = `Sample ID`,
         "Vol_L" = `Vol Filtered (L)`) %>%
  filter(str_detect(SampID, "MU")) %>%
  mutate(SampID = paste("220628_Smp_", SampID, sep = "")) %>%
  mutate(Vol_L = as.numeric(str_remove(Vol_L, "L"))) %>%
  mutate(Cruise = "KM1906")

#D1
d1.vol.filt <- read_csv(d1.vol.file) %>%
  filter(sample_type == "METABS") %>%
  select(sample_id, parent_id, station, depth_m, cast, niskin, vol_filt_l) %>%
  rename("SampID" = sample_id,
         "Vol_L" = vol_filt_l) %>%
  mutate(SampID = paste("221006_Smp_", SampID, sep = "")) %>%
  select(SampID, Vol_L) %>%
  mutate(Cruise = "RC078")

vol.filt.dat <- rbind(g4.vol.filt, g3.vol.filt, d1.vol.filt)



### combine metab area data with vol.filt.data and 
## calculate concentration in Vial using normal RF and RFratio approach 
vial.quant.dat <- read_csv(hilic.file) %>%
  rename("Name" = MF) %>%
  left_join(., vol.filt.dat) %>%
  filter(!is.na(Vol_L)) %>%
  left_join(., rf.dat) %>%
  mutate(umol.in.vial.ave = Adjusted_Area/RF/RFratio,
         umol.in.vial.max = Adjusted_Area/RFmin/RFratio,
         umol.in.vial.min = Adjusted_Area/RFmax/RFratio) 





#Quantify compounds with matched IS
hilic.is.dat <- read_csv(hilic.is.file)  %>%
  rename("IS" = Compound,
         "IS.area" = Area,
         "SampID" = Rep) %>%
  unique() %>%
  filter(!SampID == "221006_Smp_S7_C1_D1_A")



#get internal standard concentration values 
is.std <- read_csv(stds.file) %>%
  filter(Compound_Name %in% hilic.is.dat$IS) %>%
  select(Compound_Name, Concentration_uM) %>%
  rename("IS" = Compound_Name,
         "IS.Conc.uM" = Concentration_uM) 


####Pull in unnormalized hilic data 
hilic.notnorm.dat <- read_csv(hilic.file.notnorm) %>%
  rename("SampID" = Rep,
         "Name" = Compound) %>%
  select(Name, SampID, Area) %>%
  unique() %>%
  full_join(., hilic.is.dat)  %>%
  mutate(match.name = str_extract(.$IS, Name)) %>%
  filter(!is.na(match.name)) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  filter(!str_detect(.$SampID, "Std")) %>%
  filter(!str_detect(.$SampID, "Blk")) %>%
  filter(!SampID == "221006_Smp_S7_C1_D1_A")



#Quantify using internal standard
vial.is.quant.dat <- left_join(hilic.notnorm.dat, is.std) %>%
  mutate(umol.in.vial.ave = Area*IS.Conc.uM/IS.area) %>%
  left_join(., vol.filt.dat) %>%
  filter(!is.na(Vol_L)) 



smp.quant.dat.all <- vial.quant.dat %>%
  select(Name, SampID, Cruise, umol.in.vial.ave, Vol_L) %>%
  filter(!Name %in% vial.is.quant.dat$Name) %>%
  rbind(., vial.is.quant.dat %>% 
          select(Name, SampID, Cruise, umol.in.vial.ave, Vol_L)) %>%
  mutate(nM.in.smp = umol.in.vial.ave*10^-6*400/(Vol_L)*1000)%>%
  unique()




#Calculate nM C and N per sample
std.formula <- read_csv(stds.file) %>%
  select(Compound_Name, Empirical_Formula) %>%
  rename("Name" = Compound_Name) %>%
  unique() %>%
  mutate(C = ifelse(is.na(str_extract(Empirical_Formula, "^C\\d\\d")),
                    str_extract(Empirical_Formula, "^C\\d"), 
                    str_extract(Empirical_Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Empirical_Formula, "N\\D"),
                    1, str_extract(Empirical_Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", "")))


samp.quant.dat <- left_join(smp.quant.dat.all, std.formula) %>%
  unique() %>%
  mutate(nM_C = nM.in.smp*C,
         nM_N = nM.in.smp*N) %>%
  select(Name, SampID, Cruise, Vol_L, umol.in.vial.ave, nM.in.smp, nM_C, nM_N)

write_csv(samp.quant.dat, file = "Intermediates/Particulate_Quant_Output.csv")















#
