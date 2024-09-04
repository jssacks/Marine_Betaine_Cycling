

#
#This script....


#####
library(readr)
library(tidyverse)

#Inputs

#Stds and RFs and RFratios
stds.file <- "Meta_Data/Ingalls_Lab_Data/Ingalls_Lab_Standards_03172023.csv" 
rf.file <- "Intermediates/culture_Stds_RFs_RFratios.csv"

#IS
hilic.is.file <- "Intermediates/culture_IS_data_raw.csv"

### area dat
hilic.file <- "Intermediates/Culture_HILIC_Pos_BMISed_dat.csv"
hilic.file.notnorm <- "Intermediates/culture_betaine_data_raw.csv"


###Sample volume dat:
#file goes here when ready....



#Quantify values using standards 
rf.dat <- read_csv(rf.file)


##Pull in and combine volume filtered data:
#WILL DO THIS WHEN I ACTUALLY HAVE VOLUME DATA FOR EVERYTHING


### combine metab area data with vol.filt.data and 
## calculate concentration in Vial using normal RF and RFratio approach 
vial.quant.dat <- read_csv(hilic.file) %>%
  rename("Name" = MF) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  filter(!str_detect(.$SampID, "Std")) %>%
 # left_join(., vol.filt.dat) %>%
 # filter(!is.na(Vol_L)) %>%
  left_join(., rf.dat) %>%
  mutate(uM.in.vial.ave = Adjusted_Area/RF/RFratio,
         uM.in.vial.max = Adjusted_Area/RFmin/RFratio,
         uM.in.vial.min = Adjusted_Area/RFmax/RFratio) 



#Quantify compounds with matched IS
hilic.is.dat <- read_csv(hilic.is.file)  %>%
  rename("IS" = Compound,
         "IS.area" = Area,
         "SampID" = Rep) %>%
  unique() 


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
  filter(!str_detect(.$SampID, "SAR11")) %>%
  mutate(Area = replace_na(Area, 0))



#Quantify using internal standard
vial.is.quant.dat <- left_join(hilic.notnorm.dat, is.std) %>%
  mutate(uM.in.vial.ave = Area*IS.Conc.uM/IS.area)


smp.quant.dat.all <- vial.quant.dat %>%
  select(Name, SampID, Batch, uM.in.vial.ave) %>%
  filter(!Name %in% vial.is.quant.dat$Name) %>%
  rbind(., vial.is.quant.dat %>% 
          select(Name, SampID, Batch, uM.in.vial.ave)) %>%
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
  mutate(uM_C_vial = uM.in.vial.ave*C,
         uM_N_vial = uM.in.vial.ave*N) %>%
  select(Name, SampID, Batch, uM.in.vial.ave, uM_C_vial, uM_N_vial)

write_csv(samp.quant.dat, file = "Intermediates/Culture_Quant_Output.csv")
