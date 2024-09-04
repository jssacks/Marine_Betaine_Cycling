




##load packages
library(tidyverse)

##Functions:
source("R_Code/Functions.R")


##Define input files:

#G3 Particulate
g3p.file <- "Raw_Data/Distribution_Data/Particulate/G3_Underway_Particulate_HILIC_Pos.csv"

#G4 Particulate
g4p.file <- "Raw_Data/Distribution_Data/Particulate/HILIC_Pos_G4_Particulate_Underway_LTC_July2024.csv"

#D1 Particulate
d1p.file <- "Raw_Data/Distribution_Data/Particulate/Carson2022_Particulate_HILIC_Pos_LTC_July2024.csv"


#____________________________________
#G3 Dissolved
g3d.file <- "Raw_Data/Distribution_Data/Dissolved/G3_Underway_Dissolved_HILIC_Pos.csv"

#G4 Dissolved
g4d.file.0 <-  "Raw_Data/Distribution_Data/Dissolved/G4_CX_HILIC_Pos_Pooled_Output_LTC_July2024.csv"
g4d.file.1 <-  "Raw_Data/Distribution_Data/Dissolved/G4_CX_HILIC_Pos_Group_1_LTC_July2024.csv"
g4d.file.2 <-  "Raw_Data/Distribution_Data/Dissolved/G4_CX_HILIC_Pos_Group_2_LTC_July2024.csv"

#D1 Dissolved
d1d.file <- "Raw_Data/Distribution_Data/Dissolved/Carson2022_CXSPE_HILIC_Pos_LTC_July2024.csv"

#Kinetics Experiments Dissolved
ked.file <-  "Raw_Data/Distribution_Data/Dissolved/CXSPE_Kin_Samps_Oct2023_Betaines_LTC_July2024.csv"








# Organize Particulate Data -----------------------------------------------

#Organize G3-particulate
g3p <- sky_read(g3p.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "KM1906")

#Organize G4-Particulate
g4p <- sky_read(g4p.file) %>% 
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")

#Organize D1-Particulate
d1p <- sky_read(d1p.file) %>% 
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "RC078")


####Pull together all particulate data:
all.p.dat <- rbind(g3p, g4p, d1p)

is.p.dat <- all.p.dat %>%
  filter(str_detect(Compound, ", ")) %>%
  filter(!Compound %in% c("Cys-Gly, oxidized"))

betaine.p.dat <- all.p.dat %>%
  filter(Compound %in% 
               c("beta-Alaninebetaine",
                 "Glycine betaine",
                 "Proline betaine",
                 "Homarine",
                 "Trigonelline",
                 "Betonicine",
                 "Dimethylsulfonioacetate",
                 "Dimethylsulfoniopropionate",
                 "Gonyol", 
                 "Trimethylamine N-oxide"))

#make sample list:
p.smp.list <- betaine.p.dat %>%
  select(Rep, Cruise) %>%
  unique() %>%
  mutate(Injec_vol = case_when(str_detect(Rep, "Half") ~ 1,
                             TRUE ~ 2))

####Write cleaned up IS and betaine data to .csvs 
write_csv(is.p.dat, file = "Intermediates/particulate_IS_data_raw.csv")

write_csv(betaine.p.dat, file = "Intermediates/particulate_betaine_data_raw.csv")

write_csv(p.smp.list, file = "Intermediates/particulate_smp_list.csv")




# Organize Dissolved Data ---------------------------------------------

#Organize G3-dissolved
g3d <- sky_read(g3d.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "KM1906")

#Organize G4-dissolved
g4d.0 <- sky_read(g4d.file.0) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")

g4d.1 <- sky_read(g4d.file.1) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")  %>%
  filter(!str_detect(Rep, "_Std_"))

g4d.2 <- sky_read(g4d.file.2) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")%>%
  filter(!str_detect(Rep, "_Std_"))

g4d <- rbind(g4d.0, g4d.1, g4d.2)


###Organize D1 Dissolved
d1d <- sky_read(d1d.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "RC078")
#XXX



#Organize Kinetics Experiment Dissolved
ked <- sky_read(ked.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "KinExp") 

####Pull together all dissolved data:
all.d.dat <- rbind(g3d, g4d, ked, d1d)

###Grab Just IS:
is.d.dat <- all.d.dat %>%
  filter(str_detect(Compound, ", ")) %>%
  filter(!Compound %in% c("Cys-Gly, oxidized"))


###Grab Just betaines: 
betaine.d.dat <- all.d.dat %>%
  filter(Compound %in% 
           c("beta-Alaninebetaine",
             "Glycine betaine",
             "Proline betaine",
             "Homarine",
             "Trigonelline",
             "Betonicine",
             "Dimethylsulfonioacetate",
             "Dimethylsulfoniopropionate",
             "Gonyol", 
             "Trimethylamine N-oxide"))

#make sample list:
d.smp.list <- betaine.d.dat %>%
  select(Rep, Cruise) %>%
  unique() %>%
  mutate(Injec_vol = case_when(str_detect(Rep, "Half") ~ 1,
                             TRUE ~ 2))

####Write cleaned up IS and betaine data to .csvs 
write_csv(is.d.dat, file = "Intermediates/dissolved_IS_data_raw.csv")

write_csv(betaine.d.dat, file = "Intermediates/dissolved_betaine_data_raw.csv")

write_csv(d.smp.list, file = "Intermediates/dissolved_smp_list.csv")































