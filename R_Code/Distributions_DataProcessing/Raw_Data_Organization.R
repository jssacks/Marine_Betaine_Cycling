




##load packages
library(tidyverse)

##Functions:
source("R_Code/Functions.R")


##Define input files:


#________________Particulate___________
#G3 Particulate
g3p.file <- "Raw_Data/Distribution_Data/Particulate/G3_Underway_Particulate_HILIC_Pos.csv"

#G4 Particulate
g4p.file <- "Raw_Data/Distribution_Data/Particulate/HILIC_Pos_G4_Particulate_Underway_LTC_July2024.csv"

#D1 Particulate
d1p.file <- "Raw_Data/Distribution_Data/Particulate/Carson2022_Particulate_HILIC_Pos_LTC_July2024.csv"


#_________________Dissolved_____________
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


#_________________Culture_____________
bact.file <- "Raw_Data/Culture_Data/HILIC_POS_QE_CultureRerun_BacteriaBatch.csv"
cyano.file <- "Raw_Data/Culture_Data/HILIC_POS_QE_CultureRerun_CyanoBatch.csv"
diatom.file <- "Raw_Data/Culture_Data/HILIC_POS_QE_CultureRerun_DiatomBatch.csv"
dinogreen.file <- "Raw_Data/Culture_Data/HILIC_POS_QE_CultureRerun_DinoGreenBatch.csv"
hapto.file <- "Raw_Data/Culture_Data/HILIC_POS_QE_CultureRerun_HaptophytesBatch.csv"



#_________G2 Size Fraction Data
g2.file <- "Raw_Data/SizeFractionated_Data/HILIC_QE_POS_G2_Josh_July2024.csv"





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



# Organize Culture Data _____________________________________

#Organize data from each batch
bact.dat <- sky_read(bact.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Batch = "Bacteria")

cyano.dat <- sky_read(cyano.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Batch = "Cyano")

diatom.dat <- sky_read(diatom.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Batch = "Diatom")

dinogreen.dat <- sky_read(dinogreen.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Batch = "Dino_Green")

hapto.dat <- sky_read(hapto.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Batch = "Haptophyte")

####Pull together all datasets together:
all.c.dat <- rbind(bact.dat, cyano.dat, diatom.dat, dinogreen.dat, hapto.dat)

###Grab Just IS:
is.c.dat <- all.c.dat %>%
  filter(str_detect(Compound, ", ")) %>%
  filter(!Compound %in% c("Cys-Gly, oxidized"))

###Grab Just betaines: 
betaine.c.dat <- all.c.dat %>%
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
c.smp.list <- betaine.c.dat %>%
  select(Rep, Batch) %>%
  unique() %>%
  mutate(Injec_vol = case_when(str_detect(Rep, "Half") ~ 1,
                               TRUE ~ 2))

####Write cleaned up IS and betaine data to .csvs 
write_csv(is.c.dat, file = "Intermediates/culture_IS_data_raw.csv")

write_csv(betaine.c.dat, file = "Intermediates/culture_betaine_data_raw.csv")

write_csv(c.smp.list, file = "Intermediates/culture_smp_list.csv")



#________Organize Size Fractionated G2 Data

#organize data:
g2.dat <- sky_read(g2.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Batch = "G2")

###Grab Just IS:
is.g2.dat <- g2.dat %>%
  filter(str_detect(Compound, ", ")) %>%
  filter(!Compound %in% c("Cys-Gly, oxidized"))

###Grab just betaines: (no TMAO because instrument not set to detect low masses)
betaine.g2.dat <- g2.dat %>%
  filter(Compound %in% 
           c("beta-Alaninebetaine",
             "Glycine betaine",
             "Proline betaine",
             "Homarine",
             "Trigonelline",
             "Betonicine",
             "Dimethylsulfonioacetate",
             "Dimethylsulfoniopropionate",
             "Gonyol"))

#make sample list:
g2.smp.list <- betaine.g2.dat %>%
  select(Rep, Batch) %>%
  unique() %>%
  mutate(Injec_vol = case_when(str_detect(Rep, "Half") ~ 1,
                               TRUE ~ 2))

####Write cleaned up IS and betaine data to .csvs 
write_csv(is.g2.dat, file = "Intermediates/G2_IS_data_raw.csv")

write_csv(betaine.g2.dat, file = "Intermediates/G2_betaine_data_raw.csv")

write_csv(g2.smp.list, file = "Intermediates/G2_smp_list.csv")















