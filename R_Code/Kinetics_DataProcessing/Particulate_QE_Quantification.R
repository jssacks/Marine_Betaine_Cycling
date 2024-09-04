


######load in packages
library(tidyverse)
source("R_Code/Functions.R")



#####Define inputs
qe.file <- "Raw_Data/Kinetics_Data/Kinetics_Particulate_QE_skyline_output.csv"
md.file <- "Intermediates/Kinetics_Metadata_for_Part_Quant.csv"
###Need a file with Dilution Factors and Sample Volumes


#Load in a tidy up metadata:
meta.dat <- read_csv(md.file) %>%
  select(Cruise, Samp_ID, Exp, Treatment, Vol_Filt_L) %>%
  filter(Treatment %in% c("0", "0nM", "Glucose", "Con")) 


###Load in data, clean up files, add in meta data, and separate into environmental and IS dataframes:
qe.dat <- sky_read(qe.file) %>%
  select(Rep, Compound, Area) %>%
  filter(Compound %in% c("Homarine", "Homarine, 2H3", "Glycine betaine", "Glycine betaine, 13C5, 15N")) %>%
  filter(str_detect(Rep, "_Smp")) 

metab.dat <- qe.dat %>%
  filter(Compound %in% c("Homarine", "Glycine betaine")) %>%
  mutate(match = Compound) %>%
  mutate(Samp_ID = str_remove(Rep, "240807_Smp_")) %>%
  mutate(Samp_ID = str_replace(Samp_ID, "G5_UKG1_", "G5_UKG_")) %>%
  mutate(Samp_ID = str_replace(Samp_ID, "RC104_UCH1_0nM", "RC104_UCH1_Con")) %>%
  mutate(Samp_ID = str_replace(Samp_ID, "TN397_UCH1_0nM", "TN397_UCH1_Gluc")) %>%
  left_join(., meta.dat) %>%
  mutate(Cruise = case_when(str_detect(Samp_ID, "GBT-K") ~ "KM1906",
                            TRUE ~ Cruise),
         Exp = case_when(str_detect(Samp_ID, "GBT-K1") ~ "UKG1",
                         str_detect(Samp_ID, "GBT-K2") ~ "UKG2",
                         TRUE ~ Exp),
         Treatment = case_when(str_detect(Samp_ID, "GBT-K1") ~ "0nM",
                         str_detect(Samp_ID, "GBT-K2") ~ "0nM",
                         TRUE ~ Treatment),
         Vol_Filt_L = case_when(str_detect(Samp_ID, "GBT-K1") ~ 2,
                         str_detect(Samp_ID, "GBT-K2") ~ 2,
                         TRUE ~ Vol_Filt_L)) %>%
  mutate(Samp_ID = case_when(str_detect(Samp_ID, "GBT-K1") ~ "KM1906_UKG1",
                             str_detect(Samp_ID, "GBT-K2") ~ "KM1906_UKG2",
                             TRUE ~ Samp_ID))# %>%
#  mutate(dilution.factor = case_when(str_detect(Samp_ID, "RC104_UKG1") ~ 10,
#                                     str_detect(Samp_ID, "RC078_UKG1") ~ 10,
#                                     str_detect(Samp_ID, "KM1906_UKG1") ~ 10,
#                                     str_detect(Samp_ID, "KM1906_UKG2") ~ 10,
#                                     TRUE ~ 1)) %>%
#  mutate(Area = Area*dilution.factor)


is.dat <- qe.dat %>%
  filter(Compound %in% c("Homarine, 2H3", "Glycine betaine, 13C5, 15N")) %>%
  mutate(match = case_when(Compound == "Homarine, 2H3" ~ "Homarine",
                           Compound == "Glycine betaine, 13C5, 15N" ~ "Glycine betaine")) %>%
  mutate(is.conc.uM = case_when(Compound == "Homarine, 2H3" ~ 4,
                           Compound == "Glycine betaine, 13C5, 15N" ~ 1)) %>%
  rename("IS" = Compound,
         "IS.Area" = Area) %>%
  mutate(RF.IS = IS.Area/is.conc.uM)
  
  
###Quantify GBT and Homarine in Each Sample using IS values
metab.dat.quant <- left_join(metab.dat, is.dat) %>%
  mutate(Conc.vial.uM = Area/RF.IS,
         Conc.Smp.nM = Conc.vial.uM*1000*(400*1e-6)/Vol_Filt_L)


###Summarize Quantification Results:
metab.dat.sum <- metab.dat.quant %>%
  select(Compound, Cruise, Exp, Conc.Smp.nM) %>%
  group_by(Compound, Cruise, Exp) %>%
  reframe(mean.part.conc.nM = mean(Conc.Smp.nM),
          sd.part.conc.nM = sd(Conc.Smp.nM))

#write to a csv:
write_csv(metab.dat.sum, file = "Intermediates/Kinetics_QE_Particulate_Quantification_Results.csv")





















































































