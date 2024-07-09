


#packages 
library(tidyverse)

#Strategy --> Calculate batch specific RF for labeled Homarine and GBt from IS
        # --> use this RF value + EE to quantify GBT in that one sample

###Define inputs
norm.dat.file <- "Intermediates/Dissolved_HILIC_Pos_BMISed_dat.csv"
rf.file <- "Intermediates/Dissolved_Stds_RFs_RFratios.csv"
EE.file <- "Meta_Data/Ingalls_Lab_Data/CXSPE_EEs.csv"
Std.file <- "Meta_Data/Ingalls_Lab_Data/Ingalls_Lab_Standards_03172023.csv"
IS.names.file <- "Meta_Data/Ingalls_Lab_Data/CXSPE_IS_List_JSS.csv" 
HILIC.raw.dat <- "Intermediates/dissolved_betaine_data_raw.csv"
Blk.LOD.dat <- "Intermediates/Dissolved_Blk_LOD.csv"
IS.raw.file <- "Intermediates/Dissolved_IS_data_raw.csv"

#Part 1_________________Calculate batch specific RF values for IS Homarine and GBt

#Get IS information
IS_names <- read_csv(IS.names.file) %>%
  rename(Compound = Match.New,
         IS = IS.Name.New) %>%
  select(Compound, IS, Spike.Fraction, Conc.in.vial.uM)

#bring in raw IS data and combine with IS matching
IS.dat.full <- read_csv(IS.raw.file) %>%
  rename(IS = Compound,
         IS_Area = Area)

#combine IS info and data
IS.dat.named <- left_join(IS.dat.full, IS_names)

#Calculate RF values in uM/area
IS.rf.vals <- IS.dat.named %>%
  filter(!str_detect(Rep, "Std")) %>%
  filter(!str_detect(Rep, "Blk")) %>%
  filter(!str_detect(Rep, "Poo")) %>%
  filter(!str_detect(Rep, "KM1906_GBT_F2")) %>%
  group_by(Cruise, Compound, IS) %>%
  reframe(IS.rf = IS_Area/Conc.in.vial.uM,
          mean.IS.rf = mean(IS.rf),
          sd.IS.rf= sd(IS.rf)) %>%
  select(Cruise, Compound, IS, mean.IS.rf, sd.IS.rf) %>%
  unique() %>%
  filter(Compound %in% c("Glycine betaine", "Homarine"))


#Pull in EE information and combine with RF information:
EE.rf.dat <- read_csv(EE.file) %>%
  rename("EE" = `Extraction Efficiency (%)`) %>%
  select(Compound, EE) %>%
  right_join(IS.rf.vals)




#Part 2___________________Quantify G3 sample containing labeled GBT spike (KM1906_GBT_F2):

#Pull in raw data for experiment and select just the samples of interest:
raw.dat <- read_csv(HILIC.raw.dat) %>%
  filter(str_detect(Rep, "KM1906_GBT_F2")) 


#Quantify GBt in these samples:
g3.samp.quant <- raw.dat %>%
  left_join(., EE.rf.dat) %>%
 # filter(Compound == "Glycine betaine") %>%
  mutate(nmol.conc = Area/mean.IS.rf*10^-6*400/(40*10^-3)*1000,
         EE.adjust.conc = nmol.conc/(EE/100))

write_csv(g3.samp.quant, file = "Intermediates/GBT_Quant_for_KM1906_GBT_F2_edge_case.csv")






























































