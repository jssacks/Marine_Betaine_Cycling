




library(readr)
library(tidyverse)



####Source Functions
source("R_Code/Functions.R")

#Define Inputs:
part.file <- "Intermediates/particulate_betaine_data_raw.csv"

Stds.info.file <- "Meta_Data/Ingalls_Lab_Data/Ingalls_Lab_Standards_03172023.csv" 


#Load in data
part.dat <- read_csv(part.file) %>%
  rename("SampID" = Rep,
         "Name" = Compound) %>%
  filter(str_detect(.$SampID, "Std"))
  

####Load in standards
stds.dat <- part.dat %>%
  mutate(Mix = str_extract(SampID, "Mix\\d")) 

###Get Mix and Concentration info:
stds.info <- read_csv(Stds.info.file) %>%
  filter(Priority == TRUE) %>%
  select(Compound_Name, z, Column, HILIC_Mix, Concentration_uM) %>%
  rename("Name" = Compound_Name) 

###Join stuff together + remove Matrix Samples
stds.dat.info <- left_join(stds.dat, stds.info) 

##Calculate RFs
RF.dat <- stds.dat.info %>%
  filter(Mix == HILIC_Mix) %>%
  select(-Mix, -HILIC_Mix) %>%
  filter(!str_detect(.$SampID, "Matrix")) %>%
  mutate(RF = as.numeric(Area)/Concentration_uM, NA) %>%
  group_by(Name, Cruise) %>%
  summarise(RFmax = max(RF),
            RFmin = min(RF),
            RF = mean(RF, na.rm = TRUE))  %>%
  ungroup()

RFratio.dat <- stds.dat.info %>%
  filter(HILIC_Mix == Mix | is.na(Mix)) %>% 
  mutate(RunNumber = str_extract(SampID, "_\\d$")) %>%
  mutate(RunType = ifelse(str_detect(SampID, "StdsMix\\dInH2O")|
                            str_detect(SampID, "StdsInH2O"), "Std_in_h2O", 
                          ifelse(str_detect(SampID, "StdsMix\\dInMatrix") |
                                   str_detect(SampID, "StdsInMatrix"), "Std_in_matrix",
                                 "Matrix_in_h2O"))) %>%
  filter(HILIC_Mix == Mix | is.na(Mix)) %>% 
  select(-Mix, -HILIC_Mix, -SampID, -Concentration_uM) %>%
  spread(key = RunType, value = Area ) 

RF.ratios <- RFratio.dat %>%
  ungroup() %>%
  group_by(Name, RunNumber, Cruise) %>%
  summarize("Std_in_matrix" = Std_in_matrix,
            "Matrix_in_h2o" = Matrix_in_h2O,
            "Std_in_h2o" = Std_in_h2O) %>%
  mutate(RFratio = (Std_in_matrix - Matrix_in_h2o)/ Std_in_h2o) %>%
  group_by(Name, Cruise) %>%
  summarise(RFratio = mean(RFratio)) %>%
  ungroup()

###Join it all together
RF.RFratios <- left_join(RF.dat, RF.ratios)
write_csv(RF.RFratios, file = "Intermediates/Particulate_Stds_RFs_RFratios.csv")
