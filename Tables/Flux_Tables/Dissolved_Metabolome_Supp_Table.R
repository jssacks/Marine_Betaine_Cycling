


library(tidyverse)


#define inputs:
diss.file <- "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv"
lod.file <- "Intermediates/Dissolved_Blk_LOD_Concentrations.csv"


#Load in datasets:

#load in datasets and match up samples with the batch they were run 
# with to be able to match samples with batch specific LOD for each compound
diss.dat <- read_csv(diss.file) %>% 
  mutate(Batch = case_when(Cruise == "TN397" ~ "TN397",
                           Cruise == "KM1906" ~ "TN397", 
                           Cruise == "RC078" ~ "KinExp",
                           Cruise == "RC104" ~ "KinExp",
                           Cruise == "TN412" ~ "KinExp"))

lod.dat <- read_csv(lod.file) %>%
  rename("Batch" = Cruise,
         "Compound" = Name)

#combine togther and organize
diss.lod.table <- diss.dat %>%
  left_join(., lod.dat) %>%
  mutate(Mean.Diss.Conc.nM = print(formatC(signif(Mean.Diss.Conc.nM, digits=3), digits=3,format="fg", flag = "#")),
         SD.Diss.Conc.nM = print(formatC(signif(SD.Diss.Conc.nM, digits=3), digits=3,format="fg", flag = "#")),
         LOD.nM = print(formatC(signif(EE.adjust.lod, digits=3), digits=3,format="fg", flag = "#"))) %>%
  select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, LOD.nM) %>%
  rename("Experiment" = exp,
         "Average Dissolved Concentration (nM)" = Mean.Diss.Conc.nM,
         "Standard Deviation of Dissolved Concentration (nM)" = SD.Diss.Conc.nM,
         "Batch Limit of Detection (nM)" = LOD.nM)

#export:
write_csv(diss.lod.table, file = "Tables/Flux_Tables/Dissolved_Environmental_Metabolomes_LODs_Supplemental_Table.csv")

    
    
    
    
    
    






















































































