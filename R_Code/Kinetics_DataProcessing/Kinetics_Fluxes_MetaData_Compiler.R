



library(tidyverse)




####Define inputs:
flux.file <- "Intermediates/Flux_Analysis_Output.csv"
kin.file <- "Intermediates/Kin_Analysis_Output.csv"
part.conc.file <- "Intermediates/Kinetics_QE_Particulate_Quantification_Results.csv"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"
enviro.conc.file <- "Intermediates/Dissolved_Final_Quant_QCed.csv"
#pp.file <- 



##______Load in and combine flux+TT and Kinetics data:

#load in flux data
flux.dat <- read_csv(flux.file) %>%
  rename("Cruise" = cruise)
  #separate(cruise_exp, into = c("Cruise", "exp"), remove = F)

#load in kinetics data
kin.dat <- read_csv(kin.file)

#combine kinetics and flux data
flux.kin.dat <- left_join(flux.dat, kin.dat)

write_csv(flux.kin.dat, file = "Intermediates/Compiled_Kinetics_Fluxes.csv")







####Compare Fluxes with Primary Productivity and Particulate Concentrations:

##Load in primary productivity data:
pp.dat <- read_csv(meta.data.file) %>%
  select(KinExp_ID, pp_13c, pp_14c) %>% 
  mutate(pp_13c = replace_na(pp_13c, 0),
         pp_14c = replace_na(pp_14c, 0)) %>%
  mutate(pp.uM.C.perday = case_when(pp_13c > pp_14c ~ pp_13c,
                            pp_14c > pp_13c ~ pp_14c)) %>%
  select(KinExp_ID, pp.uM.C.perday)



###________Load in and combine all particulate Concentration and Contextualizing Data for Fluxes:

####load in particulate concentration data:
p.conc.dat <- read_csv(part.conc.file) %>%
  mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT")) %>%
  rename("exp" = Exp)

#Flux data:
flux.sum <- flux.kin.dat %>%
  select(cruise_exp, Cruise, exp, Compound, mm_flux_nM_day, mm_flux_sd_nM_day) %>%
  mutate(C = case_when(Compound == "Homarine" ~7,
                       Compound == "GBT" ~ 5)) %>%
  mutate(flux_nM_C_day = mm_flux_nM_day*C,
         flux_sd_C = mm_flux_sd_nM_day*C) %>%
  rename("KinExp_ID" = cruise_exp) %>%
  left_join(., p.conc.dat) %>%
  mutate(Flux_Perc_of_Part = (mm_flux_nM_day/mean.part.conc.nM)*100,
         Flux_Perc_of_Part_Error = Flux_Perc_of_Part*sqrt((mm_flux_sd_nM_day/mm_flux_nM_day)^2+(sd.part.conc.nM/mean.part.conc.nM)^2)) %>%
  left_join(., pp.dat) %>%
  mutate(flux_Perc_of_PP = (flux_nM_C_day/pp.uM.C.perday)*100,
         flux_Perc_of_PP_Error = flux_Perc_of_PP*sqrt((flux_sd_C/flux_nM_C_day)^2))

####export
write_csv(flux.sum, file = "Intermediates/flux_contextualizaing_data")






#Compile, summarize, and organize dissolved metabolome data:

#####__Pull in environmental dissolved metabolite data:
diss.dat <- read_csv(enviro.conc.file) %>%
  mutate(Diss_Samp_ID = SampID) 

exp.enviro.dat <- read_csv(match.file) %>%
  select(-Compound) %>%
 # mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
  left_join(diss.dat) %>%
 # select(KinExp_ID, Diss.Conc.nM) %>%
  group_by(KinExp_ID, Compound) %>%
  reframe(Mean.Diss.Conc.nM = mean(Diss.Conc.nM),
          SD.Diss.Conc.nM = sd(Diss.Conc.nM)) %>%
  separate(KinExp_ID, into=c("Cruise", "exp"))

write_csv(exp.enviro.dat, file = "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv")




###Explore Enviro metab data to summarize for paper:
exp.enviro.dat.tot.sum <- read_csv(match.file) %>%
  select(-Compound) %>%
  # mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
  left_join(diss.dat) %>%
  # select(KinExp_ID, Diss.Conc.nM) %>%
  group_by(Rep, KinExp_ID) %>%
  reframe(tot.conc.nM = sum(Diss.Conc.nM)) %>%
  group_by(KinExp_ID) %>%
  reframe(Mean.Diss.Conc.nM = sum(tot.conc.nM),
          SD.Diss.Conc.nM = sd(tot.conc.nM)) 

exp.enviro.dat.relcontrib.sum <- read_csv(match.file) %>%
  select(-Compound) %>%
  # mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
  left_join(diss.dat) %>%
  # select(KinExp_ID, Diss.Conc.nM) %>%
  group_by(KinExp_ID) %>%
  mutate(tot.conc.nM = sum(Diss.Conc.nM)) %>%
  mutate(group = case_when(Compound %in% c("Dimethylsulfonioacetate",
                                           "Trigonelline", "Betonicine",
                                           "Proline betaine") ~ "Minor",
                           TRUE ~ "Major")) %>%
  group_by(KinExp_ID, group) %>%
  reframe(rel_contrib = mean(sum(Diss.Conc.nM/tot.conc.nM)))

gbt.hom.comp <- exp.enviro.dat %>%
  group_by(Cruise, exp) %>%
  mutate(GBT.hom.ratio = Mean.Diss.Conc.nM[Compound == "Glycine betaine"]/Mean.Diss.Conc.nM[Compound == "Homarine"])








