



library(tidyverse)
library(scales)




#define inputs
kf.file <- "Intermediates/Compiled_Kinetics_Fluxes.csv"
f.rel.file <- "Intermediates/flux_contextualizaing_data"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"



###Make Table 1 (Kinetics Parameters)

#Get kinetics data
kf.dat <- read_csv(kf.file)

#Get Region dat:
region.dat <- read_csv(meta.data.file) %>%
  select(KinExp_ID, Region) %>%
  separate(KinExp_ID, into = c("Cruise", "exp"), remove = FALSE) %>%
  mutate(Region = case_when(Region == "Gyre" ~ "NPSG",
                            Region == "Coastal" ~ "PS",
                            Region == "Equatorial" ~ "Eq",
                            TRUE ~ Region)) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "NPSG", "Eq", "PS"))) 


#Get particualte concentrations:
p.conc.dat <- read_csv(f.rel.file) %>%
  select(Cruise, exp, mean.part.conc.nM, sd.part.conc.nM)



kin.table <- kf.dat %>%
  left_join(., p.conc.dat) %>% 
  left_join(., region.dat) %>%
  select(Region, Cruise, exp, Compound, mean.part.conc.nM, sd.part.conc.nM, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, mean_ks, sd_ks, mean_vmax, sd_vamx) %>%
  rename("sd_vmax" = sd_vamx) %>%
  mutate(Mean.Part.Conc.nM = print(formatC(signif(mean.part.conc.nM, digits=3), digits=3,format="fg")),
         SD.Part.Conc.nM = print(formatC(signif(sd.part.conc.nM, digits=2), digits=2,format="fg")),
         Mean.Diss.Conc.nM = print(formatC(signif(Mean.Diss.Conc.nM, digits=3), digits=3,format="fg")),
         SD.Diss.Conc.nM = print(formatC(signif(SD.Diss.Conc.nM, digits=2), digits=2,format="fg")),
         mean_ks = print(formatC(signif(mean_ks, digits=3), digits=3,format="fg")),
         sd_ks = print(formatC(signif(sd_ks, digits=2), digits=2,format="fg")),
         mean_vmax = print(formatC(signif(mean_vmax, digits=3), digits=3,format="fg")),
         sd_vmax = print(formatC(signif(sd_vmax, digits=2), digits=2,format="fg"))) %>%
  mutate(P.Conc.nM = paste(Mean.Part.Conc.nM, SD.Part.Conc.nM, sep = "\u00B1"),
         D.Conc.nM = paste(Mean.Diss.Conc.nM, SD.Diss.Conc.nM, sep = "\u00B1"),
         ks = paste(mean_ks, sd_ks, sep = "\u00B1"),
         vmax = paste(mean_vmax, sd_vmax, sep = "\u00B1")) %>%
  select(Region, Cruise, exp, Compound, P.Conc.nM, D.Conc.nM, ks, vmax) %>%
  rename("Experiment" = exp,
         "Dissolved Concentration (nM)" = D.Conc.nM,
         "Particulate Concentration (nM)" = P.Conc.nM) %>%
  rename("Kt (nM)" = ks,
         "Vmax nmol L^-1 day^-1" = vmax) %>%
  arrange(Compound)

#export
write_csv(kin.table, file = "Tables/Flux_Tables/Main_Text_Table1_Kinetics.csv")




###TT and FLux Table (Main Text Table 2)
flux.table <- kf.dat %>%
  left_join(., region.dat) %>%
  select(Region, Cruise, exp, Compound, mm_tt, mm_tt_sd, mm_flux_nM_day, mm_flux_sd_nM_day) %>%
  mutate(carbon = case_when(Compound == "GBT" ~ 5,
                            Compound == "Homarine" ~ 7),
         C_flux_nM_day = mm_flux_nM_day*carbon,
         C_flux_nM_day_sd = mm_flux_sd_nM_day*carbon) %>%
  mutate(TT = print(formatC(signif(mm_tt, digits=3), digits=3,format="fg")),
         TT_sd = print(formatC(signif(mm_tt_sd, digits=2), digits=2,format="fg")),
         flux = print(formatC(signif(mm_flux_nM_day, digits=3), digits=3,format="fg")),
         flux_sd = print(formatC(signif(mm_flux_sd_nM_day, digits=2), digits=2,format="fg")),
         c_flux = print(formatC(signif(C_flux_nM_day, digits=3), digits=3,format="fg")),
         c_flux_sd = print(formatC(signif(C_flux_nM_day_sd, digits=2), digits=2,format="fg"))) %>%
  mutate(TT = paste(TT, TT_sd, sep = "\u00B1"),
         Flux = paste(flux, flux_sd, sep = "\u00B1"),
         C_Flux = paste(c_flux, c_flux_sd, sep = "\u00B1")) %>%
  select(Region, Cruise, exp, Compound, TT, Flux, C_Flux) %>%
  rename("Experiment" = exp,
         "Turnover Time (hr)" = TT,
         "Flux (nmol compound L^-1 day^-1)" = Flux,
         "Carbon Flux (nmol C L^-1 day^-1)" = C_Flux) %>%
  arrange(Compound)

#export
write_csv(flux.table, file = "Tables/Flux_Tables/Main_Text_Table2_Fluxes.csv")






###Make Wright-Hobbie TT Table (Supplemental):
WH.table <- kf.dat %>%
  select(Cruise, exp, Compound, wh_mean_tt, wh_sd_tt) %>%
  mutate(TT = print(formatC(signif(wh_mean_tt, digits=3), digits=3,format="fg")),
         TT_sd = print(formatC(signif(wh_sd_tt, digits=2), digits=2,format="fg"))) %>%
  select(Cruise, exp, Compound, TT, TT_sd) %>%
  rename("Experiment" = exp,
         "WH TT (hr)" = TT,
         "WH TT StdDev (hr)" = TT_sd) %>%
  arrange(Compound) 

#export
write_csv(WH.table, file = "Tables/Flux_Tables/WH_TT_Supplemental_Table.csv")






###Make Flux in Context Table (Percent of Primary Production, etc.): 
f.rel.table <- read_csv(f.rel.file) %>%
  select(Cruise, exp, Compound, mean.part.conc.nM, sd.part.conc.nM, 
         Flux_Perc_of_Part, Flux_Perc_of_Part_Error, pp.uM.C.perday, flux_Perc_of_PP, flux_Perc_of_PP_Error) %>%
  mutate(mean.part.conc.nM = print(formatC(signif(mean.part.conc.nM, digits=3), digits=3,format="fg")),
         sd.part.conc.nM = print(formatC(signif(sd.part.conc.nM, digits=2), digits=2,format="fg")),
         Flux_Perc_of_Part = print(formatC(signif(Flux_Perc_of_Part, digits=3), digits=3,format="fg")),
         Flux_Perc_of_Part_Error = print(formatC(signif(Flux_Perc_of_Part_Error, digits=2), digits=2,format="fg")),
         pp.uM.C.perday = print(formatC(signif(pp.uM.C.perday, digits=3), digits=3,format="fg")),
         flux_Perc_of_PP = print(formatC(signif(flux_Perc_of_PP, digits=3), digits=3,format="fg")),
         flux_Perc_of_PP_Error = print(formatC(signif(flux_Perc_of_PP_Error, digits=2), digits=2,format="fg",))) %>%
  rename("Experiment" = exp,
         "Average Particulate Concentration (nM)" = mean.part.conc.nM,
         "Standard Deviation of Particulate Concentration (nM)" = sd.part.conc.nM, 
         "Daily Flux as Percentage of Particulate Pool (%)" = Flux_Perc_of_Part,
         "Daily Flux as Percentage of Particulate Pool Error (%)" = Flux_Perc_of_Part_Error,
         "Estimated Primary Production (uM C day^-1)" = pp.uM.C.perday,
         "Daily Flux as Percentage of Primary Production (%)" = flux_Perc_of_PP,
         "Daily Flux as Percentage of Primary Production Error (%)" = flux_Perc_of_PP_Error) %>%
  arrange(Compound)

#export:
write_csv(f.rel.table, file = "Tables/Flux_Tables/Flux_in_Context_Supplemental_Table.csv")
  





