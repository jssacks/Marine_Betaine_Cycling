

library(tidyverse)


#define inputs:
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"



#Meta.Data:
meta.dat <- read_csv(meta.data.file)

meta.dat.table <- meta.dat %>%
  select(KinExp_ID, Compound, Lat, Long, Region, Exp_Time_Local, Exp_Date_Local,
         Exp_Time_UTC, Exp_Date_UTC, sst, sss, chla, pc, pn, N_N, bact_c, bact_c_sd, pp_13c, pp_14c) %>%
  separate(KinExp_ID, into = c("Cruise", "Experiment")) %>%
  mutate(sst = print(formatC(signif(sst, digits=3), digits=3,format="fg", flag = "#")),
         sss = print(formatC(signif(sss, digits=3), digits=3,format="fg", flag = "#")),
         chla = print(formatC(signif(chla, digits=3), digits=3,format="fg", flag = "#")),
         pc = print(formatC(signif(pc, digits=3), digits=3,format="fg", flag = "#")),
         pn = print(formatC(signif(pn, digits=3), digits=3,format="fg", flag = "#")),
         N_N = print(formatC(signif(N_N, digits=3), digits=3,format="fg", flag = "#")),
         bact_c = print(formatC(signif(bact_c, digits=3), digits=3,format="fg", flag = "#")),
         bact_c_sd = print(formatC(signif(bact_c_sd, digits=2), digits=2,format="fg", flag = "#")),
         pp_13c = print(formatC(signif(pp_13c, digits=3), digits=3,format="fg", flag = "#")),
         pp_14c = print(formatC(signif(pp_14c, digits=3), digits=3,format="fg", flag = "#"))) %>%
  rename("Temperature (C)" = sst,
         "Salintiy (‰)" = sss,
         "Chlorophyl-a (mg/m^3)" = chla,
         "Particulate Carbon (umol C/L)" = pc,
         "Particulate Nitrogen (umol N/L)" = pn,
         "Nitrate and Nitrite (umol/L)" = N_N,
         "Bacterial Biomass (ug C/L)" = bact_c,
         "Standard Deviation of Bacterial Biomass (ug C/L)" = bact_c_sd,
         "13C Derived Primary Production Estimate (nmol C/L/day)" = pp_13c,
         "14C Derived Primary Productin Estimate (nmol C/L/day)" = pp_14c)

##export:
write_csv(meta.dat.table, file = "Tables/Flux_Tables/Meta_Data_Table.csv")
