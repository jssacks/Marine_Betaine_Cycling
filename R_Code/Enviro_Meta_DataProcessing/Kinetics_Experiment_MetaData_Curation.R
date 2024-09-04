


library(tidyverse)
library(lubridate)



### Define inputs:

#Matching Keys and time/location information for experiments
exp.file <- "Meta_Data/Ingalls_Lab_Data/Kin_Exp_MetaData.csv"

#Meta data from KM1906:
g3.file <- "Meta_Data/Enviro_Meta_Data/KM1906/G3_Kinetics_Experiments_Meta_Data.csv"

#Meta data from TN397
g4.file <- "Intermediates/G4_MetaData_Interpolated.csv"

#Meta data from TN412
g5.file <- "Intermediates/G5_MetaData_Interpolated.csv"

#Meta data from RC078
d1.file <- "Intermediates/RC078_MetaData_Compiled.csv"

#Meta data from RC104
d2.file <- "Meta_Data/Enviro_Meta_Data/RC104/RC104_kinexp_metadata.csv"

#Additional Flow Cytometry Measurements
fcm.file <- "Meta_Data/Enviro_Meta_Data/FCM/KinExpFCMData.csv"

#Puget Sound PP Estimates File:
ps.pp.file <- "Meta_Data/Data_From_Other_Studies/PS_14C_PP.csv"


#Load in Experiment Information:
exp.info <- read_csv(exp.file)



###_____FCM data
fcm.dat <- read_csv(fcm.file) %>%
  group_by(KinExp_ID) %>%
  reframe(bact_c = mean(carbon_biomass_bacteria),
          bact_c_sd = sd(carbon_biomass_bacteria))

###PS PP data
ps.pp.est <- read_csv(ps.pp.file) %>%
  rename("pp_14c" = PP_nMC_median) %>%
  select("pp_14c")



#____Match G4 and G5 Experiments to Interpolated Data:

###G4:
#convert G4 experiment times in rounded datetimes for matching with metadata
g4.info <- exp.info %>%
  filter(str_detect(KinExp_ID, "TN397")) %>%
  mutate(utc.date = mdy(Exp_Date_UTC)) %>%
  mutate(time = as_datetime(paste(utc.date, Exp_Time_UTC))) %>%
  mutate(time = round_date(time, unit = "hour"))

#load in metadata:
g4.dat <- read_csv(g4.file)

g4.dat.info <- left_join(g4.info, g4.dat) %>%
  select(KinExp_ID, sst, sss, chla_interp, pc_interp, pn_interp, N_N_interp, bact_c_interp, pp_13c_interp, pp_14c_interp) %>%
  rename(chla = chla_interp,
         pc = pc_interp,
         pn = pn_interp,
         N_N = N_N_interp,
         pp_13c = pp_13c_interp,
         pp_14c = pp_14c_interp,
         bact_c = bact_c_interp) %>%
  mutate(pp_13c = pp_13c*(1/1000)*(1/12.001)*1e6,    #convert 13C and 14C PP measurments to nmol C/L/day
         pp_14c = pp_14c*(1/1000)*(1/12.001)*1e6)


###G5:
#convert G5 experiment times in rounded datetimes for matching with metadata
g5.info <- exp.info %>%
  filter(str_detect(KinExp_ID, "TN412")) %>%
  mutate(utc.date = mdy(Exp_Date_UTC)) %>%
  mutate(time = as_datetime(paste(utc.date, Exp_Time_UTC))) %>%
  mutate(time = round_date(time, unit = "hour"))

#load in metadata:
g5.dat <- read_csv(g5.file)

g5.dat.info <- left_join(g5.info, g5.dat) %>%
  select(KinExp_ID, sst, sss, chla_interp, N_N_interp, pp_13c_interp) %>%
  rename(chla = chla_interp,
      #   pc = pc_interp,
       #  pn = pn_interp,
         pp_13c = pp_13c_interp,
         N_N = N_N_interp) %>%
  mutate(pp_13c = pp_13c*(1/1000)*(1/12.001)*1e6) %>%   #convert 13C and 14C PP measurments to nmol C/L/day
  left_join(., fcm.dat)




###___________ Match RC078 experiments to CTD bottle data 

#pull out RC078 info
d1.info <- exp.info %>%
  filter(str_detect(KinExp_ID, "RC078"))

#load in RC078 data
d1.dat <- read_csv(d1.file) %>%
  rename("RC078_CTD_matching" = parent_id)

###Combine and summarize
d1.dat.info <- left_join(d1.info, d1.dat) %>%
  select(KinExp_ID, POC_uM, PN_uM, NO3, NO2, sal, temp, Chl_fluor) %>%
  group_by(KinExp_ID) %>%
  reframe(pc = mean(POC_uM),
          pn = mean(PN_uM),
          N_N = mean(NO3+NO2),
          sss = mean(sal),
          sst = mean(temp),
          chla = mean(Chl_fluor)) %>%
  left_join(., fcm.dat) %>%
  cross_join(ps.pp.est)


###Pull in data from Boysen et al. 2022 for KM1906 data

#pull out KM1906 info
g3.info <- exp.info %>%
  filter(str_detect(KinExp_ID, "KM1906"))

#load in KM1906 data
g3.dat.means <- read_csv(g3.file) %>%
  pivot_wider(id_cols = c("Cruise", "Experiment"), names_from = Parameter, values_from = Mean) %>%
  rename("pc" = PC_uM,
         "pn" = PN_uM,
         "N_N" = DIN_uM,
         "sst" = Temp_C,
         "sss" = Sal,
         "chla" = Chl_ug_L,
         "bact_c" = Bact_Biomass_ugC_L,
         "pp_13c" = pp_13C_umolC_L_d) %>%
  mutate(pp_13c = pp_13c*1000) %>%      #convert uMC to nMC
  mutate(KinExp_ID = case_when(Experiment == "UKG1" ~ "KM1906_UKG1",
                               TRUE ~ "KM1906_UKG2")) %>%
  select(KinExp_ID, pc, pn, N_N, sst, sss, chla, bact_c, pp_13c)

g3.dat.sd <- read_csv(g3.file) %>%
  pivot_wider(id_cols = c("Cruise", "Experiment"), names_from = Parameter, values_from = SD) %>%
  mutate(KinExp_ID = case_when(Experiment == "UKG1" ~ "KM1906_UKG1",
                               TRUE ~ "KM1906_UKG2")) %>%
  select(KinExp_ID, Bact_Biomass_ugC_L) %>%
  rename("bact_c_sd" = Bact_Biomass_ugC_L)

g3.dat.info <- left_join(g3.dat.means, g3.dat.sd)



### RC104 Metadata:
d2.dat <- read_csv(d2.file) %>%
  cross_join(ps.pp.est)




###_________Combine all data together into a single document:

#########.  Need to combine all of the values for each experiment group individually and then 
#           rbind them together 
j1 <- full_join(g4.dat.info, g5.dat.info)

j2 <- full_join(g4.dat.info, d1.dat.info)

j3 <- full_join(g4.dat.info, g3.dat.info)

j4 <- full_join(g4.dat.info, d2.dat)
  
  
bind1 <- rbind(j1, j2, j3, j4) %>%
  unique()


#Combine all datasets together with experiment informaiton:
kin.dat.info <- exp.info %>%
  left_join(., bind1) 

#export
write_csv(kin.dat.info, file = "Intermediates/Kinetics_Meta_Data_Compiled.csv")















































