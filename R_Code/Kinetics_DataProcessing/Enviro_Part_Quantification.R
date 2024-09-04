


###Load Packages and Functions
library(tidyverse)
library(lubridate)
library(broom)

source("R_Code/Kinetics_DataProcessing/Functions.R")

###define inputs
cruises <- c("TN397", "RC078", "TN412", "RC104")


###Read in raw data using TQS_csv_read() function
dat.1 <- TQS_csv_read("Raw_Data/Kinetics_Data/G4_UKH4_UCH1_TQS_results.csv") #data from TN397 UKH4 and UCH1
dat.2 <- TQS_csv_read("Raw_Data/Kinetics_Data/G4_UKH5_TQS_results.csv") #data from TN397 UKH5 
dat.3 <- TQS_csv_read("Raw_Data/Kinetics_Data/G4_Kinetics_TQS_UKH123.csv") #data from TN397 UKH1, UKH2, and UKH3
dat.4 <- TQS_csv_read("Raw_Data/Kinetics_Data/TQS_HomarineFocus_SKylineReport.csv") #data from RC078 UKH1
dat.5 <- TQS_csv_read("Raw_Data/Kinetics_Data/KinSept23_GBT_output.csv") #data from G5 (TN412) UKG, RC078 UKG1, RC104 UKG1
dat.6 <- TQS_csv_read("Raw_Data/Kinetics_Data/KinSept23_Hom_output.csv") #data from RC078 UKH2, RC104 UCH1, RC104 UKH1

###Read in metadata 
mdat.raw.1 <- read_csv("Raw_Data/Kinetics_Data/Kin_metadata.csv") #metadata for TN397 
mdat.raw.2 <- read_csv("Raw_Data/Kinetics_Data/Kin_MetaData_2.csv") # metadata for TN412, RC104, and RC078





####_______ 1. Wrangle Sample Data____________________________________________

#Wrangle Sample data for Homarine experiments
dat.smp.hom <- rbind(dat.1, dat.2, dat.3, dat.4, dat.6) %>%
  filter(Compound %in% c("D3-Homarine", "Homarine")) %>% #filter for just labeled homarine 
  mutate(SampID = str_replace(SampID, "808_U", "808_Smp_TN397_U")) %>% # add in cruise info 
  mutate(SampID = str_replace(SampID, "403_Smp_U", "403_Smp_RC078_U")) %>% # add in cruise info 
  mutate(SampID = str_replace(SampID, "_G5_", "_TN412_")) %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(cruise %in% cruises) %>%
  mutate(SampID = str_replace(SampID, "OnM", "0nM")) # fix samples where Os where typed instead of 0s 

#Wrangle Sample data for GBT experiments
dat.smp.gbt <- dat.5  %>%
  filter(Compound %in% c("13C-15N Glycine betaine", "Glycine betaine")) %>% #filter for just labeled GBT
  mutate(SampID = str_replace(SampID, "_G5_", "_TN412_")) %>%
  mutate(SampID = str_remove(SampID, "-10x")) %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(cruise %in% cruises)

####Combine Homarine and GBT sample data
dat.smp <- rbind(dat.smp.hom, dat.smp.gbt)






####________ 2. Wrangle Standard Curve Data___________________________________________________

###Homarine Samples
dat.std.hom <- rbind(dat.1, dat.2, dat.3, dat.4, dat.6) %>%
  filter(str_detect(SampID, "Std")) %>%
  filter(!str_detect(SampID, "In")) %>%
  filter(!str_detect(SampID, "inM")) %>%
  filter(!str_detect(SampID, "Lab")) %>%
  mutate(SampID = str_replace(SampID, "UHC1", "UCH1")) %>%
  mutate(SampID = str_replace(SampID, "ine_UCH1", "ine_TN397_UCH1")) %>% # add in cruise info for TN397 UCH1
  mutate(SampID = str_replace(SampID, "ine_UKH", "ine_TN397_UKH")) %>% # add in cruise info for TN397 UKH1 - UKH5
  mutate(SampID = str_replace(SampID, "Field1", "RC078_UKH1_1")) %>% # add in cruise info for RC078 UKH1
  mutate(SampID = str_replace(SampID, "Field2", "RC078_UKH1_2")) %>% # add in cruise infor RC078 UKH1
  separate(col = SampID, into = c("date", "type", "spike_conc", "hom", "cruise", "exp", "cal_curve"), remove = FALSE) %>%
  filter(Compound == "D3-Homarine") %>%
  mutate(spike_conc = as.numeric(str_remove(spike_conc, "nM")))

###GBT Samples
dat.std.gbt <- dat.5 %>%
  filter(str_detect(SampID, "Std")) %>%
  filter(!str_detect(SampID, "In")) %>%
  mutate(SampID = str_replace(SampID, "_G5_", "_TN412_")) %>%
  separate(col = SampID, into = c("date", "type", "spike_conc", "hom", "cruise", "exp", "cal_curve"), remove = FALSE) %>%
  filter(Compound == "13C-15N Glycine betaine") %>%
  mutate(spike_conc = as.numeric(str_remove(spike_conc, "nM")))

#combine homarine and GBT standard curve data
dat.std <- rbind(dat.std.hom, dat.std.gbt)






#####______ 3. Calculate slopes, intercepts, and errors for standard curves using linear models______________

#####group standard curve data by cruise, experiment, and fragment
dat.cal.curve <- dat.std %>%
  ungroup() %>%
  group_by(cruise, exp, Compound, Fragment_mz) 

##Visualize calibration curves:
ggplot(dat = dat.cal.curve, aes(x = spike_conc, y = Area, color = as.factor(Fragment_mz))) +
  geom_point() +
  facet_wrap(cruise ~ exp, scales = "free") +
  geom_smooth(method = "lm")


###perform linear models across groupings
lm.out <- do(dat.cal.curve,
             tidy(
               lm(Area ~ spike_conc, data = .)))

###pull out, tidy, and rename values for standard curve slopes
lm.out.slopes <- lm.out %>%
  filter(term == "spike_conc") %>%
  select(cruise, exp, Compound, Fragment_mz, estimate, std.error) %>%
  rename("slope" = estimate,
         "std.error.slope" = std.error)

###pull out, tidy, and rename values for standard curve intercepts
lm.out.intercepts  <- lm.out %>%
  filter(term == "(Intercept)") %>%
  select(cruise, exp, Fragment_mz, estimate, std.error) %>%
  rename("intercept" = estimate,
         "std.error.intercept" = std.error)

#combine all calibration curve outputs but calculate standard error as just the standard error of the slopes 
cal.curve.output <- left_join(lm.out.slopes, lm.out.intercepts) %>%
  mutate(std_error = std.error.slope/slope) %>%
  ungroup() %>%
  select(-Compound)






####___________ 4. Wrangle Metadata___________________ 

#Convert time to date-time objects and combine metadata files
mdat.raw.1.dt <- mdat.raw.1 %>%
  select(Exp, Treatment, Samp_ID, Vol_Filt_L, Local_Date, Spike_Time, 
         Filt_time_start, Filt_time_end) %>%
  mutate(Cruise = "TN397") %>% #add in cruise identifier 
  unite("Spike_time", c(Local_Date, Spike_Time), sep = " ", remove = FALSE) %>%
  unite("Filt_time_start", c(Local_Date, Filt_time_start), sep = " ", remove = FALSE) %>%
  unite("Filt_time_end", c(Local_Date, Filt_time_end), sep = " ", remove = FALSE) %>%
  mutate("Spike_time" = dmy_hms(Spike_time),
         "Filt_time_start" =  dmy_hms(Filt_time_start), 
         "Filt_time_end" = dmy_hms(Filt_time_end)) %>%
  select(-Spike_Time)

mdat.raw.2.dt <- mdat.raw.2 %>%
  rename("treatment" = spike_conc)%>%
  select(cruise, exp, treatment, Samp_ID, Vol_Filt_L, Local_Date, Spike_time, 
         Filt_time_start, Filt_time_end) %>% 
  mutate(Samp_ID = str_replace(Samp_ID, "_G5_", "_TN412_")) %>%
  unite("Spike_time", c(Local_Date, Spike_time), sep = " ", remove = FALSE) %>%
  unite("Filt_time_start", c(Local_Date, Filt_time_start), sep = " ", remove = FALSE) %>%
  unite("Filt_time_end", c(Local_Date, Filt_time_end), sep = " ", remove = FALSE) %>%
  mutate("Spike_time" = mdy_hms(Spike_time),
         "Filt_time_start" =  mdy_hms(Filt_time_start), 
         "Filt_time_end" = mdy_hms(Filt_time_end)) %>%
  rename("Cruise" = cruise,
         "Exp" = exp,
         "Treatment" = treatment)

mdat.raw.vol <- rbind(mdat.raw.1.dt, mdat.raw.2.dt) %>%
  select(Cruise, Exp, Treatment, Samp_ID, Vol_Filt_L) %>%
  rename("SampID" = Samp_ID) %>%
  mutate(SampID = str_replace(SampID, "G5_", "TN412_"))




#####___________________  Quantification__________________________________

####Adjust peak areas for dilutions (RC104 and RC078 were diluted 10x)
dat.smp.diladj <- dat.smp %>%
  mutate(dilution.factor = case_when(str_detect(SampID, "RC104_UKG1") ~ 10,
                                     str_detect(SampID, "RC078_UKG1") ~ 10,
                                     TRUE ~ 1
  )) %>%
  mutate(Area.new = Area*dilution.factor)



#####_____Quantify Compounds in Experiments using standard curves and adjusting for blanks

#match fragments from internal standards and standard curves
dat.q <- dat.smp.diladj %>% 
  mutate(Fragment_mz = case_when(        #match unlabled fragments to labeled fragments 
    Fragment_mz == 94.1067 ~ 97.0954,
    Fragment_mz == 78.1109 ~ 78.0399,
    Fragment_mz == 58.9328 ~ 62.8977,
    Fragment_mz == 42.2168 ~ 45.1882,
    TRUE ~ Fragment_mz
  )) %>%
  mutate(SampID = str_remove(SampID, "230808_Smp_"),
         SampID = str_remove(SampID, "230403_Smp_"),
         SampID = str_remove(SampID, "230605_Smp_"),
         SampID = str_remove(SampID, "230607_Smp_"),
         SampID = str_remove(SampID, "230921_Smp_"),
         SampID = str_replace(SampID, "RC104_UKG1_ 0nM", "RC104_UKG1_0nM")) %>%
  select(SampID, cruise, exp, treatment, rep, Area, Compound, Fragment_mz) %>%
  left_join(., cal.curve.output) %>%
  left_join(., mdat.raw.vol) %>%
  mutate(nM_in_vial = Area/slope,
         nM.in.smp = nM_in_vial*(400*1e-6)/Vol_Filt_L)

######Summarize environmental concentrations of unlabeled GBT and Homarine for each experiment
dat.q.sum <- dat.q %>%
  group_by(cruise, exp, Compound, Fragment_mz) %>%
  filter(!is.na(nM.in.smp)) %>%
 # filter(Fragment_mz %in% c(97.0954, 62.8977)) %>%
  filter(!Compound %in% c("13C-15N Glycine betaine", "D3-Homarine")) %>%
  reframe(mean.enviro.nM = mean(nM.in.smp),
          sd.environ.nM = sd(nM.in.smp))


#####Export to csv:
write_csv(dat.q.sum, file = "Intermediates/Uptake_Experiments_Particulate_Enviromental_Concentrations.csv")

