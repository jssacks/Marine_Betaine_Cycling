


###Load Packages and Functions
library(tidyverse)
library(lubridate)
library(broom)

source("R_Code/Kinetics/Functions.R")


###define inputs





###defne analysis thresholds
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





####Wrangle Sample Data____________________________________________

#Wrangle Sample data for Homarine experiments
dat.smp.hom <- rbind(dat.1, dat.2, dat.3, dat.4, dat.6) %>%
  filter(Compound == "D3-Homarine") %>% #filter for just labeled homarine 
  mutate(SampID = str_replace(SampID, "808_U", "808_Smp_TN397_U")) %>% # add in cruise info 
  mutate(SampID = str_replace(SampID, "403_Smp_U", "403_Smp_RC078_U")) %>% # add in cruise info 
  mutate(SampID = str_replace(SampID, "_G5_", "_TN412_")) %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(cruise %in% cruises) %>%
  mutate(SampID = str_replace(SampID, "OnM", "0nM")) # fix samples where Os where typed instead of 0s 

#Wrangle Sample data for GBT experiments
dat.smp.gbt <- dat.5  %>%
  filter(Compound == "13C-15N Glycine betaine") %>% #filter for just labeled GBT
  mutate(SampID = str_replace(SampID, "_G5_", "_TN412_")) %>%
  mutate(SampID = str_remove(SampID, "-10x")) %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(cruise %in% cruises)

####Combine Homarine and GBT sample data
dat.smp <- rbind(dat.smp.hom, dat.smp.gbt)




####Wrangle Standard Curve Data___________________________________________________

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


#####Calculate slopes, intercepts, and errors for standard curves using linear models______________

#####group standard curve data by cruise, experiment, and fragment
dat.cal.curve <- dat.std %>%
  ungroup() %>%
  group_by(cruise, exp, Fragment_mz) 

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
  select(cruise, exp, Fragment_mz, estimate, std.error) %>%
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
  mutate(std_error = std.error.slope/slope)


####Wrangle Metadata 

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

mdat.raw <- rbind(mdat.raw.1.dt, mdat.raw.2.dt)

#Calculate incubation times from metadata to determine rates
mdat.tidy <- mdat.raw %>%
  select(Cruise, Exp, Treatment, Samp_ID, Vol_Filt_L, Spike_time, Filt_time_start, Filt_time_end) %>%
  mutate(exp_time = as.numeric(Filt_time_start - Spike_time),
         filt_time = as.numeric(Filt_time_end - Filt_time_start),
         incubation_time1 = exp_time + (filt_time/2),
         incubation_time2 = exp_time + filt_time) %>%
  rename("SampID" = Samp_ID) %>%
  mutate(SampID = str_replace(SampID, "G5_", "TN412_"))
  
####################


#####Quantification__________________________________

####Adjust peak areas for dilutions (RC104 and RC0___ were diluted 10x)
dat.smp.diladj <- dat.smp %>%
  mutate(dilution.factor = case_when(str_detect(SampID, "RC104_UKG1") ~ 10,
                                     str_detect(SampID, "RC078_UKG1") ~ 10,
                                     TRUE ~ 1
                                     )) %>%
  mutate(Area = Area*dilution.factor)



#####_____Quantify Compounds in Experiments using standard curves and adjusting for blanks
dat.q.1 <- dat.smp.diladj%>%
  select(SampID, cruise, exp, treatment, rep, Area, Fragment_mz) %>%
#  left_join(., m.dat.3) %>%
  left_join(., cal.curve.output) %>%
  mutate(nM_in_vial = Area/slope)

#####Need to adjust for dilutions

####Determine blank relationships to concentration for data sets 
dat.blks <- dat.q.1 %>%
  filter(str_detect(SampID, "Blk"))  %>%
  mutate(treatment = as.numeric(str_remove(treatment, "nM"))) 

###visualize blanks in relation to calibration curve
blk.plot.1 <- ggplot(dat.blks, aes(x = treatment, y = nM_in_vial)) +
  geom_smooth(method = "lm", alpha = 0.5) +
  geom_point(aes(color = cruise)) +
  facet_wrap(.~Fragment_mz, scales = "free_y")
blk.plot.1


dat.blks.fix <- dat.blks %>%
  filter(SampID %in% c("230921_Smp_RC104_UKH1_3000nM_Blk", "230921_Smp_RC104_UKH1_2000nM_Blk", "230921_Smp_RC104_UKH1_1000nM_Blk")) %>%
  mutate(new_SampID = case_when(SampID == "230921_Smp_RC104_UKH1_3000nM_Blk" ~ "230921_Smp_RC104_UKH1_1000nM_Blk",
                                SampID == "230921_Smp_RC104_UKH1_1000nM_Blk" ~ "230921_Smp_RC104_UKH1_2000nM_Blk",
                                SampID == "230921_Smp_RC104_UKH1_2000nM_Blk" ~ "230921_Smp_RC104_UKH1_3000nM_Blk")) %>%
  mutate(new_treatment = case_when(treatment == 3000 ~ 1000,
                                treatment == 1000 ~ 2000,
                                treatment == 2000 ~ 3000)) %>%
  select(new_SampID, new_treatment, cruise, exp, rep, Area, Fragment_mz, slope, std.error.slope, intercept, std.error.intercept, std_error, nM_in_vial) %>%
  rename("SampID" = new_SampID,
         "treatment" = new_treatment) 

dat.blks.final <- dat.blks %>%
  filter(!SampID %in% c("230921_Smp_RC104_UKH1_3000nM_Blk", "230921_Smp_RC104_UKH1_2000nM_Blk", "230921_Smp_RC104_UKH1_1000nM_Blk")) %>%
  rbind(., dat.blks.fix)

###visualize fixed blank data in relation to calibration curve
blk.plot.2 <- ggplot(dat.blks.final, aes(x = treatment, y = nM_in_vial)) +
  geom_smooth(method = "lm", alpha = 0.5) +
  geom_point(aes(color = cruise)) +
  facet_wrap(.~Fragment_mz, scales = "free_y")
blk.plot.2

###visualize fixed blank data in relation to calibration curve but with a log10 scaling on both x and y axes
blk.plot.2.log10 <- ggplot(dat.blks.final, aes(x = treatment, y = nM_in_vial)) +
  geom_smooth(method = "lm", alpha = 0.5) +
  geom_point(aes(color = cruise)) +
  facet_wrap(.~Fragment_mz, scales = "free") +
  scale_x_log10() +
  scale_y_log10()
blk.plot.2.log10


#######Determine blank correction slopes

#group data by Fragment_mz
dat.blk.slope <- dat.blks.final %>%
  ungroup() %>%
  group_by(Fragment_mz) 

###perform linear models across groupings
lm.blk.out <- do(dat.blk.slope,
             tidy(
               lm(nM_in_vial ~ treatment, data = .)))

###pull out, tidy, and rename values for blk curve slopes
lm.blk.slopes <- lm.blk.out %>%
  filter(term == "treatment") %>%
  select(Fragment_mz, estimate, std.error) %>%
  rename("slope" = estimate,
         "std.error.slope" = std.error)

###pull out, tidy, and rename values for standard curve intercepts
lm.blk.intercepts  <- lm.blk.out %>%
  filter(term == "(Intercept)") %>%
  select(Fragment_mz, estimate, std.error) %>%
  rename("intercept" = estimate,
         "std.error.intercept" = std.error)

#combine all calibration curve outputs
blk.lm.output <- left_join(lm.blk.slopes, lm.blk.intercepts)

blk.lm.slope.only <- blk.lm.output %>%
  select(Fragment_mz, slope) %>%
  rename("blk.slope" = slope)






########
####Adjust quantification using blanks to get just microbial uptake values 
dat.q.2 <- dat.smp.diladj%>%
  select(SampID, cruise, exp, treatment, rep, Area, Fragment_mz) %>%
  left_join(., blk.lm.slope.only) %>%
  left_join(., cal.curve.output) %>%
  mutate(treatment = str_replace(treatment, "OnM", "0nM")) %>%
  mutate(treatment_conc = case_when(str_detect(treatment, "nM") ~ as.numeric(str_remove(treatment, "nM")),
                                    TRUE ~ 50)) %>%
  mutate(nM.in.vial = Area/slope,
         blk.conc = treatment_conc*blk.slope,
         nM.in.vial.adj = nM.in.vial-blk.conc,
         se.nM.in.vial.adj = nM.in.vial*std_error)
###
ggplot(dat.q.2, aes(x = treatment_conc, y = nM.in.vial.adj)) +
  geom_point(aes(color = rep)) +
  geom_errorbar(aes(ymin = nM.in.vial.adj-se.nM.in.vial.adj, ymax= nM.in.vial.adj+se.nM.in.vial.adj)) +
  facet_wrap(cruise~exp, scales = "free")



####______________________Calculate Uptake Rates____________________
dat.q.m <- dat.q.2 %>%
  mutate(SampID = str_remove(SampID, "230808_Smp_"),
         SampID = str_remove(SampID, "230403_Smp_"),
         SampID = str_remove(SampID, "230605_Smp_"),
         SampID = str_remove(SampID, "230607_Smp_"),
         SampID = str_remove(SampID, "230921_Smp_"),
         SampID = str_replace(SampID, "RC104_UKG1_ 0nM", "RC104_UKG1_0nM")) %>%
  left_join(., mdat.tidy) %>%
  mutate(nM.in.smp = nM.in.vial.adj*(400*1e-6)/Vol_Filt_L,
         se.nM.in.smp = std_error*nM.in.smp,
         nM.per.hour.1 = nM.in.smp/(incubation_time1/60),
         nM.per.hour.2 = nM.in.smp/(incubation_time2/60),
         se.nM.per.hour.1 = se.nM.in.smp/(incubation_time1/60),
         se.nM.per.hour.2 = se.nM.in.smp/(incubation_time2/60))
  

###Export to csv
write_csv(dat.q.m, file = "Intermediates/quantified_uptake_rates.csv")












####Use standard curves to calculate concentrations________________________

#visualize standard curves
ggplot(dat.std, aes(x = spike_conc, y = Area, color = as.factor(Fragment_mz))) +
  facet_wrap(exp ~ cruise, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point() 




































