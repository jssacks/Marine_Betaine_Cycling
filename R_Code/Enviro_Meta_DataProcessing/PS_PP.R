

library(tidyverse)


#define inputs:
pp.file <- "Meta_Data/Data_From_Other_Studies/PS_PP_Compilation_Newton_VanVoorhis_2002.csv"
  
  
  
#read in data and convert c and chl to concentrations
pp.dat <- read_csv(pp.file) %>%
  mutate(PP_mgC_m3_d1 = PP_mg_C_m2_d1/`Euphotic Zone Depth_m`) %>%
  mutate(PP_nmolC_L_d = PP_mgC_m3_d1*(1/12.001)*(1/1000)*(1e6)) %>%
  mutate(Chl_mgChl_m3 = Chl_mg_chl_m2/`Euphotic Zone Depth_m`)

#summarize PP in monthly averages, standard deviations, max-min, and median
pp.month.sum <- pp.dat %>%
  group_by(Month) %>%
#  filter(!PP_nmolC_L_d > 100000) %>%
  reframe(PP_nMC_ave = mean(PP_nmolC_L_d),
          PP_nMC_sd = sd(PP_nmolC_L_d),
          PP_nMC_min = min(PP_nmolC_L_d),
          PP_nMC_max = max(PP_nmolC_L_d),
          PP_nMC_median = median(PP_nmolC_L_d),
          PP_mgC_ave = mean(PP_mgC_m3_d1),
          PP_mgC_sd = sd(PP_mgC_m3_d1),
          PP_mgC_min = min(PP_mgC_m3_d1),
          PP_mgC_max = max(PP_mgC_m3_d1),
          PP_mgC_median = median(PP_mgC_m3_d1))

###
pp.overall.sum <- pp.dat %>%
 # filter(!PP_nmolC_L_d > 100000) %>%
  reframe(PP_nMC_ave = mean(PP_nmolC_L_d),
          PP_nMC_sd = sd(PP_nmolC_L_d),
          PP_nMC_min = min(PP_nmolC_L_d),
          PP_nMC_max = max(PP_nmolC_L_d),
          PP_nMC_median = median(PP_nmolC_L_d),
          PP_mgC_ave = mean(PP_mgC_m3_d1),
          PP_mgC_sd = sd(PP_mgC_m3_d1),
          PP_mgC_min = min(PP_mgC_m3_d1),
          PP_mgC_max = max(PP_mgC_m3_d1),
          PP_mgC_median = median(PP_mgC_m3_d1)) %>%
  mutate(Month = "All")

#Export Overall Summary:
write_csv(pp.overall.sum, file = "Meta_Data/Data_From_Other_Studies/PS_14C_PP.csv")
  