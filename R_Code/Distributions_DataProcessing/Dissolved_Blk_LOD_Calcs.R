


#library
library(tidyverse)
####Quality Control Script

##Define inputs:
smp.file <- "Intermediates/Dissolved_HILIC_Pos_BMISed_dat.csv"



#####Blanks QC
Blank.dat <- read_csv(smp.file) %>%
  filter(str_detect(.$SampID, "Blk")) %>%
  filter(!str_detect(.$SampID, "MQBlk_")) %>%
  filter(!str_detect(.$SampID, "TN397_MQBlk")) %>%
  filter(!str_detect(.$SampID, "FilterBlk")) %>%
  filter(!str_detect(.$SampID, "BottleBlk")) %>%
  filter(!str_detect(.$SampID, "CXC_Blk")) %>%    
  mutate(Area = replace_na(Area, 0)) %>%    
  mutate(Adjusted_Area = replace_na(Area, 0)) 


##Count number of blanks and assign student's t-value
Blk.sum <- Blank.dat %>%
  select(SampID, Cruise) %>%
  unique() %>%
  group_by(Cruise) %>%
  summarize(count = n()) %>%
  mutate(t_val = case_when(count == 12 ~ 1.782,
                           count == 14 ~ 1.761,
                           count == 15 ~ 1.753)) %>%
  ungroup()
  
  
### Calculate AVERAGE BLANK VALUES AND CORRECT FOR DMSP HAVING A 
# BLANK OF 0 in KIM1906 and imput half the value of a detected peak (20000)
Blk.ave.dat <- Blank.dat %>%
  group_by(MF, Cruise) %>%
  left_join(., Blk.sum) %>%
  summarize(Blk.Av = mean(Adjusted_Area),
            Blk.sd = sd(Adjusted_Area),
            Blk.max = max(Adjusted_Area),
            Blk.LD = Blk.Av + (t_val * (Blk.sd/sqrt(count)))) %>%
  unique() %>%    
  ungroup() %>%
  mutate(Blk.LD = case_when(MF == "Dimethylsulfoniopropionate" & Cruise == "KM1906" ~ 20000,
                            TRUE ~ Blk.LD))
       
write_csv(Blk.ave.dat, file = "Intermediates/Dissolved_Blk_LOD.csv")



































