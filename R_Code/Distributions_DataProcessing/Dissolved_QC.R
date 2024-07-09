#library
library(tidyverse)


#Define inputs
raw.file <- "Intermediates/dissolved_betaine_data_raw.csv"
lod.file <- "Intermediates/Dissolved_BLk_LOD.csv"


#Define QC thresholds
min.area.threshold <- 40000 #minimum peak area threshold 
min_perc_replicates <- 0.5 #minimum percentage of replicates that a compound must pass QC in to be considered presten  


#Load in datasets
raw.dat <- read_csv(raw.file)
lod.dat <- read_csv(lod.file)


####Perform minimum area and blk_lod QC 
qc.dat <- raw.dat %>%
  left_join(., lod.dat %>%
              rename("Compound" = MF)) %>%
  filter(!str_detect(Rep, "Std")) %>%
  filter(!str_detect(Rep, "Blk")) %>%
  filter(!str_detect(Rep, "Poo")) %>%
  mutate(min.area.flag = case_when(Area <= min.area.threshold ~ "Flag",
                                   TRUE ~ NA),
         blk.lod.flag = case_when(Area <= Blk.LD ~ "Flag",
                                    TRUE ~ NA))



#wrangle and tidy sample IDs into sample names and replicates
samp.id.key <- raw.dat %>%
  filter(!str_detect(Rep, "Std")) %>%
  filter(!str_detect(Rep, "Blk")) %>%
  filter(!str_detect(Rep, "Poo")) %>%
  mutate(SampID = str_remove(Rep, "220623_Smp_"),
         SampID = str_remove(SampID, "240102_Smp_"),
         SampID = str_remove(SampID, "240415_Smp_"),
         SampID = str_remove(SampID, "220602_Smp_"),
         SampID = str_remove(SampID, "TN397_"),
         SampID = str_remove(SampID, "U_"),
         SampID = str_remove(SampID, "220602_Smp_TN397_")) %>%
  mutate(replicate = case_when(str_detect(SampID, "_A$") ~ "A",
                               str_detect(SampID, "_B$") ~ "B",
                               str_detect(SampID, "_C$") ~ "C")) %>%
  mutate(SampID = str_remove(SampID, "_A$"),
         SampID = str_remove(SampID, "_B$"),
         SampID = str_remove(SampID, "_C$")) %>%
  select(-Area)



###Perform replicated QC
###Combine replicate information with QC information
dat.qc.rep <- left_join(samp.id.key, qc.dat) %>%
  group_by(Compound, Cruise, SampID) %>%
  mutate(rep.count = n()) %>%
  ungroup() %>%
  filter(min.area.flag == "Flag" | blk.lod.flag == "Flag") %>%
  group_by(Compound, Cruise, SampID) %>%
  mutate(min.area.flag.count = sum(!is.na(min.area.flag))) %>%
  mutate(blk.flag.count = sum(!is.na(blk.lod.flag))) %>%
  mutate(min.area.flag.ratio = min.area.flag.count/rep.count,
         blk.flag.ratio = blk.flag.count/rep.count) %>%
  mutate(smp.remove = case_when(min.area.flag.ratio >= min_perc_replicates | 
                                  blk.flag.ratio >= min_perc_replicates ~ "Smp_Remove",
                                TRUE ~ NA)) %>%
  mutate(blk.imputed.value = Blk.LD/2) %>%
  ungroup()


###Create QCed dataset where samples not passing QC are removed
dat.qc.remove <- samp.id.key %>%
  left_join(., raw.dat) %>%
  left_join(., dat.qc.rep) %>%
  filter(is.na(smp.remove)) %>%
  select(Rep, SampID, replicate, Cruise, Compound, Area, min.area.flag, blk.lod.flag, Blk.Av)


###Create blk_imputed_value dataset where samples not passing QC have values equal to 1/2 the blank imputed
dat.qc.impute <- raw.dat %>%
  left_join(., samp.id.key) %>%
  left_join(., dat.qc.rep) %>%
  mutate(Area = case_when(blk.lod.flag == "Flag" | min.area.flag == "Flag" ~ blk.imputed.value,
                          TRUE ~ Area))  %>%
  select(Rep, SampID, replicate, Cruise, Compound, Area, min.area.flag, blk.lod.flag, smp.remove, Blk.Av)


###Export both QCed datasets:
write_csv(dat.qc.remove, file = "Intermediates/Dissolved_betaine_QCdat_samplesremoved.csv")
write_csv(dat.qc.impute, file = "Intermediates/Dissolved_betaine_QCdat_blanksimputed.csv")
