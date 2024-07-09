
library(tidyverse)
library(viridis)


g2.dat <- read_csv("Intermediates/G2_SizeFractionated_HILIC_Pos_BMISed_dat.csv")
g2.metadat <- read_csv("Meta_Data/Ingalls_Lab_Data/G2_SF_Sample_MetaData.csv")


g2.blks <- read_csv("Intermediates/g2_betaine_data_raw.csv") %>%
  filter(str_detect(Rep, "Blk")) %>%
  mutate(Blk_type = case_when(str_detect(Rep, "MQ") ~ "MQ",
                              TRUE ~ "FSW")) %>%
  group_by(Compound, Blk_type) %>%
  summarize(ave.blk = mean(Area)) %>%
  rename("MF" = Compound)






#Summarizing, normalizing to vol filtered, and calculating percentage of signal in small size fraction
g2.all.dat <- left_join(g2.dat, g2.metadat) %>%
  filter(!is.na(Station)) %>%
  mutate(Vol.norm.area = Adjusted_Area/Vol_L) %>%
  group_by(MF, Station, Lat, Size_Fraction) %>%
  summarize(Mean.Area = mean(Vol.norm.area)) %>%
  group_by(MF, Station) %>%
  mutate(Sml_Frac_Perc = Mean.Area[Size_Fraction == "S"]/Mean.Area[Size_Fraction == "F"])

###Combining ___ and ___ data to determine if compounds were detected or not:
g2.blk.qc <- full_join(g2.dat, g2.blks) %>%
  left_join(., g2.metadat) %>%
  ungroup() %>%
  filter(!is.na(Station)) %>%
  mutate(blk.ratio = Adjusted_Area/ave.blk)












#heatmap showing % of signal in small size fraction 
ggplot(g2.all.dat, aes(x = as.factor(Lat), y = MF, fill = Sml_Frac_Perc)) +
  geom_tile() +
  scale_fill_viridis() +
  xlab("Latitude")

#bar chart showing distributions
ggplot(g2.all.dat, aes(x = as.factor(Lat), y = Mean.Area, fill = Size_Fraction)) +
  geom_col(position = "dodge") +
  facet_wrap(.~MF, scales = "free") +
  xlab("Latitude") +
  ylab("Peak Area")









