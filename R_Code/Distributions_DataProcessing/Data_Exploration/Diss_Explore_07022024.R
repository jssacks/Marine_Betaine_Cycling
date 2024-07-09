


library(tidyverse)

diss.dat <- read_csv("Intermediates/Dissolved_Quantified_Data.csv")
lod.dat <- read_csv("Intermediates/Dissolved_Blk_LOD_Concentrations.csv")


dis.lod.dat <- left_join(diss.dat, lod.dat) %>%
  filter(!str_detect(SampID, "Poo"),
         !str_detect(SampID, "Blk"),
         !str_detect(SampID, "CXBLK"))

dis.blanks <- diss.dat %>%
  filter(str_detect(.$SampID, "Blk")) %>%
  filter(!str_detect(.$SampID, "MQBlk_")) %>%
  filter(!str_detect(.$SampID, "TN397_MQBlk")) %>%
  filter(!str_detect(.$SampID, "FilterBlk")) %>%
  filter(!str_detect(.$SampID, "BottleBlk")) %>%
  filter(!str_detect(.$SampID, "CXC_Blk")) %>%
  mutate(EE.adjust.conc = replace_na(EE.adjust.conc, 0)) %>%    
  group_by(Cruise, Name) %>%
  summarize(Average_Blank = mean(EE.adjust.conc),
            SD_Blank = sd(EE.adjust.conc),
            Max_Blank = max(EE.adjust.conc),
            Min_Blank = min(EE.adjust.conc))


ggplot(dis.lod.dat, aes(x = Name, y = EE.adjust.conc, color = Name)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  facet_wrap(.~Cruise, scales = "free") +
  scale_y_log10()

dis.lod.sum <- dis.lod.dat %>%
  group_by(Cruise, Name) %>%
  summarize(median = median(EE.adjust.conc),
            mean = mean(EE.adjust.conc),
            sd = sd(EE.adjust.conc),
            max = max(EE.adjust.conc),
            min = min(EE.adjust.conc))
