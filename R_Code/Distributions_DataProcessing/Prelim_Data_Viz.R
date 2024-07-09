


library(tidyverse)



p.q.dat <- read_csv("Intermediates/Particulate_Quant_Output.csv")
g3.samp.loc <- read_csv("Meta_Data/Ingalls_Lab_Data/G3_Samp_Locations.csv")
g4.samp.loc <- read_csv("Meta_Data/Ingalls_Lab_Data/G4_Station_Locations.csv")
D1.samp.dat <- read_csv("Meta_Data/Ingalls_Lab_Data/RC078_metadata.csv")


p.q.sum <- p.q.dat %>%
  group_by(Name, Cruise) %>%
  summarize(count = n(),
            max = max(nM.in.smp),
            mean = mean(nM.in.smp),
            median = median(nM.in.smp),
            min = min(nM.in.smp))


ggplot(p.q.dat, aes(x = Name, y = nM.in.smp, color = Name)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(.~Cruise, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_log10(limits = c(0.0001, 100))

ggplot(p.q.dat, aes(x = SampID, y = nM.in.smp, fill = Name)) +
  geom_col(position = "fill") +
  facet_wrap(.~Cruise, scales = "free")

ggplot(p.q.dat, aes(x = SampID, y = nM.in.smp, fill = Name)) +
  geom_col(position = "stack") +
  facet_wrap(.~Cruise, scales = "free")

p.q.rank <- p.q.dat %>%
  group_by(Name) %>%
  mutate(rank.order = rank(nM.in.smp))

ggplot(p.q.rank, aes(x = rank.order, y = nM.in.smp, color = Cruise)) +
  geom_jitter(alpha = 0.2, size = 2) +
  facet_wrap(.~Name, scales = "free_y") +
  scale_y_log10() 



######organize location data:
####Location data
g3.loc <- g3.samp.loc %>%
  mutate(Cruise = "KM1906") %>%
  select(Cruise, Sample_ID, everything()) %>%
  unite(Cruise:Sample_ID, col = "Samp_ID") %>%
  filter(!is.na(Lat))


g4.loc <- g4.samp.loc %>%
  mutate(Lat = Lat_N,
         Long = -1*abs(Long_W)) %>%
  select(Samp_ID, Local_Date, Local_Time, Lat, Long) %>%
  filter(!is.na(Lat))

all.loc <- rbind(g3.loc, g4.loc)

###combine location and metabolite data:
metab.loc.dat <- p.q.dat %>%
  mutate(Samp_ID = str_remove(SampID, "220628_Smp_")) %>%
  mutate(Samp_ID = str_remove(Samp_ID, "220902_Smp_TN397_")) %>%
  mutate(Samp_ID = paste(Cruise, Samp_ID, sep = "_")) %>%
  filter(!Cruise == "RC078") %>%
  left_join(., all.loc)

ggplot(metab.loc.dat, aes(x = Lat, y = nM.in.smp, color = Long)) +
  geom_point(alpha = 0.8, size = 1.5) +
  facet_wrap(.~Name, scales = "free")# +
 # scale_y_log10()
  

g3.dat <- metab.loc.dat %>%
  filter(Cruise == "KM1906")

ggplot(g3.dat, aes(x = Lat, y = nM.in.smp, color = Local_Date)) +
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(.~Name, scales = "free") +
  scale_y_log10()


####
g4.dat <- metab.loc.dat %>%
  filter(Cruise == "TN397") %>%
  filter(Lat < 30)

ggplot(g4.dat, aes(x = Lat, y = nM.in.smp, color = as.numeric(Local_Time))) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth() +
  facet_wrap(.~Name, scales = "free") #+
#  scale_y_log10()




########################



