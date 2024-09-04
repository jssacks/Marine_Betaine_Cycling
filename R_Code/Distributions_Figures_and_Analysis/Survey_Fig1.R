


#load packages:
library(tidyverse)
library(viridis)


#Define inputs:

#particulate data
p.file <- "Intermediates/Particualte_Final_Quant_QCed.csv"

#dissolved data
d.file <- "Intermediates/Dissolved_Final_Quant_QCed.csv"

#location data
g3.loc.file <- "Meta_Data/Ingalls_Lab_Data/G3_Samp_Locations.csv"
g4.loc.file <- "Meta_Data/Ingalls_Lab_Data/G4_Samp_Locations.csv"
  
#dissolved and particulate sample matching key:
key.file <- "Intermediates/Distributions_Part_Diss_Sample_Key.csv"



###Organize location data:
g3.loc <- read_csv(g3.loc.file) %>%
  mutate(Cruise = "KM1906") %>%
  rename("SampID" = Sample_ID) %>%
  mutate(SampID = str_remove(SampID, "_A"),
         SampID = str_remove(SampID, "_B"),
         SampID = str_remove(SampID, "_C")
         )

g4.loc <- read_csv(g4.loc.file) %>%
  mutate(Cruise = "TN397",
         Long = -1*abs(Long_W)) %>%
  rename("Lat" = Lat_N,
         "SampID" = Samp_ID) %>%
  select(SampID, Cruise, Lat, Long, Local_Date, Local_Time) %>%
  mutate(SampID = str_remove(SampID, "_A$"),
         SampID = str_remove(SampID, "_B$"),
         SampID = str_remove(SampID, "_C$"),
         SampID = str_remove(SampID, "U_"),
         SampID = str_remove(SampID, "TN397_"),
         SampID = str_remove(SampID, "_U$")
  )





all.samp.loc <- rbind(g3.loc, g4.loc) %>%
  filter(!is.na(SampID)) %>%
  ungroup() 


###Define gradients regions (RC078 are all "Puget_Sound":
gradients.loc <- all.samp.loc %>%
  mutate(Region = case_when(Lat > 20 & Long > -120 ~ "California_Current", 
                            Lat < 35 & Lat > 10 & Long < -120 ~ "Gyre",
                            Lat < 10 & Long < -130 ~ "Equatorial",
                            Lat > 35 & Long < -130 ~ "NPTZ" )) %>%
  unique()







#Pull in Environmental Datasets:

#particulate + dissolved
p.dat <- read_csv(p.file) %>%
  select(Rep, SampID, replicate, Cruise, Compound, nM.in.smp, nM_C, nM_N) %>%
  rename("Part.Rep" = Rep,
         "Part.Conc.nM" = nM.in.smp,
         "Part.C.nM" = nM_C,
         "Part.N.nM" = nM_N) 
  
d.dat <- read_csv(d.file) %>%
  select(Rep, SampID, replicate, Cruise, Compound, Diss.Conc.nM, Diss.Nmol.C, Diss.Nmol.N, LOD.nM, LOD.nM.blk.sub) %>%
  rename("Diss.Rep" = Rep,
         "Diss.C.nM" = Diss.Nmol.C,
         "Diss.N.nM" = Diss.Nmol.N) %>%
  mutate(remove = case_when(Compound == "Glycine betaine" & Cruise == "KM1906" ~ "yes",
                             TRUE ~ "no")) %>%
  filter(remove == "no") %>%
  select(-remove)
  

#pull in matching key:
key <- read_csv(key.file)

#Co-locate Dissolved and Particulate Data and add in region + location information + remove RC078 samples not from surface:
all.dat <- p.dat %>%
  full_join(., key) %>%
  full_join(., d.dat) %>%
  left_join(., gradients.loc) %>%
  filter(!str_detect(SampID, "MQ")) %>% 
  filter(!is.na(Compound)) %>%
  mutate(Region = case_when(Cruise == "RC078" ~ "Puget_Sound",
                            TRUE ~ Region)) %>%
  mutate(Total_nM = Diss.Conc.nM+Part.Conc.nM,
         Perc_Diss = Diss.Conc.nM/(Diss.Conc.nM+Part.Conc.nM)) %>%
  mutate(remove = case_when(Cruise == "RC078" & is.na(Diss.Conc.nM) ~ "yes",
                            TRUE ~ "no")) %>%
  filter(remove == "no") %>%
  select(-remove)

####Create boxplots:

tot.dat <- all.dat %>%
  filter(!is.na(Region)) %>%
  select(-Region) %>%
  mutate(Region = "All")

fig.dat <- all.dat %>%
  filter(!is.na(Region)) %>%
  rbind(., tot.dat) %>%
  mutate(Class = as.factor(case_when(Compound %in% c("beta-Alaninebetaine", "Glycine betaine", "Homarine",
                                           "Trigonelline", "Proline betaine", "Betonicine") ~ "Betaine",
                           Compound %in% c("Dimethylsulfoniopropionate", "Dimethylsulfonioacetate", "Gonyol") ~ "Sulfonium",
                           TRUE ~ "Amine_Oxide"))) %>%
  mutate(Class = fct_relevel(Class, c("Betaine", "Sulfonium", "Amine_Oxide"))) %>%
  mutate(Region = fct_relevel(Region, c("All", "Puget_Sound", "California_Current", "Gyre", "Equatorial", "NPTZ"))) %>%
  group_by(Compound) %>%
  mutate(mean.conc = mean(Part.Conc.nM, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Class) %>%
  mutate(rank = rank(mean.conc)) %>%
  ungroup() %>%
  mutate(remove = case_when(Compound == "Dimethylsulfoniopropionate" & Diss.Conc.nM < 1E-3 ~ "yes",
                             TRUE ~ "no")) %>%  
  filter(remove == "no") %>%
  select(-remove)



#Define color scheme:

region.pal <- c("#E1DFD0", "#6794a7", "#76c0c1", "#014d64", "#01a2d9", "#00887d")



#particulate_boxplot 
p.boxplot <- ggplot(fig.dat, aes(x = reorder(Compound, -rank), y = Part.Conc.nM, fill = Region)) +
  geom_boxplot(aes(color = Region), alpha = 0.5,  position = position_dodge(preserve = "single")) +
  geom_point(aes(fill = Region, group = Region),  alpha = 0.7, shape = 21, stroke = 0.22, size = 2, 
               position=position_jitterdodge(jitter.width = 0.15)) +
  scale_y_log10() +
  facet_grid(.~Class, scale = "free_x", space = "free")  +
  theme_bw() +
  scale_color_manual(values = region.pal) +
  scale_fill_manual(values = region.pal) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.4)) +
  xlab("Compound") +
  ylab("Particulate Concentration (nM)")

p.boxplot

ggsave(p.boxplot, dpi = 300, width = 9, height = 5, scale = 1.3, units = "in", file = "Figures/Survey_Particulate_Boxplot.png")



#dissolved_boxplot 
d.boxplot <- ggplot(fig.dat, aes(x = reorder(Compound, -rank), y = Diss.Conc.nM, fill = Region)) +
  geom_boxplot(aes(color = Region), alpha = 0.5,  position = position_dodge(preserve = "single")) +
  geom_point(aes(fill = Region, group = Region),  alpha = 0.7, shape = 21, stroke = 0.22, size = 2, 
             position=position_jitterdodge(jitter.width = 0.15)) + 
  scale_y_log10() +
  facet_grid(.~Class, scale = "free_x", space = "free")  +
  theme_bw() +
  scale_color_manual(values = region.pal) +
  scale_fill_manual(values = region.pal) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.4)) +
  xlab("Compound") +
  ylab("Dissolved Concentration (nM)")

d.boxplot

ggsave(d.boxplot, dpi = 300, width = 9, height = 5, scale = 1.3, units = "in", file = "Figures/Survey_Dissolved_Boxplot.png")




##########Dissolve and Particulate Correlation Plot:
fig.no.all <- fig.dat %>%
  select(-Local_Time) %>%
  select(-Local_Date) %>%
  unique() %>%
  filter(!Region == "All") %>%
  unique()

#make pal with no all
region.pal.no.all <- c("#6794a7", "#76c0c1", "#014d64", "#01a2d9", "#00887d")

part_diss_cor_plot <- ggplot(fig.no.all, aes(x = Part.Conc.nM, y = Diss.Conc.nM)) +
  geom_point(aes(fill = Region),  alpha = 0.7, shape = 21, stroke = 0.22, size = 2) +
  scale_fill_manual(values = region.pal.no.all) +
  geom_smooth(method = "lm", alpha = 0.5, color = "black") +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(.~Compound, scales = "free") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("Particulate Concentration (nM)") +
  ylab("Dissolved Concentration (nM)")
part_diss_cor_plot

ggsave(part_diss_cor_plot, dpi = 300, width = 9, height = 5, scale = 1.3, units = "in", file = "Figures/Survey_CorPlot.png")








##############Summary Statistics Table
fig.dat.sum <- fig.dat %>%
  select(-Local_Time) %>%
  select(-Local_Date) %>%
  unique() %>%
  group_by(Region, Class, Compound) %>%
  reframe(n = n(),
          Mean_Particulate_Conc_nM = mean(Part.Conc.nM, na.rm = TRUE),
          SD_Particulate_Conc_nM = sd(Part.Conc.nM, na.rm = TRUE),
          Mean_Dissolved_Conc_nM = mean(Diss.Conc.nM, na.rm = TRUE),
          SD_Dissolved_Conc_nM = sd(Diss.Conc.nM, na.rm = TRUE),
          Mean_Percent_Dissolved = mean(Perc_Diss, na.rm = TRUE),
          SD_Percent_Dissolved = sd(Perc_Diss, na.rm = TRUE))



########### Convergent Betaine Analysis 


#####______________________PARTICULATE_______________________
betaine.dat.p <- fig.dat %>%
  select(-Local_Time) %>%
  select(-Local_Date) %>%
  ungroup() %>%
  unique() %>%
  filter(!Region == "All") %>%
  unique() %>%
  filter(Class == "Betaine",
         !is.na(Part.Conc.nM)) %>%
  group_by(Part.Rep) %>%
  mutate(rel.conc = Part.Conc.nM/Part.Conc.nM[Compound == "Glycine betaine"],
         betaine.rank = rank(1/Part.Conc.nM)) %>%
  ungroup()

#relative concentration:
ggplot(betaine.dat.p, aes(x = rel.conc, y = reorder(Compound, rank), fill = Region)) +
  geom_jitter(height = 0.1, width = 0, alpha = 0.5, shape = 21) +
  scale_x_log10() +
  facet_wrap(.~Region)

#rank order points:
ggplot(betaine.dat.p, aes(x = betaine.rank, y = reorder(Compound, rank), fill = Region)) +
  geom_jitter(height = 0.1, width = 0.1, alpha = 0.5, shape = 21) +
  facet_wrap(.~Region)

#rank order heatmap:
betaine.dat.p.ranksum <- betaine.dat.p %>%
  group_by(Region, Compound) %>%
  mutate(num_obs = n()) %>%
  group_by(Region, Compound, betaine.rank, rank, num_obs) %>%
  reframe(count = n(),
          rel_count = count/(num_obs))

ggplot(betaine.dat.p.ranksum, aes(x = betaine.rank, y = reorder(Compound, rank), fill = rel_count)) +
  geom_tile() +
  scale_fill_viridis(option = "G", end = 0.8, begin = 0.1) +
  facet_wrap(.~Region)
  
#Particulate Betaine Convergence summary:
betaine.converg.part.sum <- betaine.dat.p %>%
  group_by(Compound, Region) %>%
  reframe(mean.rel.conc = mean(rel.conc),
          sd.rel.conc = sd(rel.conc),
          mean.rank = mean(betaine.rank),
          sd.rank = sd(betaine.rank))

betaine.overall.ratio <- betaine.dat.p %>%
  group_by(Compound) %>%
  reframe(mean.rel.conc = mean(rel.conc),
          sd.rel.conc = sd(rel.conc),
          mean.rank = mean(betaine.rank),
          sd.rank = sd(betaine.rank),
          log.rel.conc = log10(mean.rel.conc),
          log.sd.rel.conc = log10(sd.rel.conc))

#overall correlation among betaines:
betaine.dat.a <- betaine.dat.p %>%
  select(Part.Rep, Region, Compound, Part.Conc.nM) %>%
  rename("Compound_A" = Compound,
         "Conc_A" = Part.Conc.nM)

betaine.dat.b <- betaine.dat.p %>%
  select(Part.Rep, Region, Compound, Part.Conc.nM) %>%
  rename("Compound_B" = Compound,
         "Conc_B" = Part.Conc.nM)

betaine.dat.cor <- left_join(betaine.dat.a, betaine.dat.b) %>%
  ungroup() %>%
  group_by(Compound_A, Compound_B) %>%
  reframe(cor = cor(Conc_A, Conc_B))

#Cor Viz
ggplot(betaine.dat.cor, aes(x = cor, y = Compound_A, fill = Compound_B)) +
  geom_point(size = 3, shape = 21) +
  xlim(0,1)
  

#####______________________DISSOLVED_______________________
betaine.dat.d <- fig.dat %>%
  select(-Local_Time) %>%
  select(-Local_Date) %>%
  ungroup() %>%
  unique() %>%
  filter(!Region == "All") %>%
  filter(!Cruise == "KM1906") %>%
  unique() %>%
  filter(Class == "Betaine",
         !is.na(Diss.Conc.nM)) %>%
  group_by(Diss.Rep) %>%
  mutate(rel.conc = Diss.Conc.nM/Diss.Conc.nM[Compound == "Glycine betaine"],
         betaine.rank = rank(1/Diss.Conc.nM)) %>%
  ungroup()

#relative concentration:
ggplot(betaine.dat.d, aes(x = rel.conc, y = reorder(Compound, rank), fill = Region)) +
  geom_jitter(height = 0.1, width = 0, alpha = 0.5, shape = 21) +
  scale_x_log10() +
  facet_wrap(.~Region)

#rank order:
ggplot(betaine.dat.d, aes(x = betaine.rank, y = reorder(Compound, rank), fill = Region)) +
  geom_jitter(height = 0.1, width = 0.1, alpha = 0.5, shape = 21) +
  facet_wrap(.~Region)

#Particulate Betaine Convergence summary:
betaine.converg.diss.sum <- betaine.dat.d %>%
  group_by(Compound, Region) %>%
  reframe(mean.rel.conc = mean(rel.conc),
          sd.rel.conc = sd(rel.conc),
          mean.rank = mean(betaine.rank),
          sd.rank = sd(betaine.rank))

betaine.overall.ratio <- betaine.dat.d %>%
  group_by(Compound) %>%
  reframe(mean.rel.conc = mean(rel.conc),
          sd.rel.conc = sd(rel.conc),
          mean.rank = mean(betaine.rank),
          sd.rank = sd(betaine.rank),
          log.rel.conc = log10(mean.rel.conc),
          log.sd.rel.conc = log10(sd.rel.conc))

#overall correlation among betaines:
betaine.dat.a <- betaine.dat.d %>%
  select(Diss.Rep, Region, Compound, Diss.Conc.nM) %>%
  rename("Compound_A" = Compound) %>%
  mutate("Conc_A" = log10(Diss.Conc.nM)) %>%
  select(-Diss.Conc.nM)

betaine.dat.b <- betaine.dat.d %>%
  select(Diss.Rep, Region, Compound, Diss.Conc.nM) %>%
  rename("Compound_B" = Compound) %>%
  mutate("Conc_B" = log10(Diss.Conc.nM)) %>%
  select(-Diss.Conc.nM)

betaine.dat.cor <- left_join(betaine.dat.a, betaine.dat.b) %>%
  ungroup() %>%
  group_by(Compound_A, Compound_B) %>%
  reframe(cor = cor(Conc_A, Conc_B))

#Cor Viz
ggplot(betaine.dat.cor, aes(x = cor, y = Compound_A, fill = Compound_B)) +
  geom_point(size = 3, shape = 21) +
  xlim(0,1)








g.dat <-dp.dat %>%
  filter(!is.na(Lat)) %>%
  filter(!is.na(Compound)) %>%
  unique()


##make some preliminary figures:

#part_conc:
ggplot(g.dat, aes(x = Lat, y = Part.Conc.nM, color = Long)) +
  geom_point(alpha = 0.9, size = 3) +
  facet_wrap(.~Compound, scales = "free_y") +
  scale_color_viridis(option = "B", end = 0.85) #+
 # scale_y_log10()
ggsave(filename = "Figures/all_comp_particulate_gradeints_transect.png",
       height = 8, width = 13)

#diss_conc
ggplot(g.dat, aes(x = Lat, y = Diss.Conc.nM, color = Long)) +
  geom_point(alpha = 0.9, size = 3) +
  facet_wrap(.~Compound, scales = "free_y") +
  scale_color_viridis(option = "B", end = 0.85)# +
 # scale_y_log10() 
ggsave(filename = "Figures/all_comp_dissolved_gradeints_transect.png",
       height = 8, width = 13)

#combined:
d.p.transect <- ggplot(g.dat) +
  geom_jitter(aes(x = Lat, y = Diss.Conc.nM), width = 0.1,
              fill = "skyblue",
              alpha = 0.9, size = 2.5, shape = 21, stroke = 0.1) +
  geom_jitter(aes(x = Lat, y = Part.Conc.nM), width = 0.1,
              alpha = 0.9, size = 2.5, 
              fill = "#bf616a", 
              shape = 24, stroke = 0.1) +
  facet_wrap(Compound~., scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  xlab("Latitude") +
  ylab("Concentration (nM)")

ggsave(d.p.transect, filename = "Figures/all_comp_dissolved_particualte_gradeints_transect.png",
       height = 6, width = 12, scale = 1.1)



#perc_diss
ggplot(g.dat, aes(x = Lat, y = Perc_Diss, fill = Long)) +
  geom_point(alpha = 0.7, size = 3, shape = 21) +
  facet_wrap(.~Compound, scales = "free_y") +
  scale_fill_viridis(option = "B", end = 0.85) +
  ylim(0,1)
ggsave(filename = "Figures/all_comp_percent_dissolved_gradeints_transect.png",
       height = 6, width = 12, scale = 1.1)



ggplot(g.dat, aes(x = Part.Conc.nM, y = Diss.Conc.nM, color = Lat)) +

  scale_color_viridis(option = "B", end = 0.85) +
  geom_point(alpha = 0.9) +
  geom_smooth(method = "lm", alpha = 0.5) +
  facet_wrap(.~Compound, scales = "free") +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() 
  
























































