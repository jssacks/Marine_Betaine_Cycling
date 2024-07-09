


library(tidyverse)


diss.dat <- read_csv("Intermediates/Dissolved_Quantified_Data.csv")
g3.samp.loc <- read_csv("Meta_Data/Ingalls_Lab_Data/G3_Samp_Locations.csv")
g4.samp.loc <- read_csv("Meta_Data/Ingalls_Lab_Data/G4_Station_Locations.csv")
D1.samp.dat <- read_csv("Meta_Data/Ingalls_Lab_Data/RC078_metadata.csv")
p.q.dat <- read_csv("Intermediates/Particulate_Quant_Output.csv")

#LODs:
lod.dat <- read_csv("Intermediates/Dissolved_Blk_LOD_Concentrations.csv")





hom.dat <- diss.dat %>%
  filter(Name == "Homarine") %>%
  filter(Cruise == "TN397") %>%
  left_join(., lod.dat) 
ggplot(hom.dat, aes(x = SampID, y = EE.adjust.conc))



diss.dat.2 <- diss.dat %>%
  mutate(flag = case_when(Name == "Glycine betaine" & Cruise == "KM1906" ~ "flag",
                          TRUE ~ NA)) %>%
  filter(is.na(flag)) %>%
  select(-flag) %>%
  filter(!is.na(Cruise)) %>%
  filter(!str_detect(SampID, "Blk")) %>%
  filter(!str_detect(SampID, "CXBLK")) %>%
  filter(!str_detect(SampID, "Poo")) %>%
  filter(!str_detect(SampID, "MQ")) %>%
  filter(!str_detect(SampID, "Mix")) %>%
  left_join(., lod.dat) %>%
  group_by(Name) %>%
  mutate(rank.order = rank(EE.adjust.conc))

ggplot(diss.dat.2, aes(x = Name, y = EE.adjust.conc, color = Name)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(.~Cruise, scales = "free") +
  scale_y_log10() +
  geom_hline(aes(yintercept = EE.adjust.lod, color = Name))


ggplot(diss.dat.2, aes(x = rank.order, y = EE.adjust.conc, color = Cruise)) +
  geom_point(alpha = 0.6, size = 2,  color = "red") +
  geom_point(data = p.q.rank, aes(x = rank.order, y = nM.in.smp), color = "blue") +
  facet_wrap(.~Name, scales = "free_y") +
  scale_y_log10() +
  geom_hline(yintercept = 0.1)
#  geom_hline(aes(yintercept = EE.adjust.lod, color = Cruise))



#########
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
metab.loc.dat <- diss.dat.2 %>%
  mutate(Samp_ID = str_remove(SampID, "220623_Smp_")) %>%
  mutate(Samp_ID = str_remove(Samp_ID, "220602_Smp_TN397_")) %>%
  mutate(Samp_ID = paste(Cruise, Samp_ID, sep = "_")) %>%
  filter(Cruise %in% c("KM1906", "TN397")) %>%
  left_join(., all.loc)



ggplot(metab.loc.dat, aes(x = Lat, y = EE.adjust.conc)) +
  geom_point(alpha = 0.8, size = 3, stroke = 0.2, shape=21, aes(fill = Long)) +
  facet_wrap(.~Name, scales = "free_y") +
  scale_y_log10() +
  geom_hline(aes(yintercept = EE.adjust.lod, color = Cruise))



#######combine location and metabolite data:
p.metab.loc.dat <- p.q.dat %>%
  mutate(Samp_ID = str_remove(SampID, "220628_Smp_")) %>%
  mutate(Samp_ID = str_remove(Samp_ID, "220902_Smp_TN397_")) %>%
  mutate(Samp_ID = paste(Cruise, Samp_ID, sep = "_")) %>%
  filter(!Cruise == "RC078") %>%
  left_join(., all.loc)



both.metab.loc.dat <- metab.loc.dat %>%
  select(Name, Cruise, Samp_ID, EE.adjust.conc) %>%
  full_join(., p.metab.loc.dat) %>%
  mutate(perc.diss = EE.adjust.conc/(nM.in.smp+EE.adjust.conc))# %>%
#  filter(!is.na(nM.in.smp)) #%>%
  #filter(!is.na(EE.adjust.conc))

ggplot(both.metab.loc.dat, aes(x = Lat, y = perc.diss, color = Long)) +
  geom_point(alpha = 0.8, size = 1.5) +
  facet_wrap(.~Name)# +
  #scale_y_log10() 


###Export data:
dat.exp <- both.metab.loc.dat %>%
  ungroup() %>%
  # filter(Long < -135) %>%
  mutate(Region = case_when(Lat > 20 & Long > -130 ~ "California_Current", 
                            Lat < 35 & Lat > 10 & Long < -130 ~ "Gyre",
                            Lat < 10 & Long < -130 ~ "Equatorial",
                            Lat > 35 & Long < -130 ~ "NPTZ" )) %>%
  filter(!is.na(Region))  %>%
  filter(EE.adjust.conc > 0.001 | is.na(EE.adjust.conc)) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("California_Current", "Gyre", "Equatorial", "NPTZ")))
  
write_csv(dat.exp, file = "Intermediates/exported_d_p_location_gradients_data.csv")


region.pal <- c("#76c0c1", "#014d64", "#01a2d9", "#00887d")

#Dissolved Boxplot
dis.bp <- ggplot(dat.exp, aes(x = Name, y = EE.adjust.conc)) +
  geom_boxplot(aes(color = Region), alpha = 0.5,  position = position_dodge(preserve = "single")) +
  geom_point(aes(fill = Region, group = Region),  alpha = 0.7, shape = 21, stroke = 0.22, size = 2, position=position_jitterdodge(jitter.width = 0.15)) +
#  facet_wrap(.~Region) +
  scale_y_log10() +
  scale_color_manual(values = region.pal) +
  scale_fill_manual(values = region.pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.65)) +
  xlab("Compound") +
  ylab("Dissolved Concentration (nM)")
dis.bp
ggsave(dis.bp, dpi = 300, width = 6, height = 4, scale = 1.1, units = "in", file = "Figures/Dissolved_Boxplot.png")




#Particulate Boxplot
part.bp <- ggplot(dat.exp, aes(x = Name, y = nM.in.smp)) +
  geom_boxplot(aes(color = Region), alpha = 0.5,  position = position_dodge(preserve = "single")) +
  geom_point(aes(fill = Region, group = Region),  alpha = 0.7, shape = 21, stroke = 0.22, size = 2, position=position_jitterdodge(jitter.width = 0.15)) +
  #  facet_wrap(.~Region) +
  scale_y_log10() +
  scale_color_manual(values = region.pal) +
  scale_fill_manual(values = region.pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.65)) +
  xlab("Compound") +
  ylab("Particulate Concentration (nM)")
part.bp
ggsave(part.bp, dpi = 300, width = 6, height = 4, scale = 1.1, units = "in", file = "Figures/Particulate_Boxplot.png")


#Make Homarine Plot
hom.dat.exp <- dat.exp %>%
  filter(Name == "Homarine") %>%
  filter(!Region == "California_Current")

hom.plot <- ggplot(hom.dat.exp) +
  geom_jitter(aes(x = Lat, y = EE.adjust.conc), width = 0.21,
              fill = "skyblue",
             alpha = 0.9, size = 3, shape = 21, stroke = 0.25) +
  geom_jitter(aes(x = Lat, y = nM.in.smp, fill = Region), width = 0.21,
             alpha = 0.9, size = 3, 
             fill = "#bf616a", 
             shape = 24, stroke = 0.25) +
  scale_y_log10() +
  theme_bw() +
  xlab("Latitude") +
  ylab("Concentration (nM)")
ggsave(hom.plot, dpi = 300, width = 6, height = 4, scale = 1.1, units = "in", file = "Figures/hom_transect.png")


#Make Gbt Plot
gbt.dat.exp <- dat.exp %>%
  filter(Name == "Glycine betaine") %>%
  filter(!Region == "California_Current")

gbt.plot <- ggplot(gbt.dat.exp) +
  geom_jitter(aes(x = Lat, y = EE.adjust.conc), width = 0.21,
              fill = "skyblue",
              alpha = 0.9, size = 3, shape = 21, stroke = 0.25) +
  geom_jitter(aes(x = Lat, y = nM.in.smp, fill = Region), width = 0.21,
              alpha = 0.9, size = 3, 
              fill = "#bf616a", 
              shape = 24, stroke = 0.25) +
  scale_y_log10() +
  theme_bw() +
  xlab("Latitude") +
  ylab("Concentration (nM)")
ggsave(gbt.plot, dpi = 300, width = 6, height = 4, scale = 1.1, units = "in", file = "Figures/gbt_transcet.png")


#Make Gbt Plot
bab.dat.exp <- dat.exp %>%
  filter(Name == "beta-Alaninebetaine") %>%
  filter(!Region == "California_Current")

bab.plot <- ggplot(bab.dat.exp) +
  geom_jitter(aes(x = Lat, y = EE.adjust.conc), width = 0.21,
              fill = "skyblue",
              alpha = 0.9, size = 3, shape = 21, stroke = 0.25) +
  geom_jitter(aes(x = Lat, y = nM.in.smp, fill = Region), width = 0.21,
              alpha = 0.9, size = 3, 
              fill = "#bf616a", 
              shape = 24, stroke = 0.25) +
  scale_y_log10() +
  theme_bw() +
  xlab("Latitude") +
  ylab("Concentration (nM)")
bab.plot
ggsave(bab.plot, dpi = 300, width = 6, height = 4, scale = 1.1, units = "in", file = "Figures/bab_transcet.png")



####Plot of all compounds on in Dissolved vs. Particulate Space
dvp.plot <- ggplot(dat.exp, aes(x = nM.in.smp, y = EE.adjust.conc)) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(values = region.pal) +
  geom_point(aes(fill = Region), shape = 21, 
             size = 2, alpha = 0.8, stroke = 0.25) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope = 1, color = "black", linetype = "dashed") +
  theme_bw() +
  facet_wrap(.~Name, scales = "free") +
  xlab("Particulate Concentration (nM)") +
  ylab("Dissolved Concentration (nM)") +
  theme(legend.position = "bottom")
dvp.plot
ggsave(dvp.plot, dpi = 300, width = 6, height = 5, scale = 1.1, units = "in", file = "Figures/betaine_dis_v_part.png")

ggsave(dvp.plot, dpi = 300, width = 6, height = 5, scale = 1.1, units = "in", file = "Figures/betaine_dis_v_part.svg")


###Facet by region and compound:
ggplot(dat.exp, aes(x = nM.in.smp, y = EE.adjust.conc)) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(values = region.pal) +
  geom_point(aes(fill = Region), shape = 21, 
             size = 2, alpha = 0.8, stroke = 0.25) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope = 1, color = "black", linetype = "dashed") +
  theme_bw() +
  facet_grid(Region~Name, scales = "free") +
  xlab("Particulate Concentration (nM)") +
  ylab("Dissolved Concentration (nM)") +
  theme(legend.position = "bottom")








ggplot(dat.exp, aes(x = Name, y = nM.in.smp, color = Region)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(preserve = "single")) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1)) +
  #  facet_wrap(.~Region) +
  scale_y_log10()


write_csv(dat.exp, file = "Intermediates/Gradients_D_P_dat.csv")





av.perc.diss <- both.metab.loc.dat %>%
  ungroup() %>%
 # filter(Long < -135) %>%
  mutate(Region = case_when(Lat > 20 & Long > -130 ~ "California_Current", 
                            Lat < 35 & Lat > 10 & Long < -130 ~ "Gyre",
                            Lat < 10 & Long < -130 ~ "Equator",
                            Lat > 35 & Long < -130 ~ "NPTZ" ))  %>%
  group_by(Name, Region) %>%
  summarize(mean.perc.diss = mean(perc.diss),
            sd.perc.diss = sd(perc.diss),
            max.perc.diss = max(perc.diss),
            min.perc.diss = min(perc.diss),
            median.perc.diss = median(perc.diss),
            mean.nM.part = mean(nM.in.smp),
            sd.nM.part = sd(nM.in.smp),
            max.nM.part = max(nM.in.smp),
            min.nM.part = min(nM.in.smp),
            median.nM.part = median(nM.in.smp),
            mean.nM.diss = mean(EE.adjust.conc),
            sd.nM.diss = sd(EE.adjust.conc),
            max.nM.diss = max(EE.adjust.conc),
            min.nM.diss = min(EE.adjust.conc),
            median.nM.diss = median(EE.adjust.conc))

ggplot(av.perc.diss, aes(y = Name, x = mean.perc.diss, color = Region)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmax = mean.perc.diss+sd.perc.diss, xmin = mean.perc.diss-sd.perc.diss), 
                 height = 0.4, linewidth = 0.5, position = position_dodge(width = 0.5)) +
  xlim(0,1)
  
###write summarized data to csv:
write_csv(av.perc.diss, file = "Intermediates/Region_Mean_Dat_Prelim.csv")









ggplot(both.metab.loc.dat, aes(x = nM.in.smp, y = EE.adjust.conc, color = Lat)) +
  geom_point() +
  facet_wrap(.~Name, scales = "free")

ggplot(both.metab.loc.dat, aes(x = nM.in.smp, y = EE.adjust.conc, color = Lat)) +
  geom_smooth(method = "lm", color = "black", alpha = 0.8) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
 # coord_fixed(ratio = 1) +
  ylab("Dissolved Concentration (nM)") +
  xlab("Particulate Concentration (nM)") +
  facet_wrap(.~Name, scales = "free")  +
  geom_abline(intercept = 0, slope = 1, color = "red")




############DiNiMITE 1:
diss.dat.3 <- diss.dat %>%
  mutate(flag = case_when(Name == "Glycine betaine" & Cruise == "KM1906" ~ "flag",
                          TRUE ~ NA)) %>%
  filter(is.na(flag)) %>%
  select(-flag) %>%
  filter(!is.na(Cruise)) %>%
  filter(!str_detect(SampID, "Blk")) %>%
  filter(!str_detect(SampID, "CXBLK")) %>%
  filter(!str_detect(SampID, "Poo")) %>%
  filter(!str_detect(SampID, "MQ")) %>%
  left_join(., lod.dat) %>%
  filter(Cruise == "RC078") %>%
  mutate(sample_id = str_remove(SampID, "240415_Smp_")) %>%
  left_join(., D1.samp.dat) %>%
  select(Name, sample_id, parent_id, EE.adjust.conc, station, depth_m) %>%
  unique() 



###partiuclate
d1.all.dat <- p.q.dat %>%
  filter(Cruise == "RC078") %>%
  mutate(sample_id = str_remove(SampID, "221006_Smp_")) %>%
  right_join(., diss.dat.3) %>%
  mutate(perc.diss = EE.adjust.conc/(EE.adjust.conc+nM.in.smp)) %>%
  filter(!station == 4)

ggplot(d1.all.dat, aes(x = Name, y = perc.diss)) +
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, aes(color = as.factor(station)))

ggplot(d1.all.dat, aes(x = EE.adjust.conc, y = -depth_m, color = Name)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~station) +
  scale_x_log10()

ggplot(d1.all.dat, aes(x = EE.adjust.conc, y = nM.in.smp, color = as.factor(station)))+ 
  geom_point() +
  facet_wrap(.~Name, scales = "free") #+
 # scale_y_log10() +
 # scale_x_log10()

ggplot(d1.all.dat, aes(y = EE.adjust.conc, x = Name, color = Name)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  ylab("Dissolved Concentration (nM") +
  facet_wrap(.~station, scales = "free")# +
 # scale_y_log10() 




s2 <- d1.all.dat %>%
  filter(station == 2)

ggplot(s2, aes(x = EE.adjust.conc, y = nM.in.smp, color = as.factor(station)))+ 
  geom_point() +
  facet_wrap(.~Name, scales = "free") +
  scale_y_log10() +
  scale_x_log10()



















