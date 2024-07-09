
library(tidyverse)
library(viridis)

flux.in.dat <- data.frame(
  V = seq(0.001, 0.1, by=0.001)
)

Vmax.dat <- data.frame(
  Vmax = seq(0.01, 10, by=0.01)
)

Ks.dat <- data.frame(
  Ks = seq(0.1, 10, by=0.1)
)

Flux.ratio.dat <- data.frame(
  Part_Fr = c(1, 2, 3, 4, 5, 10)
)


mod.dat <- cross_join(flux.in.dat, Vmax.dat) %>%
  cross_join(., Ks.dat) %>%
  cross_join(., Flux.ratio.dat) %>%
  mutate(S_conc = (V*Ks)/(Vmax-V)) %>%
  mutate(S_V_ratio = S_conc/V) %>%
  mutate(P_conc = V/(Part_Fr/100)) %>%
  mutate(Perc_Diss = S_conc/(S_conc + P_conc))


ks.sub.mod <- mod.dat %>%
  filter(Ks %in% c(0.1, 1, 5, 10))

#ggplot(ks.sub.mod, aes(x = V, y = Vmax, fill = log10(S_conc))) +
#  geom_tile() +
#  facet_wrap(.~Ks) +
#  scale_fill_viridis()



ks.sub.mod <- mod.dat %>%
  filter(Ks %in% c(0.1, 1, 5, 10)) %>%
 # filter(V %in% c(0.001, 0.0026, 0.005, 0.075, 0.01, 0.025, 0.05, 0.075, 0.1)) %>%
  filter(Vmax %in% c(0.01, 0.025, 0.05, 0.075, 0.1, 0.5, 1, 5, 10)) %>%
  filter(S_conc > 0.0001) %>%
  filter(S_conc < 10)


ggplot(ks.sub.mod, aes(x = V, y = S_conc, color = as.factor(Vmax))) +
  geom_point() +
  facet_wrap(.~Ks) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10() +
  geom_hline(yintercept = 0.01)


ggplot(ks.sub.mod, aes(x = V, y = S_V_ratio, color = as.factor(Vmax))) +
  geom_point() +
  facet_wrap(.~Ks) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10() +
  geom_hline(yintercept = 10) +
  geom_hline(yintercept = 3, color = "red") +
  geom_hline(yintercept = 1, color = "blue")




Perc.diss.mod <- mod.dat %>%
  filter(Ks %in% c(0.1, 1, 5, 10)) %>%
  # filter(V %in% c(0.001, 0.0026, 0.005, 0.075, 0.01, 0.025, 0.05, 0.075, 0.1)) %>%
  filter(Vmax %in% c(0.01, 0.03, 0.05, 0.07, 0.11, 0.30, 0.5, 0.7, 1, 5, 10)) %>%
  filter(S_conc > 0.001) %>%
  filter(S_conc < 10) 

  
ggplot(Perc.diss.mod, aes(x = Part_Fr, y = Perc_Diss, color = log10(S_conc))) +
  geom_point() +
  facet_grid(Ks~Vmax) +
  scale_color_viridis() +
  geom_hline(yintercept = 0.9) +
  geom_hline(yintercept = 0.75, color = "red")


ggplot(Perc.diss.mod, aes(x = P_conc, y = Perc_Diss, color = as.factor(Part_Fr))) +
  geom_point() +
  facet_grid(Ks~Vmax) +
  scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0.9) +
  geom_hline(yintercept = 0.7) +
  geom_vline(xintercept = 0.01) +
  geom_vline(xintercept = 1) +
#  geom_hline(yintercept = 0.75, color = "red") +
#  geom_vline(xintercept = 0.01, color = "blue") +
#  geom_vline(xintercept = 0.1, color = "purple") +
#  geom_vline(xintercept = 1, color = "orange") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_log10() 





#######Pull out just the model data that fits into the parameters of different CompoundXregions:

comp.av.dat <- read_csv("Intermediates/Region_Mean_Dat_Prelim.csv")

comp.ranges <- comp.av.dat %>%
  mutate(perc.diss.u.b = mean.perc.diss+sd.perc.diss,
         perc.diss.l.b = mean.perc.diss-sd.perc.diss,
         part.conc.u.b = mean.nM.part+sd.nM.part/2,
         part.conc.l.b = mean.nM.part-sd.nM.part/2,
         diss.conc.u.b = mean.nM.diss+sd.nM.diss/2,
         diss.conc.l.b = mean.nM.diss-sd.nM.diss/2) %>%
  select(Name, Region, perc.diss.u.b, perc.diss.l.b, part.conc.u.b,
         part.conc.l.b, diss.conc.u.b, diss.conc.l.b)

hom.ranges <- comp.ranges %>%
  filter(Name == "Trigonelline")

hom.mod.vals <- mod.dat %>%
  filter(Perc_Diss > min(hom.ranges$perc.diss.l.b)) %>%
  filter(Perc_Diss < max(hom.ranges$perc.diss.u.b)) %>%
  filter(S_conc > min(hom.ranges$diss.conc.l.b)) %>%
  filter(S_conc < max(hom.ranges$diss.conc.u.b)) %>%
  cross_join(., hom.ranges) %>%
  filter(Perc_Diss > perc.diss.l.b & Perc_Diss < perc.diss.u.b,
         S_conc > diss.conc.l.b & S_conc < diss.conc.u.b)
#         P_conc > part.conc.l.b & P_conc < part.conc.u.b)

gyre.vals <- hom.mod.vals  %>%
  filter(Region == "Gyre")

ggplot(gyre.vals, aes(x = Ks, y = Vmax, color = P_conc)) +
  geom_point(alpha = 0.4, size = 2) +
  facet_wrap(.~as.factor(Part_Fr)) +
  #scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 1)

ggplot(hom.mod.vals, aes(x = P_conc, y = S_conc, color = Perc_Diss)) +
  geom_point() +
  facet_wrap(.~Region)

gyre.vals.2 <- hom.mod.vals  %>%
  filter(Region == "Gyre") #%>%
  #filter(Vmax < 0.1)


ggplot(gyre.vals.2, aes(x = Ks, y = Vmax, color = S_conc)) +
  geom_point(alpha = 0.3, size = 3) +
  facet_wrap(.~as.factor(Part_Fr)) +
  geom_vline(xintercept = 1) +
  #scale_x_log10() +
  scale_y_log10()



########################

#Trying to fit each point individually:
dat.gradients <- read_csv("Intermediates/Gradients_D_P_dat.csv")

dat.g.hom <- dat.gradients %>%
  filter(Name == "Homarine") %>%
  select(Name, Region, Lat, Long, EE.adjust.conc, nM.in.smp, Samp_ID) %>%
  rename(diss.conc.env = EE.adjust.conc,
         part.conc.env = nM.in.smp) %>%
  mutate(max.p = part.conc.env+(0.05*part.conc.env),
            min.p = part.conc.env-(0.05*part.conc.env),
            max.d = diss.conc.env+(0.05*diss.conc.env),
            min.d = diss.conc.env-(0.05*diss.conc.env)) %>%
  head(10)
  
  


dat.hom.ranges <- dat.g.hom %>%
  summarize(max.p = max(nM.in.smp),
            min.p = min(nM.in.smp),
            max.d = max(EE.adjust.conc),
            min.d = min(EE.adjust.conc),
            max.p.d = max(perc.diss),
            min.p.d = min(perc.diss))

hom.mod.vals <- mod.dat %>%
  filter(P_conc < dat.hom.ranges$max.p & P_conc > dat.hom.ranges$min.p) %>%
  filter(S_conc < dat.hom.ranges$max.d & S_conc > dat.hom.ranges$min.d)  %>%
  filter(Perc_Diss < dat.hom.ranges$max.p.d & Perc_Diss > dat.hom.ranges$min.p.d) %>%
  filter(Vmax < 0.2) %>%
  cross_join(., dat.g.hom) 


hom.filt.dat <- hom.mod.vals %>%
  filter(P_conc > min.p & P_conc < max.p) %>%
  filter(S_conc > min.d & S_conc < max.d) 

ggplot(hom.filt.dat, aes(x = Ks, y = Vmax, color = Samp_ID)) +
  geom_point(alpha = 0.4) +
  facet_wrap(.~Part_Fr) 



hom.min.ks.dat <- hom.filt.dat %>%
  group_by(Samp_ID, Part_Fr) %>%
  mutate(min.Ks = min(Ks)) %>%
  filter(Ks == min.Ks)

ggplot(hom.min.ks.dat, aes(x = Ks, y = Vmax, color = Perc_Diss)) +
  geom_point(size = 3) +
  facet_wrap(.~Part_Fr) +
  scale_x_log10()

ggplot(hom.min.ks.dat, aes(x = Ks, y = Vmax, color = Lat)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  facet_wrap(.~Part_Fr) +
  scale_x_log10()

ggplot(hom.min.ks.dat, aes(x = Ks, fill = Lat)) +
  geom_histogram() +
  facet_wrap(.~Part_Fr) +
  scale_x_log10()

ggplot(hom.min.ks.dat, aes(x = Vmax, fill = Lat)) +
  geom_histogram() +
  facet_wrap(.~Part_Fr)# +
 # scale_x_log10()

ggplot(hom.min.ks.dat, aes(x = Lat, y = Ks, color = Perc_Diss)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  facet_wrap(.~Part_Fr) +
  scale_y_log10()

ggplot(hom.min.ks.dat, aes(x = Lat, y = Vmax, color = Perc_Diss)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  facet_wrap(.~Part_Fr) #+
  #scale_y_log10()

ggplot(hom.min.ks.dat, aes(x = S_conc, y = Ks, color = Perc_Diss)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  facet_wrap(.~Part_Fr) +
  scale_y_log10() +
  scale_x_log10()

ggplot(hom.min.ks.dat, aes(x = V, y = Vmax, color = log10(S_conc))) +
  geom_point(size = 3) +
  scale_color_viridis() +
  facet_wrap(.~Part_Fr) +
 # scale_y_log10() +
  scale_x_log10()

sum.min.ks <- hom.min.ks.dat %>%
  group_by(Region, Part_Fr) %>%
  summarize(count = n(),
            mean_ks = mean(Ks),
            sd_ks = sd(Ks),
            mean_vmax = mean(Vmax),
            sd_vmax = sd(Vmax))



library(viridis)














































































