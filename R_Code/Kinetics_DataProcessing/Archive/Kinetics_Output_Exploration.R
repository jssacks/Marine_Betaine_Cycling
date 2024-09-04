
library(tidyverse)



flux.dat <- read_csv("Intermediates/Flux_Analysis_Output.csv") %>%
  separate(cruise_exp, into = c("Cruise", "exp"), remove = F)
kin.dat <- read_csv("Intermediates/Kin_Analysis_Output.csv")




kf.dat <- left_join(flux.dat, kin.dat) 

##Ks viz
ggplot(kf.dat, aes(x = cruise_exp, y = mean_ks)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_ks-sd_ks, ymax = mean_ks+sd_ks), width = 0.2) +
  facet_grid(.~Compound, scales = "free_x", space = "free") +
#  scale_y_log10() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) 

##Vmax viz

library(ggbreak)

ggplot(kf.dat, aes(x = cruise_exp, y = mean_vmax)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_vmax-sd_vamx, ymax = mean_vmax+sd_vamx), width = 0.2) +
  facet_grid(.~Compound, scales = "free_x", space = "free") +
  scale_y_break(c(0.75, 3), scales = "0.3", ticklabels = c(3.0, 3.5, 4.0)) +
  #  scale_y_log10() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) 







ggplot(kf.dat, aes(x = Mean.Diss.Conc.nM, y = mean_ks, color = Compound)) +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
#  coord_fixed(xlim = ) +
  scale_x_log10() +
  scale_y_log10()

ggplot(kf.dat, aes(x = Mean.Diss.Conc.nM, y = mean_ks, color = Compound)) +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(.~Compound, scales = "free")
  #  coord_fixed(xlim = ) +
 # scale_x_log10() +
 # scale_y_log10()
open_ocean <- kf.dat %>%
  filter(Cruise %in% c("KM1906", "TN397", "TN412"))

ggplot(open_ocean, aes(x = Mean.Diss.Conc.nM, y = mean_ks, color = Compound)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(.~Compound, scales = "free")



























