

library(tidyverse)



#define inputs
kin.file <- "G4_kin_sum.csv"
exp.conc.file <- "G4_homkin_evn_concs.csv"
g4.diss.file <- "G4_Homarine_Dissolved_Conc.csv"
g4.part.file <- "G4_Homarine_Particulate_Quant.csv"



###
kin.dat <- read_csv(kin.file)
exp.part.dat <- read_csv(exp.conc.file)
g4.diss.dat <- read_csv(g4.diss.file)
g4.part.dat <- read_csv(g4.part.file)



###
g4.ukh.diss <- g4.diss.dat %>%
  filter(str_detect(Samp_ID, "UKH")) %>%
  separate(Samp_ID, into = c("cruise", "exp", "treatment", "rep")) %>%
  group_by(exp) %>%
  mutate(mean.conc = mean(Conc.nM))

dat.a <- left_join(exp.part.dat, g4.ukh.diss) %>%
  mutate(perc.diss = mean.conc/(env_nM+mean.conc)) %>%
  select(exp, env_nM, mean.conc, perc.diss) %>%
  rename("part.conc" = env_nM,
         "diss.conc" = mean.conc) %>%
  unique() %>%
  pivot_longer(cols = part.conc:perc.diss, names_to = "param", values_to = "val")

ggplot(dat.a, aes(x = exp, y = val, fill = exp)) +
  geom_col(width = 0.65) +
  facet_grid(param~., scales = "free")

dat.c <- dat.a %>%
  filter(!param == "perc.diss")


ggplot(dat.c, aes(x = exp, y = val, fill = param)) +
  geom_col(width = 0.65, position = "dodge")

dat.c <- dat.a %>%
  filter(!param == "perc.diss")

##compare to kin params
kin.wide <- kin.dat %>%
  pivot_wider(names_from = parameter, values_from = Estimate) %>%
  rename("exp" = experiment)

dat.b <- left_join(exp.part.dat, g4.ukh.diss) %>%
  mutate(ratio = mean.conc/env_nM) %>%
  select(exp, env_nM, mean.conc, ratio) %>%
  rename("part.conc" = env_nM,
         "diss.conc" = mean.conc) %>%
  unique() %>%
  left_join(., kin.wide) %>%
  mutate(flux = (diss.conc/turnover_time)*24*7,
         perc_diss = diss.conc/(part.conc+diss.conc))

ggplot(dat.b, aes(x = diss.conc, y = Ks)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 2) +
  theme_bw() +
  xlab("Dissolved Concentration (nM)") +
  ylab("Community Affinity (Ks) (nM)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) 
ggsave(filename = "lin_plot.png", height = 4, width = 5)





lm.ks <- lm(Ks ~ diss.conc, data = dat.b)
summary(lm.ks)

##Flux
ggplot(dat.b, aes(x = exp, y = flux, fill = exp)) +
  geom_col(width = 0.6, color = "black", size = 0.1) +
  theme_bw() +
  ylab("Carbon Flux (nM C/day)") +
  xlab("Experiment")
ggsave(filename = "fluxfig.png", width = 5, height = 5)

##Turnover Time
ggplot(dat.b, aes(x = exp, y = turnover_time, fill = exp)) +
  geom_col(width = 0.6, color = "black", size = 0.1) +
  theme_bw() +
  ylab("Turnover Time (hr)") +
  xlab("Experiment")
ggsave(filename = "TTfig.png", width = 5, height = 5)

##Dissolved Concentration
ggplot(dat.b, aes(x = exp, y = diss.conc, fill = exp)) +
  geom_col(width = 0.6, color = "black", size = 0.1) +
  theme_bw() +
  ylab("Dissolved Concentration (nM)") +
  xlab("Experiment")

ggsave(filename = "Dissconcfig.png", width = 5, height = 5)




ggplot(dat.b, aes(x = exp, y = turnovertime, fill = exp)) +
  geom_col(width = 0.7) +
  ylab("Carbon Flux (nM C/day)")



ggplot(dat.b, aes(x = exp, y = perc_diss, fill = exp)) +
  geom_col(width = 0.7)






