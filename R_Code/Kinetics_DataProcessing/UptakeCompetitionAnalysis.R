

#packages
library(tidyverse)
library(lubridate)
library(broom)
library(rstatix)
library(DescTools)

source("R_Code/Kinetics_DataProcessing/Functions.R")

###inputs
dat.rate <- read_csv("Intermediates/quantified_uptake_rates.csv")




####Calculate Michelis Mentin Uptake Kinetic Parameters by Fitting Nonlinear Equations to Uptake Curves

##Subset data 
dat.comp <- dat.rate %>%
  filter(exp == "UCH1") %>%
  filter(!str_detect(SampID, "Blk")) %>%
  select(SampID, cruise, exp, treatment, rep, Fragment_mz, nM.per.hour.1) %>%
  filter(Fragment_mz == 78.0399) %>%
  filter(!str_detect(SampID, "TN397_UCH1_Hom")) %>%
  group_by(cruise, treatment) %>%
  mutate(mean.rate = mean(nM.per.hour.1))

ggplot(dat.comp, aes(x = treatment, y = nM.per.hour.1)) +
  geom_point() +
  facet_wrap(.~cruise, scales = "free_x")

#########create ANOVA for each experiment:

#TN397:
tn397.dat <- dat.comp %>%
  filter(cruise == "TN397")
tn397.model <- aov(nM.per.hour.1 ~ treatment, data = tn397.dat)
summary(tn397.model)

tn397.dt <- DunnettTest(nM.per.hour.1 ~ treatment, data = tn397.dat, control = "Gluc")
tn397.dt

tn397.output <- data.frame(tn397.dt[["Gluc"]]) %>%
  rownames_to_column(.,) %>%
  separate("rowname", into = c("treatment", "control")) %>%
  mutate(cruise = "TN397")



#RC104:
rc104.dat <- dat.comp %>%
  filter(cruise == "RC104")
rc104.model <- aov(nM.per.hour.1 ~ treatment, data = rc104.dat)
summary(rc104.model)

rc104.dt <- DunnettTest(nM.per.hour.1 ~ treatment, data = rc104.dat, control = "Gluc")
rc104.dt

rc104.output <- data.frame(rc104.dt[["Gluc"]]) %>%
  rownames_to_column(.,) %>%
  separate("rowname", into = c("treatment", "control")) %>%
  mutate(cruise = "RC104")


###Combined output:
all.stat.output <- rbind(tn397.output, rc104.output) %>%
  select(treatment, pval, cruise)

###combine stat data with data to plot
dat.fig <- dat.comp %>%
  left_join(., all.stat.output) %>%
  mutate(signif = case_when(pval < 0.05 ~ "Significant (p < 0.05)",
                            pval > 0.05 ~ "Not Significant (p > 0.05)",
                            is.na(pval) ~ "Glucose Control"),
        order = case_when(pval < 0.05 ~ 3,
                          pval > 0.05 ~ 2,
                          is.na(pval) ~ 1))

write_csv(dat.fig, file = "Intermediates/Uptake_Competition_Analysis_Output.csv")

####XXX####

#dat.perc.of.control <- dat.comp %>%
#  group_by(cruise, exp) %>%
#  mutate(perc.of.control = mean.rate/mean.rate[treatment == "Gluc"])
  
####

ggplot(dat.fig, aes(x = treatment, y = nM.per.hour.1)) +
  geom_col(aes(x = reorder(treatment, order), y = mean.rate, fill = signif), color = "black", size = 0.1, width = 0.6, position = "dodge") +
  scale_fill_manual(values = c("lightgray", "#71b05d" , "steelblue3")) +
  geom_point() +
  theme_test() +
  facet_grid(.~cruise, scales = "free_x", space = "free") +
  xlab("Treatment") +
  ylab("Uptake Rate (nM/hr)") +
  labs(fill = "") + 
  theme(legend.position = "bottom")
ggsave(filename = "Comp_fig.png", height = 4, width = 7, units = "in")















####Run test in separate groups, arrange groups by "control" type compound, 
# run separate fdr correction on only comparisons we are interested in 
# join together

#dat.t <- dat.comp %>%
#  arrange(desc(nM.per.hour.1)) %>%
#  group_by(cruise) %>%
#  t_test(nM.per.hour.1~treatment, ref.group = "Gluc")


#dat.gluc <- data.frame(
#  cruise = c("RC104", "TN397"),
#  treatment = c("Gluc", "Gluc"),
#  signif = c("not significant", "not significant"))

#dat.t.res <- dat.t %>%
#  select(cruise, group2, p.adj) %>%
#  mutate(signif = case_when(p.adj > 0.1 ~ "not significant",
#                            p.adj <= 0.1 ~ "significant")) %>%
#  rename("treatment" = group2) %>%
#  select(-p.adj) %>%
#  rbind(dat.gluc)







