
#packages
library(tidyverse)
library(lubridate)
library(broom)

source("R_Code/Kinetics_DataProcessing/Functions.R")

###inputs
dat.rate <- read_csv("Intermediates/quantified_uptake_rates.csv")
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"





###Pull in GBT uptake rate data from Gradients 3 cruise from Boysen et al. 2022
boysen.dat <- read_csv("Meta_Data/Data_From_Other_Studies/Boysen2022_GBT_uptake_rates.csv") %>%
  select(Replicate.Name, Treatment, replicate, correct.Prediction.nmoles.per.L.per.hr) %>%
  rename("treatment_conc" = Treatment,
         "nM.per.hour.1" = correct.Prediction.nmoles.per.L.per.hr) %>%
  mutate(cruise = "KM1906") %>%
  mutate(exp = case_when(str_detect(Replicate.Name, "K1") ~ "UKG1",
                         str_detect(Replicate.Name, "K2") ~ "UKG2")) %>%
  filter(replicate %in% c("A", "B", "C")) %>%
  rename("rep" = replicate) %>%
  mutate("Fragment_mz" = 62.8977) %>%
  select(-Replicate.Name)


####Calculate Michelis Mentin Uptake Kinetic Parameters by Fitting Nonlinear Equations to Uptake Curves

##Subset data 
dat.MM <- dat.rate %>%
  filter(!exp == "UCH1") %>%
  filter(!str_detect(SampID, "Blk")) %>%
  select(cruise, exp, treatment_conc, rep, Fragment_mz, nM.per.hour.1) %>%
  rbind(boysen.dat)

###Visualize data to identify starting parameters 
dat.MM.prelimviz <- dat.MM %>%
  group_by(cruise, exp) %>%
  mutate(max.rate = max(nM.per.hour.1, na.rm = TRUE))


###Examine raw data to define starting predictions for Vmax (a) and Ks (b) 
# the red line is the max uptake rate and likely represents a good prediction for Vmax,
ggplot(dat.MM.prelimviz, aes(x = treatment_conc, y = nM.per.hour.1)) +
  geom_point() +
  facet_wrap(cruise~exp, scales = "free") +
  geom_hline(aes(yintercept = max.rate), color = "red") 

##Inspect plots and manually enter predicted values into dataframe, 
# enter values in order of cruise, (cruise), experiment (exp) for a_pred (Vmax prediction) and b_pred (Ks prediction)
model.starting.estiamtes <- data.frame(
  cruise = c("KM1906","KM1906", "RC078", "RC078", "RC078", "RC104", "RC104", "TN397", "TN397", "TN397", "TN397", "TN397", "TN412"),
  exp = c("UKG1", "UKG2", "UKG1", "UKH1", "UKH2", "UKG1", "UKH1", "UKH1", "UKH2", "UKH3", "UKH4", "UKH5", "UKG"),
  a_pred = c(300, 10, 400, 200, 350, 200, 600, 125, 200, 300, 500, 125, 125),
  b_pred = c(0.4, 0.6, 4.1, 0.15, 0.4, 0.6, 0.35, 0.1, 0.1, 0.2, 0.3, 0.6, 0.3)
)

mod.start.estimates.MM.nls <- model.starting.estiamtes %>%
  unite(col = "exp", cruise:exp)

###
dat.MM.nls <- dat.MM %>%
  rename("treatment_nM" = treatment_conc,
         "nM_per_hour" = nM.per.hour.1)  %>%
  unite(col = "exp", cruise:exp)

#dat.test <- dat.MM.nls %>%
#  filter(exp.nls == "TN397_UKH5")

dat.exp.names <- dat.MM.nls %>%
  select(exp) %>%
  rename("exp.nls" = exp) %>%
  unique()

#out.test <- fit_nls(dat.test, mod.start.estimates.MM.nls, "TN397_UKH5")

MM.nls.out <- data.frame(
  experiment = as.character(),
  parameter = as.character(),
  Estimate = as.numeric()
)

for (i in dat.exp.names$exp.nls) {
  x <- fit_nls(dat.MM.nls, mod.start.estimates.MM.nls, i)
  MM.nls.out <- rbind(MM.nls.out, x)
}

final.nls.output <- MM.nls.out

out.hom <- final.nls.output %>%
  filter(str_detect(experiment, "UKH"))

ggplot(out.hom, aes(x = experiment, y = Estimate, fill = experiment)) +
  geom_col(width = 0.6) +
  facet_wrap(.~parameter, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))
###

out.gbt <- final.nls.output %>%
  filter(str_detect(experiment, "UKG"))

ggplot(out.gbt, aes(x = experiment, y = Estimate, fill = experiment)) +
  geom_col(width = 0.6) +
  facet_wrap(.~parameter, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

###
mm.plot <- final.nls.output %>%
  pivot_wider(names_from = parameter, values_from = Estimate) %>%
 # filter(Vmax < 2) %>%
  mutate(Compound = case_when(str_detect(experiment, "UKG") ~ "GBT",
                              str_detect(experiment, "UKH") ~ "Homarine")) %>%
  mutate(Environment = case_when(str_detect(experiment, "RC104") ~ "coastal",
                                 str_detect(experiment, "RC078") ~ "coastal",
                                 str_detect(experiment, "TN412") ~ "Gyre",
                                 str_detect(experiment, "KM1906") ~ "NPTZ",
                                 str_detect(experiment, "UKH4") ~ "Equatorial",
                                 str_detect(experiment, "UKH3") ~ "Equatorial",
                                 str_detect(experiment, "UKH5") ~ "Equatorial",
                                 str_detect(experiment, "UKH2") ~ "Gyre",
                                 str_detect(experiment, "UKH1") ~ "Gyre"))



#install.packages("ggbreak")
library(ggbreak)

ggplot(mm.plot, aes(x = Ks, y = Vmax, shape = Compound, fill = Environment)) +
  geom_point(size = 4.5) +
  ylim(0, 4.00) +
  scale_y_break(c(0.8, 3.8), scales = "0.2", ticklabels = c(3.8, 4.0)) +
  scale_shape_manual(values = c(21, 24)) +
 # scale_fill_wsj() +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("Community Affinity (Km, nM)") +
  ylab("Maximum Uptake Rate (Vmax, nM/hr)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

ggsave(filename = "Param_fig.png", height = 7, width = 7)


ggplot(mm.plot, aes(x = Ks, y = Vmax, shape = Compound, fill = Environment)) +
  geom_point(size = 4) +
 # ylim(0, 4.00) +
  scale_x_log10() +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 24)) +
  # scale_fill_wsj() +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) 



###########

nls.fun.test <- function(x) (0.4*(x))/(x+200)

nls.plot <- final.nls.output %>%
  pivot_wider(names_from = parameter, values_from = Estimate) %>%
  rename("exp" = experiment) %>%
  mutate(join.val = "x")



##Create a dataframe to plot the nonlinear equations with
dat.calc <- data.frame(treatment_nM = seq(1, 6000, by = 1)) %>%
  mutate(join.val = "x")

dat.max.treat <- dat.MM.nls %>%
  group_by(exp) %>%
  mutate(max.treat = max(treatment_nM)) %>%
  select(exp, max.treat)

t.1 <- full_join(dat.calc, nls.plot, relationship = "many-to-many") %>%
  mutate(nM_per_hour = (Vmax*treatment_nM)/(Ks+treatment_nM)) %>%
  left_join(., dat.max.treat, relationship = "many-to-many") %>%
  filter(treatment_nM <= max.treat) %>%
  unique() %>%
  left_join(., mm.plot)

t.2 <- full_join(dat.calc, nls.plot, relationship = "many-to-many") %>%
  mutate(nM_per_hour = (Vmax*treatment_nM)/(Ks+treatment_nM)) %>%
  left_join(., dat.max.treat, relationship = "many-to-many") %>%
  unique() %>%
  left_join(., mm.plot) %>%
  filter(treatment_nM <= 1000)

region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7")


###Making many Kinetics curves figure for Poster:
match.dat.t2 <- read_csv(match.file) %>%
  mutate("exp" = KinExp_ID) %>% 
  select(-Compound) %>% 
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "Gyre", "Equatorial", "Coastal"))) %>%
  right_join(., t.2)

kin.curves.plot <- ggplot(match.dat.t2, aes(x = treatment_nM, y = nM_per_hour, color = Region, group = exp)) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(Compound~.) +
  ylim(0, 4.00) +
  scale_color_manual(values = region.pal) +
 # scale_color_economist(economist_pal(fill = TRUE)(7)) +
  theme_test() +
  scale_y_break(c(0.7, 2.5), scales = "0.2", ticklabels = c(2.5, 3.0, 3.5, 4.0)) +
  ylab("Uptake Rate (nM/hr)") +
  xlab("Spike Concentration (nM)")
  
kin.curves.plot
ggsave(kin.curves.plot, width = 7, height = 4.5, scale = 1.1, units = "in", file = "Figures/Kin_Curves_Plot.jpg")
######


dat.plot <- left_join(dat.MM.nls, nls.plot) %>%
  left_join(., mm.plot) %>%
  filter(Fragment_mz %in% c(97.0954, 62.8977))

library(ggthemes)

ggplot(dat.plot, aes(y = nM_per_hour, x = treatment_nM)) +
  theme_classic() +
  scale_shape_manual(values = c(21, 24)) +
  geom_point(size = 3, aes(shape = Compound, fill = Environment)) +
  facet_wrap(.~exp, scales = "free") + 
#  scale_fill_few() +
  geom_line(data = t.1, aes(x = treatment_nM, y = nM_per_hour, color = Compound)) +
  scale_color_manual(values = c("red", "blue")) +  
  guides(fill = guide_legend(override.aes = list(shape = 21)),
  shape = guide_legend(override.aes = list(fill = "black"))) +
  theme(legend.position = c(0.6, 0.09),
       # legend.box = "horizontal",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)) +
  xlab("Spike Concentration (nM)") +
  ylab("Uptake Rate (nM/hr)")  +
  coord_cartesian(clip = "off")
ggsave(filename = "all_kinetics_exps.png", height = 6, width = 9, units = "in")




###Turnover time
  
###Thinkin about affinity as Vmax/Km 
mm.plot.2 <- mm.plot %>%
  mutate(Af = Vmax/Ks,
         Af2 = Ks/Vmax)
ggplot(mm.plot.2, aes(x = experiment, y = Af, fill = Environment)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.6) +
  facet_wrap(.~Compound, scales = "free")

#ggplot(mm.plot.2, aes(x = experiment, y = Af2, fill = Environment)) +
#  geom_col(color = "black", linewidth = 0.1, width = 0.6) +
#  facet_wrap(.~Compound, scales = "free")



################### Turnover time calculations
dat.TT <- dat.MM.nls %>%
  mutate(time_over_fract = 1/(nM_per_hour/treatment_nM),
         log10.tof = log10(time_over_fract),
         log10.t_nM = log10(treatment_nM)) %>%
  group_by(Fragment_mz, exp) %>%
  filter(!treatment_nM == 0)

ggplot(dat.TT, aes(x = treatment_nM, y = time_over_fract)) +
  geom_point() +
  facet_wrap(Fragment_mz ~ exp, scales = "free")

###perform linear models across groupings

tt.lm.out <- do(dat.TT,
             tidy(
               lm(log10.tof ~ log10.t_nM, data = .)))
###pull out, tidy, and rename values for standard curve slopes
#lm.out.slopes <- tt.lm.out %>%
#  filter(term == "spike_conc") %>%
#  select(cruise, exp, Fragment_mz, estimate, std.error) %>%
#  rename("slope" = estimate,
#         "std.error.slope" = std.error)

###pull out, tidy, and rename values for standard curve intercepts
tt.lm.out.intercepts  <- tt.lm.out %>%
  filter(term == "(Intercept)") %>%
  select(exp, Fragment_mz, estimate, std.error) %>%
  rename("intercept" = estimate,
         "std.error.intercept" = std.error) %>%
  mutate(TT = 10^intercept)

tt.intercept.fig.dat <- tt.lm.out.intercepts %>%
  filter(Fragment_mz %in% c(62.8977, 78.0399))

###Write to csv:
write_csv(tt.intercept.fig.dat, file = "Intermediates/turnover_time_estimates.csv")




ggplot(tt.intercept.fig.dat, aes(x = exp, y = TT)) +
  geom_col(width = 0.5) + 
  facet_wrap(.~Fragment_mz, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))


gbt.plot.dat <- dat.TT %>%
  filter(Fragment_mz %in% c(62.8977))

ggplot(gbt.plot.dat, aes(x = log10.t_nM, y = log10.tof)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm") +
  facet_wrap(.~exp)

ggplot(gbt.plot.dat, aes(x = treatment_nM, y = time_over_fract)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm") +
  facet_wrap(.~exp, scales = "free")

