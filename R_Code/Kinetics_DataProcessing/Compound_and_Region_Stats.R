



#load packages:
library(tidyverse)
library(rstatix)


#define inputs:
flux.kin.file <- "Intermediates/Compiled_Kinetics_Fluxes.csv"
context.file <- "Intermediates/flux_contextualizaing_data"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"
enviro.metab.file <- "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv"


#Perform nonparametric statistics on means of Kt and Vmax comparing GBT and Homarine
flux.kin.dat <- read_csv(flux.kin.file)

####Perform wilcoxon rank sum test on mean Kt values:
kt.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mean_ks ~ Compound)
kt.comp.comparison

vmax.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mean_vmax ~ Compound)
vmax.comp.comparison

tt.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mm_tt ~ Compound)
tt.comp.comparison

flux.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mm_flux_nM_day ~ Compound)
flux.comp.comparison

diss.conc.comp.comparison <- flux.kin.dat %>%
  wilcox_test(Mean.Diss.Conc.nM ~ Compound)




###Estimating free living bacteria as a % of POC and comparing to Kt:
#meta.dat <- read_csv(meta.data.file) %>%
#  rename("cruise_exp" = KinExp_ID)

#kin.flux.meta.dat <- left_join(meta.dat, flux.kin.dat) %>%
#  mutate(Bact.Perc.POC = bact_c/pc) 

#ggplot(kin.flux.meta.dat, aes(x = Bact.Perc.POC, y = mean_ks, color = Compound)) + 
#  geom_point() #+
 ## scale_y_log10() +
 # scale_x_log10()
  



#linear model of Flux vs. Part Conc.
context.dat <- read_csv(context.file)

ggplot(context.dat, aes(y = mm_flux_nM_day, x = mean.part.conc.nM)) +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", aes(color = Compound)) +
  geom_point(aes(color = Compound), size =2 ) +
  scale_y_log10() +
  scale_x_log10() 

log.flux.pconc <- context.dat %>%
  select(KinExp_ID, Cruise, exp, Compound, mm_flux_nM_day, mean.part.conc.nM) %>%
  mutate(log10_flux = log10(mm_flux_nM_day),
         log10_pconc = log10(mean.part.conc.nM))


flux.pconc.model <- lm(log10_flux~log10_pconc, data = log.flux.pconc)
summary(flux.pconc.model)

##Just homarine:
hom.flux.pconc.dat <- log.flux.pconc %>%
  filter(Compound == "Homarine")

hom.flux.pconc.model <- lm(log10_flux~log10_pconc, data = hom.flux.pconc.dat)
summary(hom.flux.pconc.model)

##Just GBT:
gbt.flux.pconc.dat <- log.flux.pconc %>%
  filter(Compound == "GBT")

gbt.flux.pconc.model <- lm(log10_flux~log10_pconc, data = gbt.flux.pconc.dat)
summary(gbt.flux.pconc.model)




#Perform regression of Kt vs. dissolved concentration 

#plot data of kt vs. dissolved concentration:
ggplot(flux.kin.dat, aes(x = Mean.Diss.Conc.nM, y = mean_ks)) +
  geom_point(aes(color = Compound)) +
 # facet_wrap(.~Compound, scales = "free") +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1)

###############

log.kt.dconc.dat <- flux.kin.dat %>%
  mutate(log10_ks = log10(mean_ks),
         log10_dconc = log10(Mean.Diss.Conc.nM))


###Overall model 
kt.dconc.lm <- lm(log10_dconc ~ log10_ks, data = log.kt.dconc.dat)
summary(kt.dconc.lm)


###Just homarine
hom.kt.dconc.dat <- log.kt.dconc.dat %>%
  filter(Compound == "Homarine")

hom.kt.dconc.lm <- lm(log10_dconc ~ log10_ks, data = hom.kt.dconc.dat)
summary(hom.kt.dconc.lm)


###Just GBT
gbt.kt.dconc.dat <- log.kt.dconc.dat %>%
  filter(Compound == "GBT")

gbt.kt.dconc.lm <- lm(log10_dconc ~ log10_ks, data = gbt.kt.dconc.dat)
summary(gbt.kt.dconc.lm)





####Look at all betaines + sulfoniums 
#meta.dat <- read_csv(meta.data.file)
region.dat <- meta.dat %>%
  select(KinExp_ID, Region) %>%
  separate(KinExp_ID, into = c("Cruise", "exp"))


enviro.metab.dat <- read_csv(enviro.metab.file)

##Remove TMAO and calculate total betaine + sulfonium concentration:
summed.comp.dat <- enviro.metab.dat %>%
  filter(!Compound == "Trimethylamine N-oxide") %>%
  group_by(Cruise, exp) %>%
  reframe(metab.tot.nM = sum(Mean.Diss.Conc.nM),
          metab.tot.nM.uncert = sqrt(sum(SD.Diss.Conc.nM)^2))

##join with kinetics data:
kin.tot.metab.dat <- left_join(summed.comp.dat, flux.kin.dat) %>%
  left_join(., region.dat)

#### Kt-divergence vs. % of total pool data
kt.conc.dat <- kin.tot.metab.dat %>%
  mutate(Kt_diff = (mean_ks - Mean.Diss.Conc.nM)/(mean_ks)*100,
         perc_diss_pool = (Mean.Diss.Conc.nM/metab.tot.nM)*100) %>%
  select(Cruise, exp, Compound, Kt_diff, perc_diss_pool)

ggplot(kt.conc.dat, aes(x = perc_diss_pool, y = Kt_diff)) +
  geom_smooth(method = "lm") +
  geom_point(size = 3, aes(color = Compound)) 
 # facet_wrap(.~Compound, scales = "free")#+
 # scale_x_log10() +
 # scale_y_log10()


kt.diff.lm <- lm(Kt_diff ~ perc_diss_pool, data = kt.conc.dat)
summary(kt.diff.lm)


ggplot(kin.tot.metab.dat, aes(x = metab.tot.nM, y = mean_ks)) +
  geom_smooth(alpha = 0.5, method = "lm") +
  geom_errorbar(aes(ymin = mean_ks-sd_ks, ymax = mean_ks+sd_ks), width = 0.03) +
  geom_errorbarh(aes(xmin = metab.tot.nM-metab.tot.nM.uncert, xmax = metab.tot.nM+metab.tot.nM.uncert), height = 0.03) +
  geom_point(aes(shape = Compound, fill = Region), size = 3) +
  scale_shape_manual(values = c(21, 23)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1) 


####Calculate log space 
kin.tot.log.dat <- kin.tot.metab.dat %>%
  mutate(log10_totmetab = log10(metab.tot.nM),
         log10_ks = log10(mean_ks))

##overall model
kin.totmetab.lm <- lm(log10_ks~log10_totmetab, data = kin.tot.log.dat)
summary(kin.totmetab.lm)

##just homarine model
hom.kin.totmetab.dat <- kin.tot.log.dat %>%
  filter(Compound == "Homarine")

hom.kin.totmetab.lm <- lm(log10_ks~log10_totmetab, data = hom.kin.totmetab.dat)
summary(hom.kin.totmetab.lm)

##just GBT model
gbt.kin.totmetab.dat <- kin.tot.log.dat %>%
  filter(Compound == "GBT")

gbt.kin.totmetab.lm <- lm(log10_ks~log10_totmetab, data = gbt.kin.totmetab.dat)
summary(gbt.kin.totmetab.lm)





#######Check on TT vs. %D
part.dat <- context.dat %>%
  select(Cruise, exp, KinExp_ID, mean.part.conc.nM)


diss.tt.dat <- flux.kin.dat %>%
  select(Cruise, exp, Compound, Mean.Diss.Conc.nM, mm_tt)

part.diss.tt.dat <- left_join(part.dat, diss.tt.dat) %>%
  mutate(Perc.Diss = Mean.Diss.Conc.nM/(Mean.Diss.Conc.nM+mean.part.conc.nM))

ggplot(part.diss.tt.dat, aes(x = mean.part.conc.nM, y = mm_tt, color = Compound)) +
  geom_point() #+ 
#  scale_y_log10()


###Comparisons of WH_TT and MM_TTs
flux.kin.dat

ggplot(flux.kin.dat, aes(x = mm_tt, y = wh_mean_tt)) + 
  geom_point()

#model of TTs
tt.lm <- lm(mm_tt~wh_mean_tt, data = flux.kin.dat)
summary(tt.lm)

cor.val <- cor(flux.kin.dat$mm_tt, flux.kin.dat$wh_mean_tt)
cor.val




























pelagic.only.dat <- flux.kin.dat %>%
  filter(Cruise %in% c("KM1906", "TN397", "TN412"))

ggplot(pelagic.only.dat, aes(x = Mean.Diss.Conc.nM, y = mean_ks)) +
  geom_point() +
  facet_wrap(.~Compound, scales = "free") +
  geom_smooth(method = "lm")


kt.diss.model <- 


############

###Contextualizaing data:
context.dat <- read_csv(context.file)

ggplot(context.dat, aes(y = mm_flux_nM_day, x = mean.part.conc.nM, color = Compound)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

ggplot(contex.dat, aes(y=mm_flux_nM_day, x = ))
  

#































































