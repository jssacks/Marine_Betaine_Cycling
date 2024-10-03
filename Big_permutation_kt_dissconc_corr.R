



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
region.dat <- meta.dat %>%
  select(cruise_exp, Region) %>%
  separate(cruise_exp, into = c("Cruise", "exp"))


enviro.metab.dat <- read_csv(enviro.metab.file)

##Remove TMAO and calculate total betaine + sulfonium concentration:
summed.comp.dat <- enviro.metab.dat %>%
  filter(!Compound == "Trimethylamine N-oxide") %>%
  group_by(Cruise, exp) %>%
  reframe(metab.tot.nM = sum(Mean.Diss.Conc.nM),
          metab.tot.nM.uncert = sqrt(sum(SD.Diss.Conc.nM)^2))


#####


dis.dat <- enviro.metab.dat %>%
  select(-SD.Diss.Conc.nM) 

dis.dat.add <- enviro.metab.dat %>%
  select(-SD.Diss.Conc.nM) %>%
  rename(Compound.add = Compound,
         Mean.Conc.add = Mean.Diss.Conc.nM)

variable.combination.diss.dat <- dis.dat %>%
  left_join(., dis.dat.add) %>%
  left_join(., dis.dat.add) %>%
  left_join(., dis.dat.add) 


######Make combination list:
comp.list <- enviro.metab.dat %>%
  select(Compound) %>%
  filter(!Compound == "Trimethylamine N-oxide") %>%
  unique()

all_combinations <- lapply(1:length(comp.list$Compound), function(size) combn(comp.list$Compound, size, simplify = FALSE))

# Flatten the list of lists into a single list
all_combinations_flat <- unlist(all_combinations, recursive = FALSE)

# Find the maximum combination size to handle varying lengths
max_size <- max(sapply(all_combinations_flat, length))

# Convert each combination into a data frame row and fill missing values with NA
combinations_df <- do.call(rbind, lapply(all_combinations_flat, function(x) {
  c(x, rep(NA, max_size - length(x)))
}))

# Convert to data.frame and set column names
combinations_df <- as.data.frame(combinations_df, stringsAsFactors = FALSE)
colnames(combinations_df) <- paste0("Var", 1:ncol(combinations_df))

###########Make data frame long format and add in concentrations

comb_df_long <- combinations_df %>%
  rownames_to_column(var = "rep") %>%
  pivot_longer(!rep, names_to = "var", values_to = "Compound") %>%
  filter(!is.na(Compound)) %>%
  select(-var) %>%
  group_by(rep) %>%
  mutate(numb_comps = n()) %>%
  ungroup()

###make dataframe with enviro conc for each ____
enviro.conc.dat <- enviro.metab.dat %>%
  select(Cruise, exp, Compound, Mean.Diss.Conc.nM) %>%
  unique()

##make dataframe of each compound to rbind as a "just one compound" dataset:
#enviro.indivd.comp <- enviro.conc.dat %>%
#         rep = "ind")

####combine dataframes for iterations and concentration
comb_df_long_concs <- comb_df_long %>%
  left_join(., enviro.conc.dat) 


###
comb_df_long_sum <- comb_df_long_concs %>%
  group_by(rep, Cruise, exp, numb_comps) %>%
  reframe(conc_sum_nM = sum(Mean.Diss.Conc.nM)) %>%
  left_join(., flux.kin.dat %>%
              select(Cruise, exp, Compound, mean_ks))

####Run correlation analyses:
comb_long_cor <- comb_df_long_sum %>%
  mutate(log10_conc = log10(conc_sum_nM),
         log10_Kt = log10(mean_ks)) %>%
  group_by(rep) %>%
  mutate(cor = cor(x = log10_conc, y = log10_Kt)) %>%
  left_join(., summed.comp.dat) %>%
  mutate(perc_pool = (conc_sum_nM/metab.tot.nM)*100) %>%
  group_by(rep) %>%
  mutate(mean_perc_pool = mean(perc_pool))

comb_fig_dat <- comb_long_cor %>%
  select(rep, numb_comps, cor, mean_perc_pool) %>%
  unique()

ggplot(comb_fig_dat, aes(x = mean_perc_pool, y = cor)) +
  geom_jitter(alpha = 0.4, width = 0.2, size = 2) + 
  theme_bw()

ggplot(comb_fig_dat, aes(x = numb_comps, y = cor)) +
  geom_jitter(alpha = 0.4, width = 0.2, size = 2) + 
  theme_bw()

##total sum
#





comp.list <- enviro.metab.dat %>%
  select(Compound) %>%
  unique() #%>%
  expand_grid(comp.2 = Compound,
              comp.3 = Compound,
              comp.4 = Compound,
              comp.5 = Compound,
              comp.6 = Compound,
              comp.7 = Compound) %>%
  filter(!Compound == comp.2,
         !Compound == comp.3,
         !Compound == comp.4,
         !Compound == comp.5,
         !Compound == comp.6,
         !Compound == comp.7,
         !comp.2 == comp.3)

install.packages("RcppAlgos")
library(RcppAlgos)

x.test <- data.frame(comboGeneral(comp.list$Compound, length(comp.list$Compound), repetition = TRUE)) %>%
  crossing()


comb.list <- comp.list %>%
  expand(crossing(comp.listCompound))

test.2 <- data.frame(combn(comp.list$Compound, 8))

##join with kinetics data:
kin.tot.metab.dat <- left_join(summed.comp.dat, flux.kin.dat) %>%
  left_join(., region.dat)

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
