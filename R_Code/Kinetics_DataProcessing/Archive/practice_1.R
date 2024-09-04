



library(tidyverse)
library(lubridate)
library(broom)

source("Functions.R")


#dat.1 <- read_csv("G4_UKH4_UCH1_TQS_results.csv") %>%
#  rename(SampID = "Replicate Name",
#         Compound = "Precursor Ion Name",
#         Fragment_mz = "Product Mz") %>%
#  select(SampID, Compound, Area, Fragment_mz)

dat.1 <- TQS_csv_read("G4_UKH4_UCH1_TQS_results.csv")
dat.2 <- TQS_csv_read("G4_UKH5_TQS_results.csv")
dat.3 <- TQS_csv_read("G4_Kinetics_TQS_UKH123.csv")
dat.4 <- read_csv("TQS_HomarineFocus_SKylineReport.csv")
dat.4 <- TQS_csv_read("TQS_HomarineFocus_SKylineReport.csv")

dat.1.1 <- rbind(dat.1, dat.2)



###Wrangle Sample Data
dat.1.smp <- dat.1.1 %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(Compound == "D3-Homarine")

dat.3.smp <- dat.3 %>%
  filter(!str_detect(SampID, "inM")) %>%
  filter(!str_detect(SampID, "Std")) %>%
  filter(Compound == "D3-Homarine") %>%
  mutate(SampID = str_replace(SampID, "808_U", "808_Smp_TN397_U")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) 

dat.4.smp <- dat.4 %>%
  mutate(SampID = str_replace(SampID, "808_U", "808_Smp_TN397_U")) %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(Compound == "D3-Homarine")

dat.smp <- rbind(dat.1.smp, dat.3.smp)


####Wrangle Standard Curve Data 
dat.1.std <- dat.1.1 %>%
  filter(str_detect(SampID, "Std")) %>%
  filter(!str_detect(SampID, "In")) %>%
  filter(!str_detect(SampID, "inM")) %>%
  mutate(SampID = str_replace(SampID, "UHC1", "UCH1")) %>%
  separate(col = SampID, into = c("date", "type", "spike_conc", "hom", "exp", "cal_curve"), remove = FALSE) %>%
  filter(Compound == "D3-Homarine") %>%
  mutate(spike_conc = as.numeric(str_remove(spike_conc, "nM")))

dat.3.std <- dat.3 %>%
  filter(str_detect(SampID, "Std")) %>%
  filter(!str_detect(SampID, "In")) %>%
  filter(!str_detect(SampID, "inM")) %>%
  separate(col = SampID, into = c("date", "type", "spike_conc", "hom", "exp", "cal_curve"), remove = FALSE) %>%
  filter(Compound == "D3-Homarine") %>%
  mutate(spike_conc = as.numeric(str_remove(spike_conc, "nM")))

dat.std <- rbind(dat.1.std, dat.3.std)






ggplot(dat.std, aes(x = spike_conc, y = Area, color = exp)) +
  geom_point() +
  geom_smooth(method = lm, alpha = 0.5) +
  facet_wrap(. ~ Fragment_mz, scales = "free")

ggplot(dat.std, aes(x = spike_conc, y = Area, color = exp)) +
  geom_point() +
  geom_smooth(method = lm, alpha = 0.5) +
  facet_grid(exp ~Fragment_mz, scales = "free_x")

ggplot(dat.std, aes(x = reorder(SampID, spike_conc), y = Area)) +
  geom_col() +
  facet_grid(Fragment_mz~exp, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))


dat.cal.curve <- dat.std %>%
  ungroup() %>%
  group_by(exp, Fragment_mz) 


lm.test <- lm(Area ~ spike_conc, data = dat.cal.curve)
summary(lm.test)

lm.test[["coefficients"]][[2]]

#####
dat.cal.curve <- dat.std %>%
  ungroup() %>%
  group_by(exp, Fragment_mz) 

lm.out <- do(dat.cal.curve,
   tidy(
     lm(Area ~ spike_conc, data = .)))

lm.out.slopes <- lm.out %>%
  filter(term == "spike_conc") %>%
  select(exp, Fragment_mz, estimate) %>%
  rename("slope" = estimate)

lm.out.intercepts  <- lm.out %>%
  filter(term == "(Intercept)") %>%
  select(exp, Fragment_mz, estimate) %>%
  rename("intercept" = estimate)

ggplot(lm.out.slopes, aes(x = exp, y = slope, fill = exp)) +
  geom_col(width = 0.7) +
  facet_wrap(.~Fragment_mz)

ggplot(lm.out.intercepts, aes(x = exp, y = intercept, fill = exp)) +
  geom_col(width = 0.7) +
  facet_wrap(.~Fragment_mz)





###############################################################
#organize metadata for volume filtered and experiment time
library(lubridate)



m.dat <- read_csv("Kin_metadata.csv") %>%
  select(Exp, Treatment, Samp_ID, Vol_Filt_L, Local_Date, Spike_Time, 
         Filt_time_start, Filt_time_end, time_elapsed, incubation_time) %>%
  unite("Spike_time", c(Local_Date, Spike_Time), sep = " ", remove = FALSE) %>%
  unite("Filt_time_start", c(Local_Date, Filt_time_start), sep = " ", remove = FALSE) %>%
  unite("Filt_time_end", c(Local_Date, Filt_time_end), sep = " ", remove = FALSE) %>%
  mutate("Spike_time" = dmy_hms(Spike_time),
         "Filt_time_start" =  dmy_hms(Filt_time_start), 
         "Filt_time_end" = dmy_hms(Filt_time_end)) %>%
  select(Exp, Treatment, Samp_ID, Vol_Filt_L, Spike_time, Filt_time_start, Filt_time_end) %>%
  mutate(exp_time = as.numeric(Filt_time_start - Spike_time),
         filt_time = as.numeric(Filt_time_end - Filt_time_start),
         incubation_time1 = exp_time + (filt_time/2),
         incubation_time2 = exp_time + filt_time)

m.dat.2 <- m.dat %>%
  rename("SampID" =  Samp_ID) %>%
  separate(SampID, into = c("cruise", "exp", "treatment", "rep"), remove = FALSE) 


m.dat.3 <- m.dat.2 %>%
  select(exp, treatment, rep, Vol_Filt_L, exp_time, filt_time, incubation_time1, incubation_time2) 

###organize uptake kinetics data
#dat.1.smp <- dat.1.1 %>%
#  filter(str_detect(SampID, "Smp")) %>%
#  mutate(SampID = str_replace(SampID, "OnM", "0nM")) %>%
#  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
#  filter(Compound == "D3-Homarine")

dat.join <- dat.smp%>%
  select(exp, treatment, rep, Area, Fragment_mz) %>%
  left_join(., m.dat.3) %>%
  left_join(., lm.out.slopes) %>%
  mutate(nM_in_vial = Area/slope,
         nM_in_smp = (nM_in_vial*400e-6)/Vol_Filt_L,
         nM_per_hour = (nM_in_smp/incubation_time1)*60) 
  
dat.join.adj <- dat.smp%>%
  select(exp, treatment, rep, Area, Fragment_mz) %>%
  left_join(., m.dat.3) %>%
  left_join(., lm.out.slopes) %>%
  mutate(nM_in_vial = Area/slope,
         nM_in_smp = (nM_in_vial*400e-6)/Vol_Filt_L,
         nM_per_hour = (nM_in_smp/incubation_time1)*60) %>% 
  mutate(treatment_nM = as.numeric(str_remove(treatment, "nM")),
         nM_in_vial_adj = nM_in_vial - 0.158*treatment_nM,
         nM_in_smp_adj = (nM_in_vial_adj*400e-6)/Vol_Filt_L,
         nM_per_hour_adj = (nM_in_smp_adj/incubation_time1)*60)



########UCH1
UCH1.1 <- dat.join %>%
  filter(exp == "UCH1")  %>%
  filter(Fragment_mz == 78.0399) %>%
  filter(!treatment == "Hom")

UCH1.1.means <- UCH1.1 %>%
  group_by(treatment, Fragment_mz) %>%
  mutate(mean_rate = mean(nM_per_hour)) %>%
  select(exp, treatment, mean_rate, Fragment_mz) %>%
  filter(Fragment_mz == 78.0399) %>%
  unique()

ggplot(UCH1.1, aes(x = treatment, y = nM_per_hour)) +
  geom_col(data = UCH1.1.means, aes(x = treatment, y = mean_rate), width = 0.5, alpha = 0.8, fill = "darkblue") + 
  geom_point(color = "darkgrey") +
  ylab("Uptake rate (nM/hr)") +
  theme_classic()
 # facet_wrap(.~Fragment_mz)

##vague statistics 



###Manual inspection of uptake curves to define starting predictions for nonlinear model fitting

##
dat.test <- dat.join.adj %>%  
  filter(exp %in% c("UKH1", "UKH2", "UKH3", "UKH4", "UKH5")) %>%
  group_by(exp) %>%
  mutate(max.rate = max(nM_per_hour_adj, na.rm = TRUE))


###Examine raw data to define starting predictions for Vmax (a) and Ks (b) 
# the red line is the max uptake rate and likely represents a good prediction for Vmax,
ggplot(dat.test, aes(x = treatment_nM, y = nM_per_hour_adj))  +
  geom_point() +
  facet_wrap(.~exp, scales = "free") +
  geom_hline(aes(yintercept = max.rate), color = "red")


##Inspect plots and manually enter predicted values into dataframe, 
# enter values in order of experiment (exp) for a_pred (Vmax prediction) and b_pred (Ks prediction)

model.starting.estiamtes <- data.frame(
  exp = c("UKH1", "UKH2", "UKH3", "UKH4", "UKH5"),
  a_pred = c(125, 200, 300, 500, 125),
  b_pred = c(0.1, 0.1, 0.25, 0.3, 0.6)
)






####UKH4

###Pull out just data from that experiment for the desired fragment, calculate 
# the moles of compound in the spike and then time/(inverse of the fraction taken up) 
# 

fit_nls <- function(dat.input.rate, dat.input.pred, exp_name) {
  
#pull out uptake rate data and starting predictions for the experiment defined in "exp_name"
  dat.exp <- dat.input.rate %>%
    filter(exp == exp_name)
  pred.exp <- dat.input.pred %>%
    filter(exp == exp_name)
  
#define model inputs
  x <- dat.exp$treatment_nM   #treatment concentrations
  y <- dat.exp$nM_per_hour_adj   #uptake rates
  a_start <- pred.exp$a_pred    #Vmax prediction
  b_start <- pred.exp$b_pred #Ks prediction

#Fit nonlinear equation to data
  nls.out <- nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
                   control = nls.control(printEval = F))

#Extract model fit parameters from nls.out 
  nls.sum <- summary(nls.out)
  nls.vals <- data.frame(nls.sum$coefficients) %>%
    select(Estimate) %>%
    rownames_to_column(var = "x") %>%
    mutate(parameter = case_when(x == "a" ~ "Ks",
                                 x == "b" ~ "Vmax")) %>%
    mutate(experiment = exp_name) %>%
    select(experiment, parameter, Estimate)
}

dat.practice <- dat.join.adj %>%
  filter(Fragment_mz == 78.0399)

nls.practice <- fit_nls(dat.practice, model.starting.estiamtes, "UKH5")



UKH4.1 <- dat.join.adj %>%
  filter(exp == "UKH4") %>%
  filter(Fragment_mz == 78.0399) %>%
  #mutate(treatment_nM = as.numeric(str_remove(treatment, "nM"))) %>%
  mutate(spike_moles = treatment_nM*Vol_Filt_L,
         inv_fraction_taken_up = (incubation_time1/60)/((nM_in_smp_adj*Vol_Filt_L)/spike_moles))

####Fit NLS model 
## nlm k w/o explicit d[gbt] ------
x <- UKH4.1$treatment_nM
y <-  UKH4.1$nM_per_hour_adj

#from this graph set approximate starting values
a_start<-500 #param a is the Ks + Sn
b_start<-0.4 #b is Vmax
#model
UKH4.nls.1 <-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
              control = nls.control(printEval = T))

UKH4.nls.1
sum.4 <- summary(UKH4.nls.1)
vals.4 <- data.frame(sum.4$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "a" ~ "Ks",
                               x == "b" ~ "Vmax")) %>%
  mutate(experiment = "UKH4") %>%
  select(experiment, parameter, Estimate)


###Adusted: a = 403.8, b = 0.318


nls.fun <- function(x) (0.318*x)/(403.862+x)

###
ggplot(UKH4.1, aes(x = treatment_nM, y = nM_per_hour_adj))  +
  geom_point() +
  geom_function(fun = nls.fun, color = "red")


UKH4.1.no0 <- UKH4.1 %>%
  filter(!treatment_nM == 0) %>%
  mutate(log10_treat_nm = log10(treatment_nM),
         log10_inv_fract = log10(inv_fraction_taken_up))

ggplot(UKH4.1.no0, aes(x = treatment_nM, y = inv_fraction_taken_up)) +
  geom_point() +
  geom_smooth(method = lm)

tt.ukh4 <- lm(inv_fraction_taken_up ~ treatment_nM, data = UKH4.1.no0)
tt.ukh4

tt.ukh4.log10 <- lm(log10_inv_fract ~ log10_treat_nm, data = UKH4.1.no0)
tt.ukh4.log10
lm.4 <- summary(tt.ukh4.log10)
tt.4 <- data.frame(lm.4$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "(Intercept)" ~ "turnover_time",
                               TRUE ~ x)) %>%
  filter(parameter == "turnover_time") %>%
  mutate(Estimate = 10^Estimate) %>%
  mutate(experiment = "UKH4") %>%
  select(experiment, parameter, Estimate)


ggplot(UKH4.1.no0, aes(x = log10_treat_nm, y = log10_inv_fract)) +
  geom_point() +
  geom_smooth(method = lm)






####UKH5
UKH5.1 <- dat.join.adj %>%
  filter(exp == "UKH5") %>%
  filter(Fragment_mz == 78.0399) %>%
  mutate(treatment_nM = as.numeric(str_remove(treatment, "nM"))) %>%
  mutate(spike_moles = treatment_nM*Vol_Filt_L,
         inv_fraction_taken_up = (incubation_time1/60)/((nM_in_smp_adj*Vol_Filt_L)/spike_moles))

####NLS function

## nlm k w/o explicit d[gbt] ------
x <- UKH5.1$treatment_nM
y <-  UKH5.1$nM_per_hour_adj

#from this graph set approximate starting values
a_start<-200 #param a is the Ks + Sn
b_start<-0.76 #b is Vmax
#model
UKH5.nls.1 <-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
                 control = nls.control(printEval = T))

UKH5.nls.1


sum.5 <- summary(UKH5.nls.1)
vals.5 <- data.frame(sum.5$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "a" ~ "Ks",
                               x == "b" ~ "Vmax")) %>%
  mutate(experiment = "UKH5") %>%
  select(experiment, parameter, Estimate)


# unadjusted a = 246, b = 0.6878
# adjusted a = 210, b = 0.607
UKH5.nls.fun.1 <- function(x) (0.607*x)/(x+210)

ggplot(UKH5.1, aes(x = treatment_nM, y = nM_per_hour_adj))  +
  geom_point() +
  geom_function(fun = UKH5.nls.fun.1, color = "red")




#####Define Substrate value:
s <- 50

UKH5.1.test <- UKH5.1 %>%
  mutate(treatment_conc_adj = treatment_nM+s)

x <- UKH5.1.test$treatment_conc_adj
y <-  UKH5.1.test$nM_per_hour


UKH5.nls.2 <- nls(y~(b*(x+s))/(a+x+s), start = list(a=a_start, b=b_start),
                  control = nls.control(printEval = T))

UKH5.nls.2

fitted(UKH5.nls.2)
fit.a.2 <- as.numeric(coefficients(UKH5.nls.2)[1])
fit.b.2 <- as.numeric(coefficients(UKH5.nls.2)[2])

UKH5.nls.fun.2 <- function(x) (fit.b.2*(x+s))/(x+fit.a.2+s)

ggplot(UKH5.1, aes(x = treatment_nM, y = nM_per_hour))  +
  geom_point() +
  geom_function(fun = UKH5.nls.fun.2, color = "red")

#

ggplot(UKH5.1, aes(x = treatment_nM, y = nM_per_hour))  +
  geom_point() +
  geom_function(fun = UKH5.nls.fun, color = "red")

UKH5.1.no0 <- UKH5.1 %>%
  filter(!treatment_nM == 0) %>%
  mutate(log10_treat_nm = log10(treatment_nM),
         log10_inv_fract = log10(inv_fraction_taken_up))


tt.ukh5 <- lm(log10_inv_fract ~ log10_treat_nm, data = UKH5.1.no0)
tt.ukh5
lm.5 <- summary(tt.ukh5)
tt.5 <- data.frame(lm.5$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "(Intercept)" ~ "turnover_time",
                               TRUE ~ x)) %>%
  filter(parameter == "turnover_time") %>%
  mutate(Estimate = 10^Estimate) %>%
  mutate(experiment = "UKH5") %>%
  select(experiment, parameter, Estimate)

#TT
#unadjusted: TT = 112 hr (4.6)
# adjusted: TT = 111 hr (~4.6 days)



ggplot(UKH5.1.no0, aes(x = treatment_nM, y = inv_fraction_taken_up)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = lm)
  
ggplot(UKH5.1.no0, aes(x = log10_treat_nm, y = log10_inv_fract)) +
  geom_point() +
  geom_smooth(method = lm)




####UKH1
UKH1.1 <- dat.join.adj %>%
  filter(exp == "UKH1") %>%
  filter(Fragment_mz == 78.0399) %>%
  mutate(treatment_nM = as.numeric(str_remove(treatment, "nM"))) %>%
  mutate(spike_moles = treatment_nM*Vol_Filt_L,
         inv_fraction_taken_up = (incubation_time1/60)/((nM_in_smp_adj*Vol_Filt_L)/spike_moles))

####Fit NLS model 
## nlm k w/o explicit d[gbt] ------
x <- UKH1.1$treatment_nM
y <-  UKH1.1$nM_per_hour_adj

#from this graph set approximate starting values
a_start<-200 #param a is the Ks + Sn
b_start<-0.1 #b is Vmax
#model
UKH1.nls.1 <-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
                 control = nls.control(printEval = T))

sum.1 <- summary(UKH1.nls.1)
vals.1 <- data.frame(sum.1$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "a" ~ "Ks",
                               x == "b" ~ "Vmax")) %>%
  mutate(experiment = "UKH1") %>%
  select(experiment, parameter, Estimate)


##unadjsuted: a = 133.6, b = 0.1493 
##adjusted: a = 60.57, b = 0.09528


nls.fun <- function(x) (0.09528*x)/(60.57+x)

###
ggplot(UKH1.1, aes(x = treatment_nM, y = nM_per_hour_adj))  +
  geom_point(size = 2) +
  geom_function(fun = nls.fun, color = "red", size = 1) +
  geom_hline(yintercept = 0.09528, color = "blue") +
  geom_segment(x = 60.57, y = 0, xend = 60.57, yend = 0.09628/2, linetype = "dashed", color = "blue") +
  geom_segment(x = 0, y = 0.09628/2, xend = 60.57, yend = 0.09628/2, linetype = "dashed", color = "blue") +
  xlab("Spike Concentration (nM)") +
  ylab("Uptake Rate (nM/hr)") +
  theme_bw()


UKH1.1.no0 <- UKH1.1 %>%
  filter(!treatment_nM == 0) %>%
  mutate(log10_treat_nm = log10(treatment_nM),
         log10_inv_fract = log10(inv_fraction_taken_up))

ggplot(UKH1.1.no0, aes(x = treatment_nM, y = inv_fraction_taken_up)) +
  geom_point() +
  geom_smooth(method = lm)

tt.ukh1 <- lm(inv_fraction_taken_up ~ treatment_nM, data = UKH1.1.no0)
tt.ukh1

tt.ukh1.log10 <- lm(log10_inv_fract ~ log10_treat_nm, data = UKH1.1.no0)
tt.ukh1.log10

lm.1 <- summary(tt.ukh1.log10)
tt.1 <- data.frame(lm.1$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "(Intercept)" ~ "turnover_time",
                               TRUE ~ x)) %>%
  filter(parameter == "turnover_time") %>%
  mutate(Estimate = 10^Estimate) %>%
  mutate(experiment = "UKH1") %>%
  select(experiment, parameter, Estimate)


ggplot(UKH1.1.no0, aes(x = log10_treat_nm, y = log10_inv_fract)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  theme_bw() +
  xlab("log10(Spike Concentration) (nM)") +
  ylab("log10(fraction of metabolite taken up per hour) (hr^-1)")

#unadjusted: TT = 158 hr = 6.5 days
#adjusted: TT = 145 hr = 6 days 




####UKH2
UKH2.1 <- dat.join.adj %>%
  filter(exp == "UKH2") %>%
  filter(Fragment_mz == 78.0399) %>%
  mutate(treatment_nM = as.numeric(str_remove(treatment, "nM"))) %>%
  mutate(spike_moles = treatment_nM*Vol_Filt_L,
         inv_fraction_taken_up = (incubation_time1/60)/((nM_in_smp_adj*Vol_Filt_L)/spike_moles))

####Fit NLS model 
## nlm k w/o explicit d[gbt] ------
x <- UKH2.1$treatment_nM
y <-  UKH2.1$nM_per_hour_adj

#from this graph set approximate starting values
a_start<-200 #param a is the Ks + Sn
b_start<-0.2 #b is Vmax
#model
UKH2.nls.1 <-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
                 control = nls.control(printEval = T))

UKH2.nls.1

nls.fun <- function(x) (0.077*x)/(83.75+x)


sum.2 <- summary(UKH2.nls.1)
vals.2 <- data.frame(sum.2$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "a" ~ "Ks",
                               x == "b" ~ "Vmax")) %>%
  mutate(experiment = "UKH2") %>%
  select(experiment, parameter, Estimate)

####unadjusted: a = 329, b = 0.17
####adjusted: a = 83.75, b = 0.077


###
ggplot(UKH2.1, aes(x = treatment_nM, y = nM_per_hour_adj))  +
  geom_point() +
  geom_function(fun = nls.fun, color = "red")


UKH2.1.no0 <- UKH2.1 %>%
  filter(!treatment_nM == 0) %>%
  mutate(log10_treat_nm = log10(treatment_nM),
         log10_inv_fract = log10(inv_fraction_taken_up))

ggplot(UKH2.1.no0, aes(x = treatment_nM, y = inv_fraction_taken_up)) +
  geom_point() +
  geom_smooth(method = lm)

tt.ukh2 <- lm(inv_fraction_taken_up ~ treatment_nM, data = UKH2.1.no0)
tt.ukh2

tt.ukh2.log10 <- lm(log10_inv_fract ~ log10_treat_nm, data = UKH2.1.no0)
tt.ukh2.log10

lm.2 <- summary(tt.ukh2.log10)
tt.2 <- data.frame(lm.2$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "(Intercept)" ~ "turnover_time",
                               TRUE ~ x)) %>%
  filter(parameter == "turnover_time") %>%
  mutate(Estimate = 10^Estimate) %>%
  mutate(experiment = "UKH2") %>%
  select(experiment, parameter, Estimate)

##unadjusted: TT = 363 hr, 15.1 days
##adjusted: TT = 331 hr, 13.7 days 


ggplot(UKH2.1.no0, aes(x = log10_treat_nm, y = log10_inv_fract)) +
  geom_point() +
  geom_smooth(method = lm)

#TT = 367 hr = 15 days




####UKH3
UKH3.1 <- dat.join.adj %>%
  filter(exp == "UKH3") %>%
  filter(Fragment_mz == 78.0399) %>%
  mutate(treatment_nM = as.numeric(str_remove(treatment, "nM"))) %>%
  mutate(spike_moles = treatment_nM*Vol_Filt_L,
         inv_fraction_taken_up = (incubation_time1/60)/((nM_in_smp_adj*Vol_Filt_L)/spike_moles))

####Fit NLS model 
## nlm k w/o explicit d[gbt] ------
x <- UKH3.1$treatment_nM
y <-  UKH3.1$nM_per_hour_adj

#from this graph set approximate starting values
a_start<-400 #param a is the Ks + Sn
b_start<-0.2 #b is Vmax
#model
UKH3.nls.1 <-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
                 control = nls.control(printEval = T))

UKH3.nls.1


sum.3 <- summary(UKH3.nls.1)
vals.3 <- data.frame(sum.3$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "a" ~ "Ks",
                               x == "b" ~ "Vmax")) %>%
  mutate(experiment = "UKH3") %>%
  select(experiment, parameter, Estimate)

###unadjusted: a = 531.5, b = 0.3296
###adjusted: a = 309.718, b = 0.2146



nls.fun <- function(x) (0.215*x)/(309.7+x)

###
ggplot(UKH3.1, aes(x = treatment_nM, y = nM_per_hour_adj))  +
  geom_point() +
  geom_function(fun = nls.fun, color = "red")


UKH3.1.no0 <- UKH3.1 %>%
  filter(!treatment_nM == 0) %>%
  mutate(log10_treat_nm = log10(treatment_nM),
         log10_inv_fract = log10(inv_fraction_taken_up))

ggplot(UKH3.1.no0, aes(x = treatment_nM, y = inv_fraction_taken_up)) +
  geom_point() +
  geom_smooth(method = lm)

tt.ukh3 <- lm(inv_fraction_taken_up ~ treatment_nM, data = UKH3.1.no0)
tt.ukh3

tt.ukh3.log10 <- lm(log10_inv_fract ~ log10_treat_nm, data = UKH3.1.no0)
tt.ukh3.log10

lm.3 <- summary(tt.ukh3.log10)
tt.3 <- data.frame(lm.3$coefficients) %>%
  select(Estimate) %>%
  rownames_to_column(var = "x") %>%
  mutate(parameter = case_when(x == "(Intercept)" ~ "turnover_time",
                               TRUE ~ x)) %>%
  filter(parameter == "turnover_time") %>%
  mutate(Estimate = 10^Estimate) %>%
  mutate(experiment = "UKH3") %>%
  select(experiment, parameter, Estimate)


ggplot(UKH3.1.no0, aes(x = log10_treat_nm, y = log10_inv_fract)) +
  geom_point() +
  geom_smooth(method = lm)

#unadjusted: TT = 368 hr = 15.3 days
#adjusted: TT = 342 hr = 14.2 days 

###Summarize G4 Homarine Kinetics Experiments

G4.tt.sum <- rbind(tt.1, tt.2, tt.3, tt.4, tt.5)
ggplot(G4.tt.sum, aes(x = experiment, y = Estimate, fill = experiment)) +
  geom_col(width = 0.7)


G4.kin.sum <- rbind(vals.1, vals.2, vals.3, vals.4, vals.5)
ggplot(G4.kin.sum, aes(x = experiment, y = Estimate, fill = experiment)) +
  geom_col(width = 0.7) +
  facet_wrap(.~parameter, scales = "free")


G4.kin.sum.all <- rbind(G4.tt.sum, G4.kin.sum)

write_csv(G4.kin.sum.all, file = "G4_kin_sum.csv")




#####Blanks
dat.blk <- dat.join %>%
  filter(rep == "Blk") %>%
  mutate(treatment_nM = as.numeric(str_remove(treatment, "nM"))) 

ggplot(dat.blk, aes(x = treatment_nM, y = nM_in_vial)) +
  geom_smooth(method = "lm") +
  geom_point(size = 3, alpha = 0.8) 

ggplot(dat.blk, aes(x = treatment_nM, y = nM_in_vial)) +
  geom_smooth(method = "lm") +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10()

lm.blk <- lm(nM_in_vial ~ treatment_nM, dat = dat.blk)
lm.blk


#############Calculate Homarine Concentrations in environmental samples
###Wrangle Sample Data
dat.1.hom <- dat.1.1 %>%
  filter(str_detect(SampID, "Smp")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) %>%
  filter(Compound == "Homarine")

dat.3.hom <- dat.3 %>%
  filter(!str_detect(SampID, "inM")) %>%
  filter(!str_detect(SampID, "Std")) %>%
  filter(Compound == "Homarine") %>%
  mutate(SampID = str_replace(SampID, "808_U", "808_Smp_TN397_U")) %>%
  separate(col = SampID, into = c("date", "type", "cruise", "exp", "treatment", "rep"), remove = FALSE) 

lm.out.slopes.2 <- lm.out.slopes %>%
  filter(Fragment_mz == 97.0954)

dat.hom.env <- rbind(dat.1.hom, dat.3.hom) %>%
  filter(Fragment_mz == 94.1067) %>%
  select(exp, treatment, rep, Area) %>%
  left_join(., m.dat.3) %>%
  left_join(., lm.out.slopes.2) %>%
  mutate(nM_in_vial = Area/slope,
         nM_in_smp = (nM_in_vial*400e-6)/Vol_Filt_L)

ggplot(dat.hom.env, aes(x = exp, y = nM_in_smp)) +
  geom_boxplot()+
  geom_point() 
  
G4.hom.env <- dat.hom.env %>%
  group_by(exp) %>%
  summarize(env_nM = mean(nM_in_smp, na.rm = TRUE)) %>%
  filter(exp %in% c("UKH1", "UKH2", "UKH3", "UKH4", "UKH5"))

write_csv(G4.hom.env, file = "G4_homkin_evn_concs.csv")













#############additional code:
#dat.2 <- read_csv("G4_UKH5_TQS_results.csv")

d3hom.1 <- dat.1 %>%
  filter(`Precursor Ion Name` == "D3-Homarine") %>%
  filter()

ggplot(d3hom.1, aes(x = `Replicate Name`, y = Area)) +
  geom_col() +
  facet_wrap(.~`Product Mz`) +
  theme(axis.text.x = element_text(angle = 90))

d3hom.2 <- dat.2 %>%
  filter(Compound == "D3-Homarine",
         Fragment_mz == 78.0399)

ggplot(d3hom.2, aes(x = SampID, y = Area)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90))