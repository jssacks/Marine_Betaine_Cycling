


#Function to read-in, clean-up, and rename targeted metabolomics data collected using PRM from TQS instrument 
TQS_csv_read <- function(file_name) {
  x <- read_csv(file_name) %>%
    rename(SampID = "Replicate Name",
           Compound = "Precursor Ion Name",
           Fragment_mz = "Product Mz") %>%
    select(SampID, Compound, Area, Fragment_mz)
  x
}




####function to fit nls equation to uptake kinetic curves
fit_nls <- function(dat.input.rate, dat.input.pred, exp_name) {
  
  #pull out uptake rate data and starting predictions for the experiment defined in "exp_name"
  dat.exp <- dat.input.rate %>%
    filter(exp == exp_name)
  pred.exp <- dat.input.pred %>%
    filter(exp == exp_name)
  
  #define model inputs
  x <- dat.exp$treatment_nM   #treatment concentrations
  y <- dat.exp$nM_per_hour   #uptake rates
  s <- dat.exp$diss_conc_nM %>%
    unique()
  a_start <- pred.exp$a_pred    #Vmax prediction
  b_start <- pred.exp$b_pred #Ks prediction
  
  #Fit nonlinear equation to data
  nls.out <- nls(y~(b*(x+s))/(a+x+s),start=list(a=a_start,b=b_start),
                 control = nls.control(maxiter = 500, warnOnly = T))
  
  #Extract model fit parameters from nls.out 
  nls.sum <- summary(nls.out)
  nls.vals <- data.frame(nls.sum$coefficients) %>%
    select(Estimate) %>%
    rownames_to_column(var = "x") %>%
    mutate(parameter = case_when(x == "a" ~ "Ks",
                                 x == "b" ~ "Vmax")) %>%
    mutate(experiment = exp_name,
           s_val = s) %>%
    mutate(sigma = nls.sum$sigma,
           cor = cor(y, predict(nls.out)),
           iterations = nls.out$convInfo$finIter) %>%
    select(experiment, parameter, s_val, Estimate, sigma, cor, iterations)
}


