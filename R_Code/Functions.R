###
#
#
#
#
#Functions for ProMo Analysis
#
#
library(tidyverse)
library(readr)
#
#
#
#
#####Define Skyline read data read in function
sky_read <- function(file) {
  output <- read_csv(file) %>%
    mutate(Rep =`Replicate Name`,
           Compound = `Precursor Ion Name`,
           RT = as.numeric(`Retention Time`),
           Area = as.numeric(`Area`),
           Background = as.numeric(`Background`),
           Height = as.numeric(`Height`),
           Mass_Error_PPM = as.numeric(`Mass Error PPM`)
    ) %>%
    select(Rep, Compound, RT, Area, Background, Height, Mass_Error_PPM)
  print(output)
}

####Define Skyline transition list read in fucntion
tl_read <- function(file) {
  output <- read_csv(file, col_names = FALSE) %>%
    mutate(Mass = as.numeric(X2),
           Compound = X4) %>%
    select(Compound, Mass)
}

###Define function to read in MS-DIAL Data
#mode should be the analytical fraction ("HILICPos", "HILICNeg", or "RP")

MSDIAL_read <- function(file1, Mode) {
  read_delim(file1,
             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = Mode) %>% 
  mutate("ID" = as.character(`Alignment ID`)) %>%
  rename("RT" = `Average Rt(min)`,
         'mz' = `Average Mz`,
         "Name" = `Metabolite name`, 
         "Adduct" = `Adduct type`,
         "Note" = `Post curation result`,
         "Fill" = `Fill %`,
         "MS2" = `MS/MS assigned`,
         "Ref_RT" = `Reference RT`,
         "Ref_mz" = `Reference m/z`,
         "RT_matched" = `RT matched`,
         "mz_matched" = `m/z matched`,
         "MS2_matched" = `MS/MS matched`,
         "SN_ave" = `S/N average`)
}




####Adduct, isotope, and pos/neg match finder function
###adduct finding function
find_adducts <- function(ID_num, dataset, adduct_list) {
  
  dataset.2 <- dataset %>%
    mutate(merge.key = "x")  
  
  MF.limits <- dataset %>%
    filter(ID == ID_num) %>%
    mutate(RT_low = RT - RT_ad_tol,
           RT_high = RT + RT_ad_tol) 
  ad.limits <- cbind(MF.limits, adduct_list) %>%
    mutate(adduct.mass = ((mz) - (1.007276*polarity))/abs(Charge) + Mass.change,
           adduct.mass.high = adduct.mass+(adduct.error.ppm/10^6)*adduct.mass,
           adduct.mass.low = adduct.mass-(adduct.error.ppm/10^6)*adduct.mass) %>%
    mutate(merge.key = "x") %>%
    select(merge.key, Ion, RT_low, RT_high, adduct.mass.low, adduct.mass.high, ad.polarity) %>%
    rename("polarity" = ad.polarity)
  
  ad.detect <- full_join(dataset.2, ad.limits) %>%
    select(!merge.key) %>%
    rowwise() %>%
    filter(RT >= RT_low & RT <= RT_high) %>%
    filter(mz >= adduct.mass.low & mz <= adduct.mass.high) %>%
    rename("adduct.ID" = ID) %>%
    mutate(ID = ID_num) 
  
  print(ad.detect)
}

###Function
find_sirius_match <- function(MF.file, SIRIUS.file, RT_tol, search.error.ppm) {
  
  MF.limits <- MF.file %>%
    mutate(RT_low = RT.seconds - RT_tol,
           RT_high = RT.seconds + RT_tol) %>%
    mutate(search.mass = mz,
           search.mass.high = search.mass+(search.error.ppm/10^6)*search.mass,
           search.mass.low = search.mass-(search.error.ppm/10^6)*search.mass) %>%
    mutate(match.key = "x")
  
  sirius.MFs <- SIRIUS.file %>%
    select(id, ionMass, retentionTimeInSeconds) %>%
    mutate(match.key = "x")
  
  MF.match <- full_join(MF.limits, sirius.MFs) %>%
    rowwise() %>%
    filter(retentionTimeInSeconds >= RT_low & retentionTimeInSeconds <= RT_high) %>%
    filter(ionMass >= search.mass.low & ionMass <= search.mass.high) 
  
  MF.annotate.full <- MF.match %>%
    select(MF, Name, mz, RT, RT.seconds, id) %>%
    left_join(., SIRIUS.file)
  
  print(MF.annotate.full)
  
}







