


#load packages:
packages.required <- c("rerddap","plotdap","rerddapXtracto", "lubridate", "maps",
                       "mapdata", "mapproj", "ncdf4", "data.table", "tidyverse",
                       "sp", "ggplot2", "gridExtra", "cmocean", "ggnewscale")

lapply(packages.required, require, character.only = TRUE)



##Define inputs:
#sample.loc.file <- 
#exp.loc.file <- 



#Define time, lat, and lon range:

#Time:
tcoord <- c("2021-11-19", "2021-12-15")

#Lat and Long
xcoord <- c(-175, -110)
ycoord <- c(-10, 50)
ttext<-paste(paste(abs(xcoord), collapse="-"),"W, ", paste(ycoord, collapse="-"),"N")



###Get MODIS data:
dataInfo <- rerddap::info('erdMH1chlamday')
dataInfo

# Extract the parameter name from the metadata in dataInfo
parameter <- dataInfo$variable$variable_name

#Extract the start and end times of the dataset from the metadata in dataInfo
global <- dataInfo$alldata$NC_GLOBAL

# Run rxtracto_3D
chlMODIS<-rxtracto_3D(dataInfo,parameter=parameter,
                      tcoord=tcoord,
                      xcoord=xcoord,ycoord=ycoord)

chlMODIS$avgmap <- apply(chlMODIS$chlorophyll,c(1,2),function(x) mean(x,na.rm=TRUE))

#extract chlorophyll data and pair with lat and long

#chlorophyll 
chl.wide <- as.data.table(chlMODIS$avgmap) 

chl.long <- tibble(melt(chl.wide)) %>%
  select(value) %>%
  rename("Chl" = value)

#lat + lon
lat <- as.data.table(tibble(lat = chlMODIS$latitude)) 
lon <- as.data.table(tibble(lon = chlMODIS$longitude))

chlmap <- cross_join(lat, lon) %>%
  cbind(., chl.long) 


###########XXXXXXX













