####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 11/14/2015
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - The script does not use bacteria data at the current time.  
#          - Add spacetime object functions for later
#TO DO:
# - Select 2 inches rainfall events and correlates with bacteria data
# - Compute accumulated rain over several days using time series functions
# - Make a movie sequence later on using animation package in R
#
#################################################################################################

###Loading R library and packages                                                      

library(raster)                 # loading the raster package
library(gtools)                 # loading R helper programming tools/functions
library(sp)                     # spatial objects in R
library(gplots)                 # plotting functions such as plotCI
library(rgdal)                  # gdal driver for R
library(RColorBrewer)           # color scheme, palettes used for plotting
library(gdata)                  # read different format (including .xlsx)
library(plotrix)                # plot options and functions 
library(rasterVis)              # raster visualization
library(colorRamps)             # contains matlab.like palette
library(zoo)                    # time series objects and methods
library(xts)                    # extension of time series objects
library(BMS)                    # contains hex2bin and bin2hex
library(bitops)                 # bit operations
library(gtools)                 #
library(maptools)               #
library(rgeos)                  # spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(sphet)                  # spatial analyis, regression eg.contains spreg for gmm estimation
library(forecast)               # arima and other time series methods
library(lubridate)              # date and time handling tools

###### Functions used in this script sourced from other files

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_functions.R" #PARAM 1
#script_path <- "/home/bparmentier/~/Google Drive/NEST/R_NEST" #path to script #PARAM 2
script_path <- "/home/parmentier/Data/rainfall/NEST"
source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.

##### Functions used in this script 

create_dir_fun <- function(outDir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}


#####  Parameters and argument set up ###########

#in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50
in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS

#proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"NEST_prism_11052015" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
num_cores <- 11 #PARAM 14

rainfall_dir <- "prism_rain"
station_data_fname <- file.path(in_dir, "WQ_TECS_Q.txt")
#PROJECT_NAME	DMR_TRIP_IDENTIFIER	TRIP_START_DATE	TRIP_COMMENTS	COLLECTOR_INITIALS	EFFORT_START_TIME	LOCATION_ID	LATITUDE_DECIMAL	LONGITUDE_DECIMAL	LOCATION_NAME	LOCATION_DESCRIPTION	LOCATION_TYPE	DIRECTIONS	GROWING_AREA	OPEN_CLOSED_FLAG	WIND_DIRECTION	TIDE_STAGE	CURRENT_CLASSIFICATION_CODE	CATEGORY	DMR_CATCH_IDENTIFIER	SAMPLE_METHOD	TEMP_C	STRATEGY	ADVERSITY	FLOOD	DMR_SAMPLE_IDENTIFIER	LAB	INITIATED_BY	INITIATED_DATE	EXAM_DATE_TIME	SALINITY_PCT	COL_METHOD	COL_SCORE	RAW_COL_SCORE	DELIVERY_TEMP_C

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=file.path(in_dir,rainfall_dir),full.names=T))) #rainfall time series stack

#### Make this a function...that will run automatically the predictions

plot(r_rainfall,y=1)
data <- read.table(station_data_fname,sep=",",header=T) #this is T-mode using cor matrix

dat_stat <- subset(data, !is.na(LONGITUDE_DECIMAL) & !is.na(LATITUDE_DECIMAL))
coords <- dat_stat[,c('LONGITUDE_DECIMAL','LATITUDE_DECIMAL')]
coordinates(dat_stat) <- coords
proj4string(dat_stat) <- CRS_WGS84 #this is the WGS84 projection

## No spatial duplicates
dat_stat <- remove.duplicates(dat_stat)

## No duplicates in attributes
#dat_stat[which(!duplicated(dat_stat$id)), ]

## Combination
#pts[which(!duplicated(as.data.frame(pts))), ]

idx <- seq(as.Date('2014-01-01'), as.Date('2014-12-31'), 'day')
date_l <- strptime(idx[1], "%Y%m%d") # interpolation date being processed
dates_l <- format(idx, "%Y%m%d") # interpolation date being processed

r_rainfall <- setZ(r_rainfall, idx) #for now, this can also be made into a spacetime object

#x <- zApply(r_rainfall, by=as.yearqtr, fun=mean, name="quarters") #aggregate times series by quarter
#names(SISmm) <- month.abb
#x <- zApply(r_rainfall, by=as.yearmon, fun=mean, name="month") #aggregate time series by month
#x <- zApplyr_rainfall, by="month",fun=mean,name="month") #overall montlhy mean mean

x <- zApply(r_rainfall, by="day",fun=mean,name="overall mean") #overall mean
raster_name <- paste("day","_","overall_mean",file_format,sep="")
writeRaster(x, file=raster_name,overwrite=T)
plot_to_file(raster_name) #quick plot of raster to disk

#x <- zApply(r_rainfall, by=c(1,24),fun=mean,name="overall mean") #overall mean
r_date<-getZ(r_rainfall)
#x <- apply.daily(ndvi_ts,FUN=mean) does not work
#plot(x)

## Plot mosaics for Maine for daily predictions in 2014
## Get pixel time series at centroids of tiles used in the predictions

df_ts_pixel <- extract(r_rainfall,dat_stat,df=T,sp=F)
#df_ts_pixel <- cbind(summary_metrics_v,df_ts_pixel)
  
#d_z <- zoo(df_ts_pixel,idx) #make a time series .

res_pix <- 480

col_mfrow <- 2
row_mfrow <- 1

png(filename=paste("Figure1_","time_series_step_in_raster_mosaics",dates_l[1],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_rainfall,y=1)
plot(dat_stat,add=T)
text(dat_stat,dat_stat$tID,cex=1.4)

dev.off()

#make a function later on?

#list_selected_pix

#df_ts_pix <- subset(df_ts_pixel,)
df_ts_pix <- df_ts_pixel
list_selected_pix <- 1:2
list_pix <- vector("list",length=length(list_selected_pix))
for(i in 1:length(list_selected_pix)){
  
  selected_pix <- list_selected_pix[i]
  data_pixel <- subset(df_ts_pix,ID==selected_pix)
  #pix <- t(data_pixel[1,24:388])#can subset to range later
  pix <- t(data_pixel)#can subset to range later
  
  d_z <- zoo(pix,idx) #make a time series ...
  list_pix[[i]] <- pix
  
  res_pix <- 480
  
  col_mfrow <- 2
  row_mfrow <- 1
  
  #  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
  #    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  png(filename=paste("Figure2_","pixel_profile_",selected_pix,"_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  
  plot(d_z,type="b")
  title(paste("Pixel time series",selected_pix,sep=" "))
  
  dev.off()
}

data_dz <- do.call(cbind,list_pix)
colnames(data_dz) <- list_selected_pix
data_dz <- zoo(data_dz,idx)

## IDENTIFY 2 inches events?


############################ END OF SCRIPT #######################
