####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 11/09/2015
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - 
#         - 
#TO DO:
#
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
script_path <- "/home/bparmentier/~/Google Drive/NEST/R_NEST" #path to script #PARAM 2
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

in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS

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
station_data_fname <- file.path(in_dir, "WQ_TECS_Q.xlsx")

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
data <-read.xls(station_data_fname, sheet=1) #this is T-mode using cor matrix

idx <- seq(as.Date('2014-01-01'), as.Date('2014-12-31'), 'day')
date_l <- strptime(idx[1], "%Y%m%d") # interpolation date being processed
dates_l <- format(idx, "%Y%m%d") # interpolation date being processed

## Plot mosaics for Maine for daily predictions in 2014

res_pix <- 480

col_mfrow <- 2
row_mfrow <- 1

#  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
#    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
png(filename=paste("Figure1_","time_series_step_in_raster_mosaics",dates_l[11],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r1)
plot(summary_metrics_v,add=T)
text(summary_metrics_v,summary_metrics_v$tile_id,cex=1.4)

dev.off()

## Get pixel time series at centroids of tiles used in the predictions

df_ts_pixel <- extract(r_stack,summary_metrics_v,df=T,sp=F)

df_ts_pixel <- cbind(summary_metrics_v,df_ts_pixel)

#make a function later on?

#inputs
list_selected_pix <- c("tile_4","tile_6","tile_8","tile_11","tile_14","tile_3","tile_5","tile_7","tile_38","tile_12")
list_pix <- vector("list",length=length(list_selected_pix))
#idx <- seq(as.Date('2010-01-01'), as.Date('2010-12-31'), 'day')
#df_ts_pix

#Select one pix to profile/plot

df_ts_pix <- subset(df_ts_pixel,pred_mod=="mod1")

for(i in 1:length(list_selected_pix)){
  
  selected_pix <- list_selected_pix[i]
  data_pixel <- subset(df_ts_pix,tile_id==selected_pix)
  pix <- t(data_pixel[1,24:388])
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

list_selected_pix <- c("tile_4","tile_6","tile_8","tile_11","tile_14","tile_3","tile_5","tile_7","tile_38","tile_12")
df_ts_pix2 <- subset(df_ts_pix,tile_id%in% list_selected_pix)

pix_data <- t(df_ts_pix2[,24:388])

#d_z2 <- zoo(pix_data,idx)
#names(d_z2)<-
  
