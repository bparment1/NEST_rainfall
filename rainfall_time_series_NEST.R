####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 11/21/2015
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

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_11192015.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 2
#script_path <- "/home/parmentier/Data/rainfall/NEST"
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

CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"NEST_prism_11192015" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
num_cores <- 11 #PARAM 14

rainfall_dir <- "prism_rain"
station_data_fname <- file.path(in_dir, "WQ_TECS_Q.txt")

#start_date <- "2014-01-01"
#end_date <- "2014-12-31"
start_date <- "2010-01-01"
end_date <- "2010-12-31"
var_name <- "COL_SCORE"
var_ID <- "LOCATION_ID"

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

dates_TRIP_START <- gsub(" 0:00:00","",data$TRIP_START_DATE)
data$TRIP_START_DATE_f <- as.Date(strptime(dates_TRIP_START,"%m/%d/%Y"))
data$TRIP_START_DATE_month <- strftime(data$TRIP_START_DATE_f , "%m") # current month of the date being processed
data$TRIP_START_DATE_year <- strftime(data$TRIP_START_DATE_f , "%Y")
data$TRIP_START_DATE_day <- strftime(data$TRIP_START_DATE_f , "%d")

idx <- seq(as.Date(start_date), as.Date(end_date), 'day')
#date_l <- strptime(idx[1], "%Y%m%d") # 
dates_l <- format(idx, "%Y%m%d") #  date being processed

data_subset <- data[data$TRIP_START_DATE_f >= as.Date(start_date) & data$TRIP_START_DATE_f <= as.Date(end_date), ]
data_subset$LOCATION_ID <- as.character(data_subset$LOCATION_ID)

dat_stat <- subset(data_subset, !is.na(LONGITUDE_DECIMAL) & !is.na(LATITUDE_DECIMAL))
coords <- dat_stat[,c('LONGITUDE_DECIMAL','LATITUDE_DECIMAL')]
coordinates(dat_stat) <- coords
proj4string(dat_stat) <- projection(r_rainfall) #this is the NAD83 latitude-longitude

## Remove duplicates rows from stations to identify uniques sations
dat_stat <- remove.duplicates(dat_stat)
#dat_stat$LOCATION_ID <- as.character(dat_stat$LOCATION_ID)
nrow(dat_stat)==length(unique(dat_stat$LOCATION_ID)) #Checking that we have a unique identifier for each station

r_rainfall <- setZ(r_rainfall, idx) #for now, this can also be made into a spacetime object

#x <- zApply(r_rainfall, by=as.yearqtr, fun=mean, name="quarters") #aggregate times series by quarter
#x <- zApply(r_rainfall, by=as.yearmon, fun=mean, name="month") #aggregate time series by month
#x <- zApplyr_rainfall, by="month",fun=mean,name="month") #overall montlhy mean mean

x <- zApply(r_rainfall, by="day",fun=mean,name="overall mean") #overall mean
raster_name <- paste("day","_","overall_mean",file_format,sep="")
writeRaster(x, file=raster_name,overwrite=T)
plot_to_file(raster_name) #quick plot of raster to disk

#x <- zApply(r_rainfall, by=c(1,24),fun=mean,name="overall mean") #overall mean
r_date <- getZ(r_rainfall)
#x <- apply.daily(ndvi_ts,FUN=mean) does not work
#plot(x)

## Plot mosaics for Maine for daily predictions in 2014
## Get pixel time series at centroids of tiles used in the predictions

df_ts_pixel <- extract(r_rainfall,dat_stat,df=T,sp=T)
#df_ts_pixel <- cbind(summary_metrics_v,df_ts_pixel)
r_ts_name <- names(r_rainfall)
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

#ADD selection function by ID
#list_selected_pix: by LOCATION ID
#
plot(table(data_subset$LOCATION_ID),type="h")
freq_station <- sort(table(data_subset$LOCATION_ID),decreasing=T) # select top 2 stations in term of availability
list_selected_ID <- names(freq_station)[1:2]

#selected_pix <- list_selected_pix[i]
#data_pixel <- df_ts_pix[selected_pix,]
##id_name <- data_pixel[[var_ID]]
#id_selected <- data_subset[[var_ID]]==list_selected_ID
#test <- data_subset[data_subset$LOCATION_ID %in% list_selected_ID,]
#df_ts_pix <- subset(df_ts_pixel,)

##This will be a function later on...
df_ts_pix <- df_ts_pixel
#list_selected_pix <- 11:14
list_pix <- vector("list",length=length(list_selected_ID))
#threshold <- 

for(i in 1:length(list_selected_ID)){
  
  id_name <- list_selected_ID[i]

  #id_selected <- data_subset[[var_ID]]==id_name
  id_selected <- df_ts_pix[[var_ID]]==id_name
  #df_ts_pix[id_selected,]
  data_pixel <- df_ts_pix[id_selected,]
  #id_name <- data_pixel[[var_ID]]
  var_pix <- subset(as.data.frame(data_subset[id_selected,c(var_name,"TRIP_START_DATE_f")])) #,select=var_name)
  var_pix_ts <- t(as.data.frame(subset(data_pixel,select=var_name)))
  #pix <- t(data_pixel[1,24:388])#can subset to range later
  pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
  
  d_z <- zoo(pix_ts,idx) #make a time series ...
  names(d_z)<- "rainfall"
  d_var <- zoo(var_pix,var_pix$TRIP_START_DATE_f)
  
  d_z2 <- merge(d_z,d_var)
  d_z2$TRIP_START_DATE_f <- NULL
  
  list_pix[[i]] <- pix_ts
  
  res_pix <- 480
  
  col_mfrow <- 2
  row_mfrow <- 1
  
  #  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
  #    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  png(filename=paste("Figure2a_","pixel_profile_",id_name,"_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(d_z2,type="b",main="",pch=10)
  title(paste("Pixel time series",id_name,sep=" "))
  
  dev.off()
  
  ###
  #Figure 2b
  png(filename=paste("Figure2b_","pixel_profile_var_combined_",id_name,"_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(d_z,lty=2,ylab="rainfall",xlab="Time",main="")
  points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
  title(paste("Pixel time series",id_name,sep=" "))
  
  dev.off()
  
  ## Now correlation.
  
  ## Cumulated precip and lag?
  
  ## Threshold?
  
  
}

data_dz <- do.call(cbind,list_pix)
colnames(data_dz) <- list_selected_ID
data_dz <- zoo(data_dz,idx)

## IDENTIFY 2 inches events?


############################ END OF SCRIPT #######################
