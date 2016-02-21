##############################################  NEST Beach closure project  #######################################
###########################################  Figures and Analyses production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 02/20/2016
#Version: 2
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
library(parallel)               # access to parallelization functions

###### Functions used in this script sourced from other files

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_02212016.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
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

in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50 , param 1
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS, param 
#in_dir_rainfall <- "/home/bparmentier/Google Drive/NEST_Data/"
  
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 #param 2
CRS_reg <- CRS_WGS84 # PARAM 3

file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_02172016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

rainfall_dir <- "/home/bparmentier/Google Drive/NEST_Data" #PARAM 10
station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST_Data/", "WQ_TECS_Q.txt") #PARAM 11

start_date <- "2012-01-01" #PARAM 12
end_date <- "2012-12-31" #PARAM 13
var_name <- "COL_SCORE" #PARAM 14
var_ID <- "LOCATION_ID" #PARAM 15
year_processed <- "2012" #PARAM 16
threshold_val <- 2*25.4 #PARAM 17, in inches or mm
convert_to_inches <- FALSE #PARAM 18
units_val <- "mm"

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

list_dir_rainfall <- list.dirs(path=rainfall_dir,full.names=T)
#remove non relevant directories

#### Part 1: read in and combine the information ####

data <- read.table(station_data_fname,sep=",",header=T,fill=T,stringsAsFactors = F) #bacteria measurements
data <- read.table(station_data_fname,sep=",",header=T,stringsAsFactors = F) #bacteria measurements

in_dir_rst <- list_dir_rainfall[2]

#
#> data <- read.table(station_data_fname,sep=",",header=T) #this is T-mode using cor matrix
#Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : 
#                line 47 did not have 35 elements

debug(combine_stations_data_raster_ts_fun)
df_combined <- combine_stations_data_raster_ts_fun(data,convert_to_inches,in_dir_rst,start_date,end_date,out_dir,out_suffix)


###### Part 2: plot information ####

plot(r_rainfall,y=1)

res_pix <- 480
col_mfrow <- 1
row_mfrow <- 1

png(filename=paste("Figure1_","histogram","_station_coliform_measurements_frequency_",year_processed,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(table(data_subset$LOCATION_ID),type="h", main="Number of measurements",
     ylab="Frequency of coliform measurements",xlab="Station ID")

dev.off()

#hist(table(data_subset$COL_SCORE))
#Could output in a textfile
range(table(dat_stat$LOCATION_ID)) #1 to 60
mean(table(dat_stat$LOCATION_ID)) #8
median(table(data_subset$LOCATION_ID)) #8
hist(table(data_subset$LOCATION_ID), main="Frequency of Number of measurements by stations")

ref_pol <- create_polygon_from_extent(dat_stat,outDir=out_dir,outSuffix=out_suffix)
  
res_pix <- 480
col_mfrow <- 1.5
row_mfrow <- 1

png(filename=paste("Figure2a_","rainfall_map_",dates_l[1],"_and_stations_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_rainfall,y=1)
plot(dat_stat,add=T)
text(dat_stat,dat_stat$tID,cex=1.4)
plot(ref_pol,border="red",add=T)
legend("topright",legend=c("stations"), 
       cex=1.2, col="black",pch =3,bty="n")

dev.off()

col_mfrow <- 2
row_mfrow <- 1

png(filename=paste("Figure2b_","rainfall_map_",dates_l[1],"_and_stations_zoom_window_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_rainfall,y=1,ext=extent(dat_stat))
plot(dat_stat,add=T,pch=3)
text(dat_stat,dat_stat$tID,cex=1.4)
legend("topright",legend=c("stations"), 
       cex=1.2, col="black",pch =3,bty="n")

dev.off()

#
#make a function later on?

#ADD selection function by ID
#list_selected_pix: by LOCATION ID
#


#colnames(data_df) <- list_selected_ID
#data_dz <- zoo(data_dz,idx)

plot(log(data_df$rainfall),log(data_df$COL_SCORE))
plot(data_df$rainfall,data_df$COL_SCORE)

## IDENTIFY 2 inches events?


############################ END OF SCRIPT #######################
