##############################################  NEST Beach closure project  #######################################
###########################################  Figures and Analyses production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 03/27/2016
#Version: 3
#PROJECT: NEST beach closures            

#
#COMMENTS: - The script does not use bacteria data at the current time.  
#          - Add spacetime object functions for later
#          - Add datatype in name
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
library(hydrostats)

###### Functions used in this script sourced from other files

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_03272016.R" #PARAM 1
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
out_suffix <-"NEST_prism_03272016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

rainfall_dir <- "/home/bparmentier/Google Drive/NEST_Data" #PARAM 10
station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST_Data/", "WQ_TECS_Q.txt") #PARAM 11
station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST/", "MHB_data_2006-2015.csv") #PARAM 11

years_to_process <- 2003:2016
#start_date <- "2012-01-01" #PARAM 12
#end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
#var_name <- "COL_SCORE" #PARAM 14
var_name <- "CONCENTRATION" #PARAM 14, MH data
#var_ID <- "LOCATION_ID" #PARAM 15
var_ID <- NULL #PARAM 15 if null then create a new ID for stations!!
year_processed <- "2012" #PARAM 16
threshold_val <- 2*25.4 #PARAM 17, in inches or mm
convert_to_inches <- FALSE #PARAM 18
units_val <- "mm"
data_type <- "MHB" #for Maine Healthy Beaches
#data_type <- "DMR" #for Maine Department of Marine Resources

coord_names <- c("SITE.LONGITUDE..UTM.","SITE.LATITUDE..UTM.") #MH beach bacteria dataset
#coord_names <- c("LONGITUDE_DECIMAL","LATITUDE_DECIMAL") #cloroforms beach bacteria dataset

SMAZones_fname <- "/home/bparmentier/Google Drive/NEST/NEST_stations_s02/data/SMAZoneDissolve.shp"
SMAZones <- readOGR(dirname(SMAZones_fname),sub(".shp","",basename(SMAZones_fname)))

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
list_dir_rainfall <- grep("prism_ppt*", list_dir_rainfall,value=T)

#remove non relevant directories

#### Part 1: read in and combine the information ####

data <- read.table(station_data_fname,sep=",",header=T,fill=T,stringsAsFactors = F) #bacteria measurements
data$FID <- 1:nrow(data)
#data$ID <- paste(data[[coord_names[1]]],data[[coord_names[2]]],sep="_")

#data <- read.table(station_data_fname,sep=",",header=T,stringsAsFactors = F) #bacteria measurements
#> data <- read.table(station_data_fname,sep=",",header=T) #this is T-mode using cor matrix
#Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : 
#                line 47 did not have 35 elements
### Before combining data get unique station 

#data$LOCATION_ID <- as.character(data$LOCATION_ID)
#length(unique(data$LOCATION_ID)) #There are 2851 stations

##################
##Make this part a function later on...

### Select unique stations and make a SPDF
#data$FID <- 1:nrow(data)
data[[var_name]] <- as.numeric(data[[var_name]]) #make sure we have a numeric field
data[[coord_names[1]]] <- as.numeric(data[[coord_names[1]]])
data[[coord_names[2]]] <- as.numeric(data[[coord_names[2]]])
data <- subset(data, !is.na(data[[coord_names[1]]]) & !is.na(data[[coord_names[2]]]))
data$id_coord <- paste(data[[coord_names[1]]],data[[coord_names[2]]],sep="_")

## Remove duplicates rows from stations to identify uniques sations

id_coord <- unique(data$id_coord)
ID_stat <- 1:length(id_coord)
dat_ID <- data.frame(id_coord,ID_stat)
data <- merge(data,dat_ID,by="id_coord",all=T,suffixes=c("","_y"))
#data_merged <- merge(data,dat_stat,by="FID",all=T,suffixes=c("","_y"))
coords <- data[,coord_names]
coordinates(data) <- coords  
dat_stat <- remove.duplicates(data[,c("ID_stat","id_coord",coord_names[[1]],coord_names[[2]])])
#dat_stat <- remove.duplicates(data[,c("ID_stat","id_coord")])#,coord_names[[1]],coord_names[[2]])])
#coords <- (data[,coord_names])
#coords <- (dat_stat[,coord_names])
#coords[,1] <- as.numeric(coords[,1])
#coords[,2] <- as.numeric(coords[,2])
#coords <- as.matrix(coords)
#coordinates(dat_stat) <- coords  
#dat_stat <- remove.duplicates(dat_stat)
#dat_stat <- remove.duplicates(dat_stat)
#dat_stat <- subset(dat_stat, !is.na(dat_stat$SITE.LONGITUDE..UTM.) & !is.na(dat_stat$SITE.LATITUDE..UTM.))
#coords$LONGITUDE_DECIMAL <- as.numeric(coords$LONGITUDE_DECIMAL)
#coords$LATITUDE_DECIMAL <- as.numeric(coords$LATITUDE_DECIMAL)
#this needs to be changed to be general!!
#dat_stat <- subset(dat_stat, !is.na(dat_stat$SITE.LONGITUDE..UTM.) & !is.na(dat_stat$SITE.LATITUDE..UTM.))
#coords <- (dat_stat[,coord_names])
#coords[,1] <- as.numeric(coords[,1])
#coords[,2] <- as.numeric(coords[,2])
#coords <- as.matrix(coords)
#coordinates(dat_stat) <- coords  

#dat_stat <- subset(data_subset, !is.na(coord_names[1]) & !is.na(coord_names[2]))

#dat_stat <- subset(data, !is.na(data[[coord_names[1]]] & !is.na(data[[coord_names[2]]])))

if(data_type=="MHB"){
  proj_str <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #This will need to be added in the parameters
  proj4string(dat_stat) <- proj_str #this is the NAD83 latitude-longitude
  dat_stat<-spTransform(dat_stat,CRS(CRS_WGS84))     #Project from WGS84 to new coord. system
}
if(data_type=="DMR"){
  proj4string(dat_stat) <- CRS_WGS84 #this is the NAD83 latitude-longitude
}

## Remove duplicates rows from stations to identify uniques sations

if(is.null(var_ID)){
  dat_stat$LOCATION_ID <- dat_stat$ID_stat
  var_ID <- "LOCATION_ID"
}

### Write out basic informaiotn before processing by years

write.table(as.data.frame(dat_stat),file=file.path(out_dir,paste0("dat_stat_location",".txt")),sep=",")
writeOGR(dat_stat,dsn= out_dir,layer= "dat_stat_location", 
          driver="ESRI Shapefile",overwrite_layer="TRUE")
write.table(as.data.frame(data),file=file.path(out_dir,paste0("data",".txt")),sep=",")

data <- merge(data,dat_ID,by="id_coord",all=T,suffixes=c("","_y"))

#Process year by year:

#i<-10 #for 2012
for(i in 1:length(years_to_process)){
  year_processed <- years_to_process[i]
  in_dir_rst <- grep(paste0("prism_ppt_",year_processed), list_dir_rainfall,value=T)
  #in_dir_rst <- list_dir_rainfall[11]
  #start_date <- "2012-01-01" #PARAM 12
  start_date <- paste0(year_processed,"-01-01") #PARAM 12
  end_date <- paste0(year_processed,"-12-31") #PARAM 13 #should process by year!!!  
  #end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
  
  #debug(combine_stations_data_raster_ts_fun)
  df_combined <- combine_stations_data_raster_ts_fun(data,dat_stat,convert_to_inches,in_dir_rst,start_date,end_date,data_type,coord_names,out_dir,out_suffix)
  
}

# ###### Part 2: plot information ####
# 
# #Make this a function ??
# 
# ##Need to change here to match to the correct year...
# in_dir_rst <- list_dir_rainfall[11]
# 
# #
# 
# plot(r_rainfall,y=1)
# 
# res_pix <- 480
# col_mfrow <- 1
# row_mfrow <- 1
# 
# png(filename=paste("Figure1_","histogram","_station_coliform_measurements_frequency_",year_processed,"_",out_suffix,".png",sep=""),
#     width=col_mfrow*res_pix,height=row_mfrow*res_pix)
# 
# plot(table(data_subset$LOCATION_ID),type="h", main="Number of measurements",
#      ylab="Frequency of coliform measurements",xlab="Station ID")
# 
# dev.off()
# 
# #hist(table(data_subset$COL_SCORE))
# #Could output in a textfile
# range(table(dat_stat$LOCATION_ID)) #1 to 60
# mean(table(dat_stat$LOCATION_ID)) #8
# median(table(data_subset$LOCATION_ID)) #8
# hist(table(data_subset$LOCATION_ID), main="Frequency of Number of measurements by stations")
# 
# ref_pol <- create_polygon_from_extent(dat_stat,outDir=out_dir,outSuffix=out_suffix)
#   
# res_pix <- 480
# col_mfrow <- 1.5
# row_mfrow <- 1
# 
# png(filename=paste("Figure2a_","rainfall_map_",dates_l[1],"_and_stations_",out_suffix,".png",sep=""),
#     width=col_mfrow*res_pix,height=row_mfrow*res_pix)
# 
# plot(r_rainfall,y=1)
# plot(dat_stat,add=T)
# text(dat_stat,dat_stat$tID,cex=1.4)
# plot(ref_pol,border="red",add=T)
# legend("topright",legend=c("stations"), 
#        cex=1.2, col="black",pch =3,bty="n")
# 
# dev.off()
# 
# col_mfrow <- 2
# row_mfrow <- 1
# 
# png(filename=paste("Figure2b_","rainfall_map_",dates_l[1],"_and_stations_zoom_window_",out_suffix,".png",sep=""),
#     width=col_mfrow*res_pix,height=row_mfrow*res_pix)
# 
# plot(r_rainfall,y=1,ext=extent(dat_stat))
# plot(dat_stat,add=T,pch=3)
# text(dat_stat,dat_stat$tID,cex=1.4)
# legend("topright",legend=c("stations"), 
#        cex=1.2, col="black",pch =3,bty="n")
# 
# dev.off()
# 
# #
# #make a function later on?
# 
# #ADD selection function by ID
# #list_selected_pix: by LOCATION ID
# #
# 
# 
# #colnames(data_df) <- list_selected_ID
# #data_dz <- zoo(data_dz,idx)
# 
# plot(log(data_df$rainfall),log(data_df$COL_SCORE))
# plot(data_df$rainfall,data_df$COL_SCORE)
# 
# ## IDENTIFY 2 inches events?


############################ END OF SCRIPT #######################
