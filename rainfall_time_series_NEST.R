####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 11/05/2015
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

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_05172015_functions.R" #PARAM 1
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

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

format_s <-function(s_ID){
  #Format station ID in a vector format/tuple that is used in a psql query.
  # Argument 1: vector of station ID
  # Return: character of station ID
  tx2<-s_ID
  tx2<-as.character(tx2)
  stat_list<-tx2
  temp<-shQuote(stat_list)
  t<-paste(temp, collapse= " ")
  t1<-gsub(" ", ",",t)
  sf_ID<-paste("(",t1,")",sep="") #vector containing the station ID to query
  return(sf_ID)
}

download_files_with_pattern_ftp<-function(url_dir,out_path,file_pattern){
  #sadfd
  #adsfd
  #library(RCurl) #modify later using require
  # FTP 
  url<-url_dir
  filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE) 
  file_pattern_rx<-glob2rx(file_pattern)
  # Deal with newlines as \n or \r\n. (BDR) 
  # Or alternatively, instruct libcurl to change \n's to \r\n's for us with crlf = TRUE 
  # filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE) 
  
  filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "") 
  filenames_all<-filenames
  #Now subset filenames to match wanted files...
  i_file<-grep(file_pattern_rx,filenames_all,value=TRUE) # using grep with "value" extracts the matching names
  pos<-match(i_file,filenames_all)
  
  filenames<-filenames_all[pos]53)})) 
  for(i in 1:length(filenames)){
    d
    out_file<-basename(filenames[i])
    out_file_path<-file.path(out_path,out_file)
    download.file(filenames[i],destfile=out_file_path)
    #creating the file names
  }
}

##### STEP 1: Select station in the study area

filename<-sub(".shp","",infile1)             #Removing the extension from file.
interp_area <- readOGR(".",filename)
CRS_interp<-proj4string(interp_area)         #Storing the coordinate information: geographic coordinates longlat WGS84

#Read in GHCND database station locations
dat_stat <- read.fwf(infile2, 
                     widths = c(11,9,10,7,3,31,4,4,6),fill=TRUE)
colnames(dat_stat)<-c("STAT_ID","lat","lon","elev","state","name","GSNF","HCNF","WMOID")
coords<- dat_stat[,c('lon','lat')]
coordinates(dat_stat)<-coords
proj4string(dat_stat)<-CRS_locs_WGS84 #this is the WGS84 projection
#proj4string(dat_stat)<-CRS_interp
dat_stat2<-spTransform(dat_stat,CRS(CRS_interp))         # Project from WGS84 to new coord. system

# Spatial query to find relevant stations
inside <- !is.na(over(dat_stat2, as(interp_area, "SpatialPolygons")))  #Finding stations contained in the current interpolation area
stat_reg<-dat_stat2[inside,]              #Selecting stations contained in the current interpolation area

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
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

