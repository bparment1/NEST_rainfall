##############################################  NEST Beach closure project  #######################################
#################################  Shiny Application for exploration of station measurements  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/10/2016 
#DATE MODIFIED: 03/28/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -   
#          - 
#TO DO:

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
library(shiny)

#library(shiny)
#runExample()
#runExample("07_widgets")

###### Functions used in this script sourced from other files

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_03272016.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#in_dir <- "/home/bparmentier/Dropbox/Data/NEST/NEST_stations_s02"
script_path <- "/home/bparmentier/Dropbox/Data/NEST/NEST_stations_s03" #path to script #PARAM 
setwd(script_path)
#script_path <- "." #path to script #PARAM 

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

#in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50 , param 1
in_dir <- "/home/bparmentier/Dropbox/Data/NEST/NEST_stations_s03"
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS, param 
#in_dir_rainfall <- "/home/bparmentier/Google Drive/NEST_Data/"

CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 #param 2
CRS_reg <- CRS_WGS84 # PARAM 3

file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_03282016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=FALSE #PARAM8
num_cores <- 4 #PARAM 9

#rainfall_dir <- "/home/bparmentier/Google Drive/NEST_Data" #PARAM 10
rainfall_dir <- "./data" #PARAM 10
#station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST_Data/", "WQ_TECS_Q.txt") #PARAM 11
#station_data_fname <- file.path("data", "MHB_data_2006-2015.csv") #PARAM 11

station_measurements_MHB_data_fname <- file.path("data", "data_df_rainfall_and_measurements_MHB.txt") #PARAM 11 
#This will change
station_measurements_DMR_data_fname <- file.path("data", "data_df_rainfall_and_measurements_MHB.txt") #PARAM 11 

start_date <- "2012-01-01" #PARAM 12, user define? this is the default value...
end_date <- "2012-12-31" #PARAM 13
#var_name <- "COL_SCORE" #PARAM 14
var_name <- "CONCENTRATION" #PARAM 14, MHB data, need to add DMR
var_ID <- "LOCATION_ID" #PARAM 15
year_processed <- "2012" #PARAM 16
threshold_val <- 2*25.4 #PARAM 17, in inches or mm
convert_to_inches <- FALSE #PARAM 18
units_val <- "mm"
data_type <- "MH" #for Maine beach health

## Change coordinates to x and y and lat long!!!
coord_names <- c("SITE.LONGITUDE..UTM.","SITE.LATITUDE..UTM.") #MH beach bacteria dataset
#coord_names <- c("LONGITUDE_DECIMAL","LATITUDE_DECIMAL") #cloroforms beach bacteria dataset
data_df_fname <- "./data/df_ts_pix_2012.txt"
#/home/bparmentier/Google Drive/NEST/NEST_stations_s02/data

SMAZones_fname <- "SMAZoneDissolve.shp"

dat_stat_location_MHB_fname <- "dat_stat_location_MHB.shp"
dat_stat_location_DMR_fname <- "dat_stat_location_DMR.shp"
dat_stat_location_DMR_fname <- "dat_stat_location_MHB.shp" #Use this for now

################# START SCRIPT ###############################

#set up the working directory
#Create output directory

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
setwd(in_dir)
#if(create_out_dir_param==TRUE){
#  out_dir <- create_dir_fun(out_dir,out_suffix_s)
#  setwd(out_dir)
#}else{
#  setwd(out_dir) #use previoulsy defined directory
#}

list_dir_rainfall <- list.dirs(path=rainfall_dir,full.names=T)
#remove non relevant directories
list_dir_rainfall <- grep("prism_ppt*", list_dir_rainfall,value=T)

#### Part 1: read in combined information by stations ####

data_df_MHB <- read.table(station_measurements_MHB_data_fname,sep=",",header=T,fill=T,stringsAsFactors = F) #bacteria measurements
data_df_DMR <- read.table(station_measurements_MHB_data_fname,sep=",",header=T,fill=T,stringsAsFactors = F) #bacteria measurements

dat_stat_location_MHB <- readOGR("./data",sub(".shp","",dat_stat_location_MHB_fname))
dat_stat_location_DMR <- readOGR("./data",sub(".shp","",dat_stat_location_MHB_fname))

### Part 2: read in raster rainfall data and SMA zones

SMAZones <- readOGR("./data",sub(".shp","",SMAZones_fname))

#for(i in 1:length(list_dir_rainfall)){
#  in_dir_rst <- list_dir_rainfall[11]
#  r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=in_dir_rst,full.names=T))) #rainfall time series stack
#}

if (convert_to_inches==TRUE){
  r_rainfall <- r_rainfall/25.4 #improve efficiency later? YES!!
}

#data_df <- read.table(data_df_fname)
#coord_names <- c("SITE.LONGITUDE..UTM.","SITE.LATITUDE..UTM.")
#idx <- seq(as.Date(start_date), as.Date(end_date), 'day')
#date_l <- strptime(idx[1], "%Y%m%d") # 
#dates_l <- format(idx, "%Y%m%d") #  date being processed

proj_str <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #This will need to be added in the parameters

#test <- over( SMAZones , ind_adm , fn = NULL) 
#test <- over( data_df  , SMAZones , fn = NULL) 

########################### End of script #####################################

#http://shiny.rstudio.com/tutorial/lesson7/
#http://shiny.rstudio.com/articles/shinyapps.html
