##############################################  NEST Beach closure project  #######################################
###########################################  Figures and Analyses production  #######################################
#This script prepares data for analyes of rainfall events and beach closures due to bacteria outbreaks 
#in Maine. There are two datasets provided containing information on bacteria measurements:
#DMR:
#MHB:
#The rainfall data is currently obtained from PRISM.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 06/29/2016
#Version: 5
#PROJECT: NEST beach closures            
#
#
#COMMENTS: 
#          - Add spacetime object functions for later
#          - Add datatype in name
#TO DO:
# - Select 2 inches rainfall events and correlates with bacteria data
# - Compute accumulated rain over several days using time series functions
# - Make a movie sequence later on using animation package in R
#

###Stages:
#Stage 1: downloading_NEST_data_01232016.R
#Stage 1: processing_NEST_data_06272016.R
#Stage 1: processing_NEST_data_function_06272016.R
#Stage 2: rainfall_time_series_NEST_06292016.R
#Stage 2: rainfall_time_series_NEST_function_06292016.R
#Stage 3: Visualization and exploration Shiny: UI_04262016.R
#Stage 3: Visualization and exploration Shiny: global_04262016.R
#Stage 3: Visualization and exploration Shiny: server_04262016.R
#Stage 4: Analyses: analyses_NEST_06232016b.R
#Stage 4: Analyses: analyses_NEST_functions_06232016b.R

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

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_06292016.R" #PARAM 1
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
out_suffix <-"NEST_prism_06292016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8
num_cores <- 4 #PARAM 9

rainfall_dir <- "/home/bparmentier/Google Drive/NEST_Data" #PARAM 10
station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST_Data/", "WQ_TECS_Q.txt") #PARAM 11,DMR
#station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST/", "MHB_data_2006-2015.csv") #PARAM 11

years_to_process <- 2003:2016
#start_date <- "2012-01-01" #PARAM 12
#end_date <- "2012-12-31" #PARAM 13 #should process by year!!!
var_name <- "COL_SCORE" #PARAM 14, Name of variable of interest: bacteria measurement (DMR data)
#var_name_raster <- rainfall
#var_name <- "CONCENTRATION" #PARAM 14, MHB data
var_ID <- "LOCATION_ID" #PARAM 15 #this is for DMR data
#var_ID <- NULL #PARAM 15 if null then create a new ID for stations!!, not null for DMR
year_processed <- "2012" #PARAM 16
threshold_val <- 2*25.4 #PARAM 17, in inches or mm
convert_to_inches <- FALSE #PARAM 18
units_val <- "mm" #PARAM 19
#data_type <- "MHB" #for Maine Healthy Beaches
data_type <- "DMR" #for Maine Department of Marine Resources #PARAM 20
x_coord_range <- c(-180,180) #PARAM 21
y_coord_range <- c(-90,90) #PARM 22
ref_raster_name <- NULL#
## Add ref_raster_name ,if not null use that for the coordinates limit of the study area

#coord_names <- c("SITE.LONGITUDE..UTM.","SITE.LATITUDE..UTM.") #MHB, beach bacteria dataset
coord_names <- c("LONGITUDE_DECIMAL","LATITUDE_DECIMAL") #DMR, cloroforms beach bacteria dataset #PARAM 23

SMAZones_fname <- "/home/bparmentier/Google Drive/NEST/NEST_stations_s02/data/SMAZoneDissolve.shp" #PARAM 24
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

###########################
#### Part 1: read in and combine the information ####

data <- read.table(station_data_fname,sep=",",header=T,fill=T,stringsAsFactors = F) #bacteria measurements
data$FID <- 1:nrow(data)

##################
##Make this part a function later on...

### Select unique stations and make a SPDF
#data$FID <- 1:nrow(data)
data[[var_name]] <- as.numeric(data[[var_name]]) #make sure we have a numeric field
data[[coord_names[1]]] <- as.numeric(data[[coord_names[1]]])
data[[coord_names[2]]] <- as.numeric(data[[coord_names[2]]])
#Screen for coordinates that are out of the range!!
index_x <- data[[coord_names[1]]] < x_coord_range[1]
data[[coord_names[1]]][index_x] <- NA 
#Screen for coordinates that are out of the range!!
index_x <- data[[coord_names[2]]] > x_coord_range[2]
data[[coord_names[2]]][index_x] <- NA 
index_y <- data[[coord_names[2]]] < y_coord_range[1]
data[[coord_names[2]]][index_y] <- NA 
index_y <- data[[coord_names[2]]] > y_coord_range[2]
data[[coord_names[2]]][index_y] <- NA 

#> dim(data)
#[1] 106262     36

##Now remove all the NA in coordinates (data were set NA because out of range or already NA)
data <- subset(data, !is.na(data[[coord_names[1]]]) & !is.na(data[[coord_names[2]]]))

#> dim(data)
#[1] 83187    37

#### Now check that dates are available for all records!!

#MHB data from Maine Healthy Beaches
if(data_type=="MHB"){
  dates_TRIP_START <- unlist(lapply(strsplit(data$SAMPLE.DATE," "),function(x){x[1][1]}))
  #dates_TRIP_START <- gsub(" 0:00:00","",data$SAMPLE.DATE)
  data$TRIP_START_DATE_f <- as.Date(strptime(dates_TRIP_START,"%m/%d/%Y"))
  data$TRIP_START_DATE_year <- four.digit.year(data$TRIP_START_DATE_f , year=1968)
  #format(data$TRIP_START_DATE_f , format="%m-%d-%Y")
  data$TRIP_START_DATE_month <- strftime(data$TRIP_START_DATE_f , "%m") # current month of the date being processed
  data$TRIP_START_DATE_day <- strftime(data$TRIP_START_DATE_f , "%d")
}
#DMR data from Deparment of Marine Resources
if(data_type=="DMR"){
  dates_TRIP_START <- gsub(" 0:00:00","",data$TRIP_START_DATE)
  data$TRIP_START_DATE_f <- as.Date(strptime(dates_TRIP_START,"%m/%d/%Y"))
  data$TRIP_START_DATE_month <- strftime(data$TRIP_START_DATE_f , "%m") # current month of the date being processed
  data$TRIP_START_DATE_year <- strftime(data$TRIP_START_DATE_f , "%Y")
  data$TRIP_START_DATE_day <- strftime(data$TRIP_START_DATE_f , "%d")
}

data$TRIP_START_DATE_f <- paste0(data$TRIP_START_DATE_year,data$TRIP_START_DATE_month,data$TRIP_START_DATE_day)
data$TRIP_START_DATE_f <- as.Date(strptime(data$TRIP_START_DATE_f,"%Y%m%d"))

##### Now remove data that have no dates in the measurements!
##for DMR removing 7 rows
data <- subset(data, !is.na(data[["TRIP_START_DATE_f"]]) )

#> dim(data)
#[1] 83180    41

#### FIRST SCREEN BASED ON THE ID

if(data_type=="DMR"){
  #remove items that do not follow "EA007.00" i.e. 8 characters
  
  index_ID <- (nchar(data[[var_ID]]) != 8)
  #index_y <- data[[coord_names[2]]] > y_coord_range[2]
  #test <- data
  #test[[var_ID]][index_ID] <- NA 
  data[[var_ID]][index_ID] <- NA 
  #index_y <- data[[coord_names[2]]] > y_coord_range[2]
  #data[[coord_names[2]]][index_y] <- NA 
  ##Now remove all the NA in coordinates (data were set NA because out of range or already NA)
  data <- subset(data, !is.na(data[[var_ID]]))
  #> dim(data)
  #[1] 83170    42
}

### SECOND SCREEEN BASED ON COORDINATES

## Remove duplicates rows from stations to identify uniques sations
##Create unique identifier from coordinates first

data$id_coord <- paste(data[[coord_names[1]]],data[[coord_names[2]]],sep="_")

id_coord <- unique(data$id_coord)
ID_stat <- 1:length(id_coord)
dat_ID <- data.frame(id_coord,ID_stat)
data <- merge(data,dat_ID,by="id_coord",all=T,suffixes=c("","_y"))
#data_merged <- merge(data,dat_stat,by="FID",all=T,suffixes=c("","_y"))
coords <- data[,coord_names]
coordinates(data) <- coords  

##Create a unique set of items based on location only
dat_stat <- remove.duplicates(data[,c("ID_stat","id_coord",coord_names[[1]],coord_names[[2]])])
#> dim(dat_stat)
#[1] 2245    4
#> length(unique(dat_stat$ID_stat))
#[1] 2245
#> length(unique(data$LOCATION_ID))
#[1] 2238
#> length(unique(data$ID_stat))
#[1] 2238

if(data_type=="MHB"){
  proj_str <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #This will need to be added in the parameters
  proj4string(dat_stat) <- proj_str #this is the NAD83 latitude-longitude
  dat_stat<-spTransform(dat_stat,CRS(CRS_WGS84))     #Project from WGS84 to new coord. system
}
if(data_type=="DMR"){
  proj4string(dat_stat) <- CRS_WGS84 #this is the NAD83 latitude-longitude
}

## Remove duplicates rows from stations to identify uniques sations

##MODIFY HERE, FIND MATCHES between LOCATION_ID and ID_stat
#2251 for ID_stat
#2247 unique val for LOCATION_ID
#sort(unique(data$LOCATION_ID))
#> xtabs(data$ID_stat ~data$LOCATION_ID,data)
#data$LOCATION_ID
#2.00     25.10     29.00     31.00      4.00     86.00     A1COL  EA001.00  EA002.00  EA003.00  EA003.50 
#1640      1945      1597      1703      1795      1700      1924     32886    108813     90882     16464 
#EA004.00  EA005.00  EA006.00  EA007.00  EA008.00  EA009.00  EA010.00  EA011.00  EA012.00  EA013.00  EA014.00 
#53380     33798      3888     33054     33376     28368     41740     42606     52930     58534     60984 
#There are some issues with the station identifier. FIx this later, use ID_stat as the identifier 
#for processing for now.

if(is.null(var_ID)){
  dat_stat$LOCATION_ID <- dat_stat$ID_stat
  var_ID <- "LOCATION_ID"
}#else{
#  test <- merge(dat_stat,data,by="id_coord",all=T,suffixes=c("","_y"))
#}
#crosstab(Survey, row.vars = "Age", col.vars = "Sex", type = "f")

#xtabs(data$ID_stat,data$LOCATION_ID)
#xtabs(data$ID_stat ~data$LOCATION_ID,data)
### Write out basic informaiotn before processing by years

write.table(as.data.frame(dat_stat),file=file.path(out_dir,paste0("dat_stat_location_",data_type,".txt")),sep=",")
writeOGR(dat_stat,dsn= out_dir,layer= paste0("dat_stat_location_",data_type), 
          driver="ESRI Shapefile",overwrite_layer="TRUE")
write.table(as.data.frame(data),file=file.path(out_dir,paste0("data_measurements_",data_type,".txt")),sep=",")

#data <- merge(data,dat_ID,by="id_coord",all=T,suffixes=c("","_y"))

##########################
#### Part2: Process year by year

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
  #df_combined <- combine_stations_data_raster_ts_fun(data,dat_stat,convert_to_inches,in_dir_rst,start_date,end_date,data_type,coord_names,out_dir,out_suffix)
  df_combined <- combine_stations_data_raster_ts_fun(data,dat_stat,convert_to_inches,in_dir_rst,start_date,end_date,year_processed,data_type,coord_names,num_cores,out_dir,out_suffix)
    
}
#For DMR data, it took 1h43 minutes...

###########################
#### Part 3: combine all data together

#Don't combine ? It's too large for DMR: 676002 rows per year or 13*676002=8,788,026
#list_df_fname <- list.files(path=out_dir,pattern="data_df_.*._NEST_prism_03272016.txt",full.names=T)
#list_df <- lapply(list_df_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
#tb <- do.call(rbind,list_df)

#write.table(tb,file=file.path(out_dir,paste0("data_df_rainfall_and_measurements_",data_type,".txt")),sep=",")

############################ END OF SCRIPT #######################
