##############################################  NEST Beach closure project  #######################################
###########################################  Data preparation and download #######################################
#This script download data and processes rasters to fit the Maine study area.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 01/24/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - Running processing of prism for 2008 to 2009 included to checkt the bug!!
#          - 
#TO DO: DSS
# - make this callable from shell using optparse package with documentation with usage
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
require(RCurl)
require(stringr)
require(XML)

###### Functions used in this script sourced from other files

#function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_12112015.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#script_path <- "/home/parmentier/Data/rainfall/NEST"
#source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.

##### Functions used in this script 

create_dir_fun <- function(outDir,out_suffix=NULL){
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

## Add processing of downloaded files!!
## 1) crop and reproject if needed
## 2) creation of TimeRaster object with summaries?
## 3)
### Other functions ####

function_processing_NEST_data <- "processing_NEST_data_function_01232016c.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#script_path <- "/home/parmentier/Data/rainfall/NEST"
source(file.path(script_path,function_processing_NEST_data)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50 , param 1
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS, param 
out_dir <- "/home/bparmentier/Google Drive/NEST/" #param 2

start_date <- "2008-01-01" # param 3
end_date <- "2009-12-31" # param 4

var_name <- "ppt" #tmin,tmax #param 5
num_cores <- 4 #param 6
file_format <- ".tif" #param 7
NA_value <- -9999 # param 8
NA_flag_val <- NA_value 
out_suffix <-"NEST_prism_01232016" #output suffix for the files and ouptu folder #param 9
create_out_dir_param=TRUE # param 10
ref_rast_name <- "/home/bparmentier/Google Drive/NEST/prism_rain/prismrain2012/prismrain_20120101.tif" #param 10
agg_param <- c(FALSE,NULL,"mean") #False means there is no aggregation!!! #param 11

######### PART 0: Set up the output dir ################

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

########## PART 1: Download and unzip ##############

#turn this into a function and download and process by year
# do 2003 to 2014

#debug(download_and_process_prism_data) 
lf_reg <- download_and_process_prism_data(in_dir,out_dir,start_date,end_date,var_name,ref_rast_name,
                      agg_param,num_cores=1,create_out_dir_param=FALSE,NA_value=-9999,out_suffix="",
                      file_format=".tif",download_file=T,unzip_files=T,match_file=T,lf_zip=NULL,lf_r=NULL)

  

#################### END OF SCRIPT #####################




