##############################################  NEST Beach closure project  #######################################
###########################################  Data preparation and download #######################################
#This script download data and processes rasters to fit the Maine study area.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 01/23/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -   
#          - 
#TO DO: DSS
# -make this callable from shell?
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

function_processing_NEST_data <- "processing_NEST_data_function_01232016.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#script_path <- "/home/parmentier/Data/rainfall/NEST"
source(file.path(script_path,function_processing_NEST_data)) #source all functions used in this script 1.


############################################################################
#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50 , param 1
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS, param 
out_dir <- "/home/bparmentier/Google Drive/NEST/" #param 2

start_date <- "2003-12-17" # param 3
end_date <- "2004-01-05" # param 4

var_name <- "ppt" #tmin,tmax #param 5
num_cores <- 4 #param 6
file_format <- ".tif" #param 7
NA_value <- -9999 # param 8
NA_flag_val <- NA_value 
out_suffix <-"NEST_prism_01212016" #output suffix for the files and ouptu folder #param 9
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



#debug(downloading_prism_product) 


download_and_process_prism_data <- function(in_dir,out_dir,start_date,end_date,var_name,num_cores,file_format,
         NA_value,out_suffix,create_out_dir_param,ref_rast_name,agg_param,download_file=T,unzip_files=T,match_file=T,lf_zip=NULL){
  
  #INPUTS
  #1)in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50 , param 1
  #2)out_dir <- "/home/bparmentier/Google Drive/NEST/" #param 2
  #3)start_date <- "2003-12-17" # param 3
  #4)end_date <- "2004-01-05" # param 4
  #5)var_name <- "ppt" #tmin,tmax #param 5
  #6)num_cores <- 4 #param 6
  #7)file_format <- ".tif" #param 7
  #8)NA_value <- -9999 # param 8
  #9)out_suffix <-"NEST_prism_01212016" #output suffix for the files and ouptu folder #param 9
  #10)create_out_dir_param=TRUE # param 10
  #11)ref_rast_name <- "/home/bparmentier/Google Drive/NEST/prism_rain/prismrain2012/prismrain_20120101.tif" #param 10
  #12)agg_param <- c(FALSE,NULL,"mean") #False means there is no aggregation!!! #param 11
  
  
  NA_flag_val <- NA_value 
  
  ######### PART 1 ######
  if(download_file==T){
    setwd(in_dir)
    download_obj <- downloading_prism_product(start_date, end_date,var_name,num_cores,out_dir=NULL)
    setwd(out_dir)
  }else{
    download_obj<-try(load_obj(file.path(out_dir,paste("downloaded_obj _",var_name,".RData",sep=""))))
  }
  
  
  ####### PART 2 #######

  ## run by year!!
  ## loop through year...or make this a function?
  if(zip_file==T){
    if(!is.null(lf_zip)){
      lf_r <- lapply(lf_zip, unzip,exdir= out_dir)
      lf_r <- list.files(pattern="*bil.bil$",path=out_dir,full.names = T)
    }else{
      lf_zip <- download_obj$lf_zip
      lf_r <- lapply(lf_zip, unzip,exdir= out_dir)
      lf_r <- list.files(pattern="*bil.bil$",path=out_dir,full.names = T)
    }
  }
  
  ####### PART 3 #######  
  
  list_lf_r_reg <- vector("list",length(test))
  for(i in 1:length(test)){
    file_zip_year <- test[[i]]$file_zip
    out_dir_year <- unique((test[[i]]$dir))
    lf_zip <- file.path(out_dir_year,file_zip_year)
    lf_r <- lapply(lf_zip, unzip,exdir= out_dir_year)
    lf_r <- list.files(pattern="*bil.bil$",path=out_dir_year,full.names = T)
    
    ########## PART 2: Match to study area ##############
    
    ## Match to the study area...
    ## Now crop and reproject if necessary
    #Use function above
    #
    r_ref <- raster(ref_rast_name)
    plot(r_ref)
    agg_param <- c(FALSE,NULL,"mean") #False means there is no aggregation!!!
    #use r_ref as reference...
    out_rast_name <- NULL
    list_param_create_region <- list(as.list(lf_r),
                                     r_ref, out_rast_name,agg_param,
                                     file_format,NA_flag_val,
                                     input_proj_str=NULL,out_suffix="",out_dir_year)
    names(list_param_create_region) <- c("raster_name",
                                         "reg_ref_rast", "out_rast_name","agg_param",
                                         "file_format","NA_flag_val",
                                         "input_proj_str","out_suffix","out_dir")
    #debug(create__m_raster_region)
    #test <- create__m_raster_region(1,list_param=list_param_create_region)
    
    #test <- create__m_raster_region(1,list_param=list_param_create_region)
    lf_r_reg <- mclapply(1:length(lf_r),
                         FUN=create__m_raster_region,
                         list_param=list_param_create_region,
                         mc.preschedule=FALSE,
                         mc.cores = num_cores)
    #Currently we use only this:
    
    #rainfall <- stack(unlist(lf_r_reg ))
    #plot(rainfall)
    list_lf_r_reg[[i]] <- lf_r_reg
  }
  
  
}
  

########## PART 3: Match to study area ##############



#################### END OF SCRIPT #####################

# use strict;
# use warnings;
# use DateTime;
# day <- 
# my $clim_var = 'ppt';
# my $base_url = 'http://services.nacse.org/prism/data/public/4km';
# my $start = DateTime->new( day => 1, month => 10, year => 1999 );
# my $stop = DateTime->new( day => 30, month => 9, year => 2000 );
# while($start <= $stop) {
#   $day = $start->strftime('%Y%m%d'); #place date in proper form
#   system("wget --content-disposition $base_url/$clim_var/$day");
#   sleep 2; #to be nice to our server
#   $start->add(days => 1);
# }
# 
# 
# dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates
# 



