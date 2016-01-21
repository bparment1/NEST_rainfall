##############################################  NEST Beach closure project  #######################################
###########################################  Data preparation and download #######################################
#This script download data and processes rasters to fit the Maine study area.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 01/21/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -   
#          - 
#TO DO:
# -make this callable from shell?
#
#################################################################################################

###Loading R library and packages                                                      

library(parallel)
library(sp)
library(raster)
library(rgdal)
require(rgeos)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
require(RCurl)
require(stringr)
require(XML)
library(lubridate)

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

download_prism <- function(date_selected,var_name,out_dir){
  current_dir <- getwd()
  setwd(out_dir)
  url_product <- "http://services.nacse.org/prism/data/public/4km" #URL is a constant...
  #var_name <- "tmin" #tmax,ppt
  #url_product_var <- paste(url_product,var_name,sep="")
  url_product_var <- file.path(url_product,var_name,date_selected)
  #system("wget --content-disposition $base_url/$clim_var/$day");
  cmd_str <- paste("wget --content-disposition ",url_product_var,sep="")
  system(cmd_str)
  df <- data.frame(var=var_name,date=date_selected,url=url_product)
  setwd(current_dir)
  return(df)
}

downloading_prism_product <- function(start_date, end_date,var_name,num_cores,out_dir=NULL){
  #Function to download prism dataset for range of dates and a type of variable(tmin,tmax,precip).
  #Inputs
  #1) start_date: start of date to download
  #2) end_date: end date to download
  #3) var_name: type of variable i.e. tmin, tmax or precip
  #4) num_cores: number of cores used to download
  #5) out_dir: directory where to place downloaded files
  #Output: 
  #df_downloaded: data.frame of downloaded files
  #
  #TO DO: add option for different url product (other than 4km daily)
  #
  
  ##### Start of Script #####
  
  start_date <- as.Date(start_date,format="%Y-%m-%d") #start date
  end_date <- as.Date(end_date,format="%Y-%m-%d") #end date
  dates_range <- seq.Date(start_date, end_date, by="1 day") #sequence of dates
  #dates_range_format <- as.Date(dates_range,format="%Y%m%d") #end date
  
  ##Check if the range is over multiple years
  date_year <- strftime(dates_range, "%Y")
  date_month <- strftime(dates_range , "%m") # current month of the date being processed
  date_day <- strftime(dates_range , "%d")
  dates_range_prism_format <- paste(date_year,date_month,date_day,sep="")
  #20090405
  
  if (is.null(out_dir)){
    out_dir <- file.path(getwd(),paste("prism_",var_name,date_year,sep=""))
    create_dir_fun(out_dir)
  }

  #s_obj <- download_prism(dates_range_prism_format[1],var_name)
  l_df_download <- mclapply(dates_range_prism_format,
                            FUN=download_prism,
                            var_name=var_name,
                            out_dir=out_dir,
                            mc.preschedule=FALSE,
                            mc.cores = num_cores)
  #Currently we use only this:
  #wget http://services.nacse.org/prism/data/public/4km/tmin/20090405
  
  df_downloaded <- do.call(rbind, l_df_download)
  df_downloaded$date <- as.character(df_downloaded$date)
  lf_zip <- unlist(lapply(df_downloaded$date,function(x){list.files(pattern=paste(x,".*.zip$",sep=""),path=out_dir)}))
  df_downloaded$file_zip <- lf_zip
  df_downloaded$dir <- out_dir
  return(df_downloaded)
}

## Add processing of downloaded files!!
## 1) crop and reproject if needed
## 2) creation of TimeRaster object with summaries?
## 3)

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/NEST/" #local bpy50 , param 1
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS, param 
out_dir <- "/home/bparmentier/Google Drive/NEST/prism_tmin"
start_date <- "2011-01-01" #PARAM 12
end_date <- "2011-12-31" #PARAM 13

var_name <- "ppt" #tmin,tmax
num_cores <- 4
file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_01142016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=FALSE #PARAM8

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#debug(downloading_prism_product) 
test <- downloading_prism_product(start_date, end_date,var_name,num_cores,out_dir=NULL)
  
##### Now unzip function
unzip_lf_fun <- function(lf, out_dir){
  #unzip PRISM_ppt_stable_4kmD2_20110102_bil.zip
  unzip(lf,exdir=out_dir)
}

la

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



