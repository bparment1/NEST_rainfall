##############################################  NEST Beach closure project  #######################################
###########################################  Data preparation and download #######################################
#This script download data and processes rasters to fit the Maine study area.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 01/15/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -   
#          - 
#TO DO:
# -
#################################################################################################

###Loading R library and packages                                                      


library(sp)
library(raster)
library(rgdal)
require(rgeos)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
require(RCurl)
require(stringr)
require(XML)

###### Functions used in this script sourced from other files

#function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_12112015.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#script_path <- "/home/parmentier/Data/rainfall/NEST"
#source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.

##### Functions used in this script 

downloading_prism_product <- function(start_date, end_date,var_name,num_cores){
  #
  download_prism <- function(date_selected,var_name){
    url_product <- "http://services.nacse.org/prism/data/public/4km/" #URL is a constant...
    var_name <- "tmin" #tmax,ppt
    url_product_var <- paste(url_product,var_name,sep="")
    #system("wget --content-disposition $base_url/$clim_var/$day");
    cmd_str <- paste("wget --content-disposition ",url_product_var,date_selected,sep="")
    system(cmd_str)
    df <- data.frame(var=var_name,date=date_selected,url=url_product)
    return(df)
  }
  
  start_date <- as.Date(start_date,format="%Y-%m-%d") #start date
  end_date <- as.Date(end_date,format="%Y-%m-%d") #end date
  dates_range <- seq.Date(start_date, end_date, by="1 day") #sequence of dates
  dates_range_format <- format="%Y%m%d") #end date
  
  #20090405
  
  download_prism(dates_range[1],var_name)
  l_df_download <- mclapply(dates_range,FUN=download_prism,var_name=var_name,
                            mc.preschedule=FALSE,mc.cores = num_cores)
  
  #wget http://services.nacse.org/prism/data/public/4km/tmin/20090405
  #
  #
}

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
start_date <- "2012-01-01" #PARAM 12
end_date <- "2012-12-31" #PARAM 13

var_name <- "tmin"
num_cores <- 4
file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_01142016" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8


out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

debug(downloading_prism_product) 
test <- downloading_prism_product(start_date, end_date,var_name,num_cores)
  


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



