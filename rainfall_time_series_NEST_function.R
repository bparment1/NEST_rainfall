####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script contains functions to process and explore the correlation between rainfall events 
#and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/09/2015 
#DATE MODIFIED: 11/14/2015
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

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_functions.R" #PARAM 1
script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

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
  
  filenames<-filenames_all[pos]
  for(i in 1:length(filenames)){

    out_file<-basename(filenames[i])
    out_file_path<-file.path(out_path,out_file)
    download.file(filenames[i],destfile=out_file_path)
    #creating the file names
  }
}

plot_to_file <- function(raster_name,res_pix=480,out_suffix=NULL,out_dir=NULL){
  #Quick utility function to plot raster to png file for the workflow
  #This is useful to visually check the outputs from the workflow.
  #INPUT arguments:
  #raster_name: input raster object or file name of the raster object
  #out_suffix: output suffix
  #out_dir: output directory
  #OUTPUT:
  # png_file: name of the file containing the figure/plot
  #
  # Authors: Benoit Parmentier
  # Created: 11/01/2015
  # Modified: 11/02/2015
  # To Do: 
  # - Add option to choose plot format e.g. jpeg, tif etc.
  #
  ################################
  ## Create raster object if not already present
  if(class(raster_name)!="RasterLayer"){
    layer_rast<-raster(raster_name)
  }else{
    layer_rast <- raster_name
    raster_name <- filename(layer_rast)
  }
  
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  if(is.null(out_dir)){
    out_dir <- getwd()
  }
  
  #Extract name
  extension_str <- extension(raster_name)
  raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
  
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  #Might change file format later...
  png_file <- file.path(out_dir,paste("figure_",raster_name_tmp,"_",out_suffix,".png",sep="")) #may be changed later to other format
  png(filename=png_file,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  plot(layer_rast)
  title(raster_name_tmp)
  dev.off()
  return(png_file)
}


########################### End of script #####################################