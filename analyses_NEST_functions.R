##############################################  NEST Beach closure project  #######################################
#################################  Analyses for exploration of station measurements  #######################################
#This contains functions used to explore the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#It uses time series processed earlier in R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 06/15/2016 
#DATE MODIFIED: 06/15/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - Fixing problem with extent of input stations data 

#
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

###### Functions used in this script sourced from other files

#function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_03272016.R" #PARAM 1
#script_path <- "." #path to script #PARAM 

#script_path <- "/home/bparmentier/Google Drive/NEST/NEST_stations_s06/" #path to script #PARAM 
#script_path <- "/home/benoit/data/NEST_stations_s06" #on SSI server
#setwd(script_path)

#script_path <- "/home/parmentier/Data/rainfall/NEST"
#source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.

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

read_select_station <- function(file_name,selected_val,selected_col){
   tb <- read.table(file_name,sep=",")
   df <- subset(tb,tb[[selected_col]]==selected_val)
   #df <- subset(tb,tb$ID_stat==selected_ID)
   return(df)
}

run_simple_lm <- function(df,y_var_name,x_var_name,log_val=T){
  #
  #
  #test <- lm(log1p(df$COL_SCORE) ~ log1p(df$rainfall),df)
  if(log_val==T){
    mod <- lm(log1p(df[[y_var_name]]) ~ log1p(df[[x_var_name]]),df)
  }
  if(log_val!=T){
    mod <- lm(df[[y_var_name]] ~ df[[x_var_name]],df)
  }

  #plot(log1p(df[[y_var_name]]) ~ log1p(df[[x_var_name]])) #only two data points for 2003!!!
  summary_mod_tb <- (summary(mod))
  tb_coefficients <- as.data.frame(summary_mod_tb$coefficients)
  #plot(mod)
  #change name of rows and add n columns for the number of inputs in the model!

  return(tb_coefficients)
}

########################### End of script #####################################
