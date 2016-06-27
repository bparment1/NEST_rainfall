##############################################  NEST Beach closure project  #######################################
###########################################  Data preparation and download #######################################
#This script download data and processes rasters to fit the Maine study area.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 06/27/2016
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

downloading_prism_product <- function(start_date, end_date,var_name,num_cores,in_dir=".",out_dir=NULL){
  #Function to download prism dataset for range of dates and a type of variable(tmin,tmax,precip).
  #Inputs
  #1) start_date: start of date to download
  #2) end_date: end date to download
  #3) var_name: type of variable i.e. tmin, tmax or precip
  #4) num_cores: number of cores used to download
  #5) in_dir: input directory
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
  
  df_dates_range <- data.frame(dates_range_prism_format=dates_range_prism_format,
                               dates_range=dates_range,
                               year=date_year)
  #20090405
  ##Now subset...by year
  list_year <- unique(date_year)
  
  ##First create output dir by year!
  list_out_dir <- vector("list",length=length(list_year))
  if (is.null(out_dir)){
    for(i in 1:length(list_year)){
      out_dir_tmp <- file.path(in_dir,paste("prism_",var_name,"_",list_year[i],sep=""))
      create_dir_fun(out_dir_tmp)
      list_out_dir[[i]] <- out_dir_tmp
    }
    names(list_out_dir) <- list_year
  }

  ## loop through year..and download
  list_df_downloaded <- vector("list",length=length(list_year))
  for(i in 1:length(list_year)){
    year_tmp <- list_year[[i]] #year being dowloaded
    df_tmp <- subset(df_dates_range,df_dates_range$year==year_tmp)
    df_tmp$dates_range_prism_format <- as.character(df_tmp$dates_range_prism_format)
    #s_obj <- download_prism(dates_range_prism_format[1],var_name)
    out_dir_tmp <- list_out_dir[[i]]
    dates_range_prism_format_tmp <- df_tmp$dates_range_prism_format
    l_df_download <- mclapply(dates_range_prism_format_tmp,
                              FUN=download_prism,
                              var_name=var_name,
                              out_dir=out_dir_tmp,
                              mc.preschedule=FALSE,
                              mc.cores = num_cores)
    #Currently we use only this:
    #wget http://services.nacse.org/prism/data/public/4km/tmin/20090405
    
    df_downloaded <- do.call(rbind, l_df_download)
    df_downloaded$date <- as.character(df_downloaded$date)
    lf_zip <- unlist(lapply(df_downloaded$date,function(x){list.files(pattern=paste(x,".*.zip$",sep=""),
                                                                      path=out_dir_tmp)}))
    df_downloaded$year <- list_year[[i]]
    df_downloaded$file_zip <- lf_zip
    df_downloaded$dir <- out_dir_tmp
    df_downloaded_fname <- file.path(out_dir_tmp,paste("df_downloaded","_",var_name,"_",year_tmp,".txt",sep=""))
    write.table(df_downloaded,file=df_downloaded_fname,sep=",")
    list_df_downloaded[[i]] <- df_downloaded
  }
  ##
  names(list_df_downloaded)<- paste("year_",list_year,sep="")
  
  list_file_zip_year <- lapply(list_df_downloaded,function(x){x$file_zip})
  list_out_dir_year <- unique(lapply(list_df_downloaded,function(x){x$dir}))
  #out_dir_year <- unique((test[[i]]$dir))
  #file.path(list_out_dir_year,list_file_zip_year)
  lf_zip <- lapply(1:length(list_out_dir_year),function(i,x,y){file.path(x[[i]],y[[i]])},x=list_out_dir_year,y=list_file_zip_year)
  names(lf_zip) <- paste("year_",list_year,sep="")
  
  downloaded_obj <- list(list_df_downloaded,lf_zip)
  names(downloaded_obj) <-c("list_df","lf_zip")
  save(downloaded_obj,file= file.path(in_dir,paste("downloaded_prism_data_",var_name,".RData",sep="")))
  return(downloaded_obj)
}

#This function creates a spatial polygon data frame object for the extent matching a raster input
create_polygon_from_extent<-function(reg_ref_rast,outDir=NULL,outSuffix=NULL){
  #This functions returns polygon sp from input rast
  #Input Arguments: 
  #reg_ref_rast: input ref rast
  #outDir : output directory, if NULL then the current dir in used
  #outSuffix: output suffix used for the naming of the shapefile
  #Output: 
  #reg_outline_poly: spatial polygon data.frame
  #
  if(is.null(outDir)){
    outDir=getwd()
  }
  if(is.null(outSuffix)){
    outSuffix=""
  }
  ref_e <- extent(reg_ref_rast) #extract extent from raster object
  reg_outline_poly <- as(ref_e, "SpatialPolygons") #coerce raster extent object to SpatialPolygons from sp package 
  reg_outline_poly <- as(reg_outline_poly, "SpatialPolygonsDataFrame") #promote to spdf
  proj4string(reg_outline_poly) <- projection(reg_ref_rast) #Assign projection to spdf
  infile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="") #name of newly crated shapefile with the extent
  writeOGR(reg_outline_poly,dsn= outDir,layer= sub(".shp","",infile_reg_outline), 
           driver="ESRI Shapefile",overwrite_layer="TRUE")
  
  return(reg_outline_poly) #return spdf
}

#Function to aggregate from fine to coarse resolution, this will change accordingly once the input raster ref is given..
#
aggregate_raster <- function(agg_fact,r_in,reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL){
  #Aggregate raster from raster input and reference file
  #INPUT arguments:
  #agg_fact: factor to aggregate
  #agg_fun: default is mean
  #out_suffix: output suffix
  #file_Format: raster format used e.g. .tif
  #reg_ref_rast: reference raster to match in resolution, if NULL then send a message
  #out_dir: output directory
  #OUTPUT:
  # raster_name: name of the file containing the aggregated raster
  #
  # Authors: Benoit Parmentier
  # Created: 10/15/2015
  # Modified: 11/01/2015
  # To Do: 
  # - Add option to disaggregate
  #
  ################################
  
  if(is.null(agg_fact)){
    res_ref <- res(reg_ref_rast)[1] #assumes square cells, and decimal degrees from WGS84 for now...
    res_in <- res(r_in)[1] #input resolution, assumes decimal degrees
    agg_fact <-round(res_ref/res_in) #find the factor needed..
    #fix this to add other otpions e.g. aggregating down
  }
  
  #Default values...
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  raster_name <- file.path(out_dir,paste("r_agg_",agg_fact,out_suffix,file_format,sep="")) #output name for aster file
  r_agg <- aggregate(r_in, fact=agg_fact,FUN=agg_fun,filename=raster_name,overwrite=TRUE)
  
  return(raster_name)
  
}

### This is a very general function to process raster to match a raster of reference,
create__m_raster_region <-function(j,list_param){
  #This processes a list of raster to match to a region of interest defined by a reference raster
  #INPUT Arguments: raster name of the file,reference file with
  # j: file to be processed with input parameters
  # raster_name: list of raster to process i.e. match the region of interest
  # reg_ref_rast: reference raster used to defined the region of interest and spatial parameters
  # out_rast_name: output raster name, if NULL then use out_suffix to the input name to generate output name
  # agg_param: aggregation parameters: this is a vector used in in the aggregate function. It has three items:
  #                                     -TRUE/FALSE: if true then aggregate
  #                                     -agg_fact: aggregation factor, if NULL compute on the fly
  #                                     -agg_fun: aggregation function to use, the default is mean
  # file_format: output format used in the raster e.g. .tif, .rst
  # NA_flag_val: flag value used for no data
  # input_proj_str: defined projection,default null in which case it is extract from the input raster
  # out_suffix : output suffix added to output names if no output raster name is given
  # out_dir:  <- list_param$out_dir
  # Output: spatial grid data frame of the subset of tiles
  #
  # Authors: Benoit Parmentier
  # Created: 10/01/2015
  # Modified: 01/20/2016
  #TODO:
  # - Add option to disaggregate...
  # - Modify agg param to be able to use different ones by file j for the mcapply function
  #
  ################################################
  ## Parse input arguments
  raster_name <- list_param$raster_name[[j]] #list of raster ot project and crop, this is a list!!
  reg_ref_rast <- list_param$reg_ref_rast #This must have a coordinate system defined!!
  out_rast_name <- list_param$out_rast_name[j] #if NULL then use out_suffix to add to output name
  agg_param <- list_param$agg_param #TRUE,agg_fact,agg_fun
  file_format <- list_param$file_format #.tif, .rst
  NA_flag_val <- list_param$NA_flag_val #flag value used for no data
  input_proj_str <- list_param$input_proj_str #default null?
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  ## Start #
  
  ## Create raster object if not already present
  if(class(raster_name)!="RasterLayer"){
    layer_rast<-raster(raster_name)
  }else{
    layer_rast <- raster_name
    raster_name <- filename(layer_rast)
  }
  
  ## Create output raster name if out_rast_name is null
  if(is.null(out_rast_name)){
    extension_str <- extension(raster_name)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
    if(out_suffix!=""){
      out_rast_name <- file.path(out_dir,paste(raster_name_tmp,"_crop_proj_reg_",out_suffix,file_format,sep="")) #for use in function later...
    }else{
      out_rast_name <- file.path(out_dir,paste(raster_name_tmp,"_crop_proj_reg",out_suffix,file_format,sep="")) #for use in function later...
    }
  }
  
  ## Get the input raster projection information if needed
  if(is.null(input_proj_str)){
    input_proj_str <-projection(layer_rast)   #Extract current coordinates reference system in PROJ4 format
  }else{
    projection(layer_rast) <- input_proj_str #assign projection info
  }
  region_temp_projected <- projectExtent(reg_ref_rast,CRS(input_proj_str))     #Project from ref to current region coord. system
  
  layer_crop_rast <- crop(layer_rast, region_temp_projected) #crop using the extent from the region tile
  #layer_projected_rast<-projectRaster(from=layer_crop_rast,crs=proj4string(reg_outline),method="ngb")
  if(agg_param[1]==TRUE){
    agg_fact <- as.numeric(agg_param[2]) #in case we have a string/char type
    agg_fun <- agg_param[3]
    #debug(aggregate_raster)
    r_agg_raster_name <- aggregate_raster(reg_ref_rast, #reference raster with the desired resolution
                                          agg_fact=agg_fact, #given aggregation factor
                                          r_in=layer_crop_rast, #raster to be aggregated
                                          agg_fun="mean", #aggregation function
                                          out_suffix=out_suffix,
                                          file_format=".tif",
                                          out_dir=out_dir)
    layer_crop_rast <- raster(r_agg_raster_name)
  }
  #Should check if different projection!!!
  layer_projected_rast <- projectRaster(from=layer_crop_rast,to=reg_ref_rast,method="ngb",
                                        filename=out_rast_name,overwrite=TRUE)
  
  #Need cleanup of tmp files here!!! building up to 19gb!
  removeTmpFiles(h=0)
  
  return(out_rast_name)
}

download_and_process_prism_data <- function(in_dir,out_dir,start_date,end_date,var_name,ref_rast_name,
                                            agg_param,num_cores=1,create_out_dir_param=FALSE,NA_value=-9999,out_suffix="",
                                            file_format=".tif", download_file=T,unzip_files=T,match_file=T,lf_zip=NULL,lf_r=NULL){
  
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
  
  ######### PART 1: downloading prsim product ######
  if(download_file==T){
    setwd(in_dir)
    #debug(downloading_prism_product)
    download_obj <- downloading_prism_product(start_date, end_date,var_name,num_cores,in_dir=in_dir,out_dir=NULL)
    setwd(out_dir)
  }else{
    download_obj<-try(load_obj(file.path(in_dir,paste("downloaded_prism_data_",var_name,".RData",sep=""))))
  }
  
  ####### PART 2: extracting files #######
  
  ## run by year!!
  ## loop through year...or make this a function?
  if(unzip_files==T){
    if(is.null(lf_zip)){
      lf_zip <- download_obj$lf_zip
    }
    nb_year <- length(lf_zip)
    list_lf_r <- vector("list",length=nb_year)
    for(i in 1:nb_year){
      out_dir_year <- unique(dirname(lf_zip[[i]]))
      lf_r <- lapply(lf_zip[[i]], unzip,exdir= out_dir_year)
      lf_r <- list.files(pattern="*bil.bil$",path=out_dir_year,full.names = T)
      list_lf_r[[i]] <- lf_r
    }
  }else{
    #find where the data is
    #
    #start_date <- "2012-01-01" # param 3
    #end_date <- "2012-12-31" # param 4
    start_date <- as.Date(start_date,format="%Y-%m-%d") #start date
    end_date <- as.Date(end_date,format="%Y-%m-%d") #end date
    dates_range <- seq.Date(start_date, end_date, by="1 day") #sequence of dates
    #dates_range_format <- as.Date(dates_range,format="%Y%m%d") #end date
    
    ##Check if the range is over multiple years
    date_year <- strftime(dates_range, "%Y")
    date_month <- strftime(dates_range , "%m") # current month of the date being processed
    date_day <- strftime(dates_range , "%d")
    #lf_r <- list.files(pattern="*bil.bil$",path=out_dir_year,full.names = T)
    nb_year <- length(unique(date_year))
    list_year <- unique(date_year)
    list_lf_r <- vector("list",length=nb_year)
    for(i in 1:nb_year){
      
      out_dir_year <- file.path(in_dir,paste("prism_",var_name,"_",list_year[i],sep=""))
      lf_r <- list.files(pattern="*bil.bil$",path=out_dir_year,full.names = T)
      #lf_r <- lapply(lf_zip[[i]], unzip,exdir= out_dir_year)
      #lf_r <- list.files(pattern="*bil.bil$",path=out_dir_year,full.names = T)
      list_lf_r[[i]] <- lf_r
    }
    
  }
  
  ########## PART 3: Match to study area ##############
  
  ## Match to the study area...
  ## Now crop and reproject if necessary
  #Use function above
  r_ref <- raster(ref_rast_name)
  list_lf_r_reg <- vector("list",length(list_lf_r))
  for(i in 1:length(list_lf_r)){
    lf_r <- list_lf_r[[i]]
    out_dir_year <- unique(dirname(lf_r))
    #plot(r_ref)
    #agg_param <- c(FALSE,NULL,"mean") #False means there is no aggregation!!!
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
  
  ##########
  return(list_lf_r_reg)
}

############## End of script ################


