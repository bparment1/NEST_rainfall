####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script contains functions to process and explore the correlation between rainfall events 
#and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/09/2015 
#DATE MODIFIED: 12/11/2015
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

#This function creates a spatial polygon data frame object for the extent matching a raster input
create_polygon_from_extent<-function(reg_ref_rast,outDir=NULL,outSuffix=NULL){
  #This functions returns polygon sp from input rast
  #Input Arguments: 
  #reg_ref_rast: input ref rast or spatial object (polygons, points df)
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
  
  ref_e <- extent(reg_ref_rast) #extract extent from raster object/SpatialPointsDataframe/SpatialPolygonsDataframe
  reg_outline_poly <- as(ref_e, "SpatialPolygons") #coerce raster extent object to SpatialPolygons from sp package 
  reg_outline_poly <- as(reg_outline_poly, "SpatialPolygonsDataFrame") #promote to spdf
  proj4string(reg_outline_poly) <- projection(reg_ref_rast) #Assign projection to spdf
  infile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="") #name of newly crated shapefile with the extent
  writeOGR(reg_outline_poly,dsn= outDir,layer= sub(".shp","",infile_reg_outline), 
           driver="ESRI Shapefile",overwrite_layer="TRUE")
  
  return(reg_outline_poly) #return spdf
}

plotting_coliform_and_rainfall <- function(i,df_ts_pix,data_var,list_selected_ID,plot_fig=T){
  
  # Input arguments:
  # i : selected station
  # df_ts_pix_data : data extracted from raster layer
  # data_var : data with coliform measurement
  # list_selected_ID : list of selected station
  # plot_fig : if T, figures are plotted11111
  # Output
  #
  
  ##### START FUNCTION ############
  
  #get the relevant station
  id_name <- list_selected_ID[i] # e.g. WS037.00
  id_selected <- df_ts_pix[[var_ID]]==id_name
  
  ### Not get the data from the time series
  data_pixel <- df_ts_pix[id_selected,]
  var_pix_ts <- t(as.data.frame(subset(data_pixel,select=var_name)))
  #pix <- t(data_pixel[1,24:388])#can subset to range later
  pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
  
  ## Process the coliform data
  
  #there are several measurements per day for some stations !!!
  #id_name <- data_pixel[[var_ID]]
  
  df_tmp  <-data_subset[data_subset$LOCATION_ID==id_name,]
  #aggregate(df_tmp
  var_pix <- aggregate(COL_SCORE ~ TRIP_START_DATE_f, data = df_tmp, mean) #aggregate by date
  #length(unique(test$TRIP_START_DATE_f))
  
  #var_pix <- subset(as.data.frame(data_subset[id_selected,c(var_name,"TRIP_START_DATE_f")])) #,select=var_name)
  
  d_z <- zoo(pix_ts,idx) #make a time series ...
  names(d_z)<- "rainfall"
  d_var <- zoo(var_pix,var_pix$TRIP_START_DATE_f)
  #plot(d_var,pch=10)
  
  d_z2 <- merge(d_z,d_var)
  d_z2$TRIP_START_DATE_f <- NULL
  
  df2 <- as.data.frame(d_z2)
  df2$date <- rownames(df2)
  rownames(df2) <- NULL
  df2$COL_SCORE <- as.numeric(as.character(df2$COL_SCORE))
  df2$rainfall <- as.numeric(as.character(df2$rainfall))
  df2$LOCATION_ID <- id_name
    
  #plot(df2$rainfall)
  #list_pix[[i]] <- pix_ts
  
  if(plot_fig==T){
    
    res_pix <- 480
    col_mfrow <- 2
    row_mfrow <- 1
    
    ###
    #Figure 3b
    png(filename=paste("Figure3b_","pixel_profile_var_combined_",id_name,"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    #plot(d_z,lty=2,ylab="rainfall",xlab="Time",main="")
    #points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
    plot(d_z,lty=2,ylab="rainfall",xlab="Time",main="")
    abline(h=threshold_val,col="green")
    
    par(new=TRUE)              # key: ask for new plot without erasing old
    #plot(x,y,type="l",col=t_col[k],xlab="",ylab="",lty="dotted",axes=F) #plotting fusion profile
    plot(df2$COL_SCORE,pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
    #points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
    legend("topleft",legend=c("stations"), 
           cex=1.2,col="red",pch =10,bty="n")
    
    axis(4,cex=1.2)
    mtext(4, text = "coliform scores", line = 3)
    
    title(paste("Station time series",id_name,sep=" "))
    
    dev.off()
    
    #Figure 3c
    png(filename=paste("Figure3c_","pixel_profile_var_combined_log_scale_",id_name,"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    #plot(d_z,lty=2,ylab="rainfall",xlab="Time",main="")
    #points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
    plot(d_z,lty=2,ylab="rainfall",xlab="Time",main="")
    abline(h=threshold_val,col="green")
    par(new=TRUE)              # key: ask for new plot without erasing old
    #plot(x,y,type="l",col=t_col[k],xlab="",ylab="",lty="dotted",axes=F) #plotting fusion profile
    plot(log(df2$COL_SCORE),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
    
    #points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
    legend("topleft",legend=c("stations"), 
           cex=1.2,col="red",pch =10,bty="n")
    
    axis(4,cex=1.2)
    mtext(4, text = "coliform scores", line = 3)
    
    title(paste("Station time series",id_name,sep=" "))
    
    dev.off()
    
    ####Histogram of values
    
    res_pix <- 480
    col_mfrow <- 2
    row_mfrow <- 1
    
    png(filename=paste("Figure4_","histogram_coliform_measurements_",year_processed,"_",id_name,"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    hist_val <- hist(df2$COL_SCORE,main="",xlab="COLIFORM SCORES")
    title(paste("Histrogram of coliform scores for station",id_name,"in",year_processed,sep=" "))
    #abline(v=threshold_val,col="green" )
    legend("topright",legend=c("treshold val"), 
           cex=1.2, col="green",lty =1,bty="n")  
    
    y_loc <- max(hist_val$counts)/2
    
    #text(threshold_val,y_loc,paste(as.character(threshold_val)),pos=1,offset=0.1)
    
    dev.off()
    
    #res_pix <- 480
    #col_mfrow <- 2
    #row_mfrow <- 1
    
    #png(filename=paste("Figure4_","histogram_coliform_measurements_",year_processed,"_",id_name,"_",out_suffix,".png",sep=""),
    #    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    plot(df2$rainfall)
    plot(df2$rainfall,df2$COL_SCORE)
    plot(df2$rainfall,log(df2$COL_SCORE))
    plot(log(df2$rainfall),log(df2$COL_SCORE))
    
    
  }
  
  ## Now correlation.
  #sum(is.na(df2$rainfall))
  #[1] 0
  nb_zero <- sum((df2$rainfall==0)) #203
  nb_NA <- sum(is.na(df2$COL_SCORE))
  ## Cumulated precip and lag?
  #Keep number of  0 for every year for rainfall
  #summarize by month
  #Kepp number of NA for scores... 
  #Summarize by season...
  ## Threshold?
  station_summary_obj <- list(nb_zero,nb_NA,df2)
  names(station_summary_obj) <- c("nb_zero","nb_NA","df_combined")
  return(station_summary_obj)
}


########################### End of script #####################################