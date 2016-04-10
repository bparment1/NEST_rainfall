####################################  NEST Beach closure project  #######################################
###########################################  Figures production  #######################################
#This script contains functions to process and explore the correlation between rainfall events 
#and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/09/2015 
#DATE MODIFIED: 04/10/2016
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

#function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_functions.R" #PARAM 1
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
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

plotting_measurements_and_rainfall <- function(i,dates_val,df_ts_pix,data_var,list_selected_ID,r_ts_name,var_name,dates_str,plot_fig=T){
  
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
  data_pixel <- as.data.frame(data_pixel)
  
 
  pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
  #pix_ts <- subset(as.data.frame(pix_ts),select=r_ts_name)
  pix_ts <- (as.data.frame(pix_ts))
  ## Process the coliform data
  
  #there are several measurements per day for some stations !!!
  #id_name <- data_pixel[[var_ID]]
  
  #df_tmp  <-data_var[data_var$LOCATION_ID==id_name,]
  df_tmp <- subset(data_var,data_var$LOCATION_ID==id_name)
  #aggregate(df_tmp
  if(nrow(df_tmp)>1){
    
    formula_str <- paste(var_name," ~ ","TRIP_START_DATE_f",sep="")
    #var_pix <- aggregate(COL_SCORE ~ TRIP_START_DATE_f, data = df_tmp, mean) #aggregate by date
    var_pix <- aggregate(as.formula(formula_str), data = df_tmp, FUN=mean) #aggregate by date
    #length(unique(test$TRIP_START_DATE_f))
    #var_pix_ts <- t(as.data.frame(subset(data_pixel,select=var_name)))
    #pix <- t(data_pixel[1,24:388])#can subset to range later
  }
  
  #var_pix <- subset(as.data.frame(data_subset[id_selected,c(var_name,"TRIP_START_DATE_f")])) #,select=var_name)
  
  #Create time series object from extract pixel time series
  d_z <- zoo(pix_ts,dates_val) #make a time series ...
  names(d_z)<- "rainfall"
  #Create date object for data from stations
  d_var <- zoo(var_pix,var_pix$TRIP_START_DATE_f)
  #plot(d_var,pch=10)
  
  d_z2 <- merge(d_z,d_var)
  ##Now subset?
  d_z2 <- window(d_z2,start=dates_val[1],end=dates_val[length(dates_val)])
  
  d_z2$TRIP_START_DATE_f <- NULL
  
  df2 <- as.data.frame(d_z2)
  df2$date <- rownames(df2)
  rownames(df2) <- NULL
  df2[[var_name]] <- as.numeric(as.character(df2[[var_name]]))
  
  #df2$COL_SCORE <- as.numeric(as.character(df2$COL_SCORE))
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
    plot(df2[[var_name]],pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
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
    #plot(log(df2$COL_SCORE),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
    plot(log(df2[[var_name]]),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
    
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
    
    hist_val <- hist(df2[[var_name]],main="",xlab="COLIFORM SCORES")
    #hist_val <- hist(df2$COL_SCORE,main="",xlab="COLIFORM SCORES")
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
    #plot(df2$rainfall,df2$COL_SCORE)
    #plot(log(df2$rainfall),log(df2$COL_SCORE))
    plot(df2$rainfall,df2[[var_name]])
    plot(df2$rainfall,log(df2[[var_name]]))

    
    
  }
  
  ## Now correlation.
  #sum(is.na(df2$rainfall))
  #[1] 0
  nb_zero <- sum((df2$rainfall==0)) #203
  #nb_NA <- sum(is.na(df2$COL_SCORE))
  nb_NA <- sum(is.na(df2[[var_name]]))
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

combine_stations_data_raster_ts_fun <- function(data,dat_stat,convert_to_inches,in_dir_rst,start_date,end_date,data_type,coord_names,out_dir,out_suffix){
  #data,convert_to_inches,in_dir_rst,start_date,end_date,data_type="MH",coord_names,out_dir,out_suffix
  #data,dat_stat,convert_to_inches,in_dir_rst,start_date,end_date,data_type,coord_names,out_dir,out_suffix  
  #
  ###Add documentation here...
  #data
  #dat_stat: unique station data
  #convert_to_inches
  #in_dir_rst
  #start_date
  #end_date
  #data
  #year_processed <- "2012" #PARAM 16
  #threshold_val <- 2*25.4 #PARAM 17, in inches or mm
  #units_val <- "mm"
  #out_dir
  #out_suffix
  
  #### Start script ###
  
  if(data_type=="MHB"){
    dates_TRIP_START <- unlist(lapply(strsplit(data$SAMPLE.DATE," "),function(x){x[1][1]}))
    #dates_TRIP_START <- gsub(" 0:00:00","",data$SAMPLE.DATE)
    data$TRIP_START_DATE_f <- as.Date(strptime(dates_TRIP_START,"%m/%d/%Y"))
    data$TRIP_START_DATE_year <- four.digit.year(data$TRIP_START_DATE_f , year=1968)
    #format(data$TRIP_START_DATE_f , format="%m-%d-%Y")
    data$TRIP_START_DATE_month <- strftime(data$TRIP_START_DATE_f , "%m") # current month of the date being processed
    data$TRIP_START_DATE_day <- strftime(data$TRIP_START_DATE_f , "%d")
  }
  if(data_type=="DMR"){
    dates_TRIP_START <- gsub(" 0:00:00","",data$TRIP_START_DATE)
    data$TRIP_START_DATE_f <- as.Date(strptime(dates_TRIP_START,"%m/%d/%Y"))
    data$TRIP_START_DATE_month <- strftime(data$TRIP_START_DATE_f , "%m") # current month of the date being processed
    data$TRIP_START_DATE_year <- strftime(data$TRIP_START_DATE_f , "%Y")
    data$TRIP_START_DATE_day <- strftime(data$TRIP_START_DATE_f , "%d")
  }
  
  data$TRIP_START_DATE_f <- paste0(data$TRIP_START_DATE_year,data$TRIP_START_DATE_month,data$TRIP_START_DATE_day)
  data$TRIP_START_DATE_f <- as.Date(strptime(data$TRIP_START_DATE_f,"%Y%m%d"))
  

  r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=in_dir_rst,full.names=T))) #rainfall time series stack
  
  if (convert_to_inches==TRUE){
    r_rainfall <- r_rainfall/25.4 #improve efficiency later? YES!!
  }
  
  #dat_stat$LOCATION_ID <- as.character(dat_stat$LOCATION_ID)
  #nrow(dat_stat)==length(unique(dat_stat$LOCATION_ID)) #Checking that we have a unique identifier for each station
  
  ## Plot mosaics for Maine for daily predictions in 2014
  ## Get pixel time series at centroids of tiles used in the predictions
  
  df_ts_pixel <- extract(r_rainfall,dat_stat,df=T,sp=T)
  #test_tmp <- merge(df_ts_pixel,data_subset,by="LOCATION_ID")
  #df_ts_pixel <- test_tmp
  #df_ts_pixel <- cbind(summary_metrics_v,df_ts_pixel)
  r_ts_name <- names(r_rainfall)
  #d_z <- zoo(df_ts_pixel,idx) #make a time series .
  
  #freq_station <- sort(table(data_subset$LOCATION_ID),decreasing=T) # select top 2 stations in term of availability
  #list_selected_ID <- names(freq_station)[1:25] #select top 25
  #list_selected_ID <- names(freq_station) #select top 25
  
  #View(freq_station)
  
  ##This will be a function later on...
  df_ts_pix <- df_ts_pixel#this contains the pixels with extracted pixels
  #list_selected_pix <- 11:14
  year_processed <- year(start_date)
  write.table(df_ts_pix,paste("df_ts_pix_",year_processed,".txt",sep=""))
  #test <- merge(dat_stat,data,by="FID")
  #data_merged <- merge(data,dat_stat,by="FID",all=T,suffixes=c("_x","_y"))
  #data_merged <- merge(data,dat_stat,by="i",all=T,suffixes=c("","_y"))
  
  data_merged <- merge(data,dat_stat,by="ID_stat",all=T,suffixes=c("","_y"))
  #data_merged <- data
  #### NOW SELECT RELEVANT DATES
  #need to figure out the join!!
  #paste(coord_names_tmp,"_y"), drop all the y
  #data$id_test <- paste(data[[coord_names[1]]],data[[coord_names[2]]],sep="_")
  #coord_names <- c("SITE.LONGITUDE..UTM.","SITE.LATITUDE..UTM.")
  dates_val <- seq(as.Date(start_date), as.Date(end_date), 'day')
  #date_l <- strptime(idx[1], "%Y%m%d") # 
  dates_str <- format(dates_val, "%Y%m%d") #  date being processed
  
  ###Get the relevant dates from the original dataset
  #data_subset <- data_merged[data_merged$TRIP_START_DATE_f >= as.Date(start_date) & data_merged$TRIP_START_DATE_f <= as.Date(end_date), ]
  #data_subset$LOCATION_ID <- as.character(data_subset$LOCATION_ID)
  #data_subset[[coord_names[1]]] <- as.numeric(data_subset[[coord_names[1]]])
  #data_subset[[coord_names[2]]] <- as.numeric(data_subset[[coord_names[2]]])
  
  #list_selected_ID <- unique(data_subset$LOCATION_ID)
  list_selected_ID <- unique(data_merged$LOCATION_ID)
  
  #list_selected_ID <- df_ts_pix$LOCATION_ID
  
  list_pix <- vector("list",length=length(list_selected_ID))
  
  #debug(plotting_measurements_and_rainfall)
  #i <- 1
  
  #test <- lapply(1:1,FUN=plotting_measurements_and_rainfall,
  #               df_ts_pix=df_ts_pix,data_var=data_merged,list_selected_ID=list_selected_ID,var_name=var_name,r_ts_name=r_ts_name,dates_val,plot_fig=F)
  
  #Takes 5mintues or less on bpy50 laptop
  num_cores <- 4
  list_df_combined <- mclapply(1:length(list_selected_ID),FUN=plotting_measurements_and_rainfall,
                               df_ts_pix=df_ts_pix,data_var=data_merged,list_selected_ID=list_selected_ID,var_name=var_name,r_ts_name=r_ts_name,dates_val,plot_fig=F,
                               mc.preschedule=FALSE,mc.cores= num_cores)
  #list_df_combined <- mclapply(1:4,FUN=plotting_measurements_and_rainfall,
  #                             df_ts_pix=df_ts_pix,data_var=data_merged,list_selected_ID=list_selected_ID,var_name=var_name,r_ts_name=r_ts_name,dates_val,plot_fig=F,
  #                             mc.preschedule=FALSE,mc.cores= num_cores)
  
  save(list_df_combined,file= file.path(out_dir,paste("list_df_combined_obj",year_processed,"_",
                                                      out_suffix,".RData",sep="")))
  
  #l_png_files <- mclapply(1:length(unlist(lf_mean_mosaic)),FUN=plot_mosaic,
  #                        list_param= list_param_plot_mosaic,
  #                        mc.preschedule=FALSE,mc.cores = num_cores)
  
  list_cleaning_df <- lapply(1:length(list_df_combined),FUN=function(i,x){x[[i]]$df_combined},x=list_df_combined)
  data_df <- do.call(rbind,list_cleaning_df)
  write.table(data_df,file= file.path(out_dir,paste("data_df","_",year_processed,"_",
                                                      out_suffix,".txt",sep="")),sep=",")
  return(data_df)
}


########################### End of script #####################################