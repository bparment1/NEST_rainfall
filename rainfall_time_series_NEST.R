##############################################  NEST Beach closure project  #######################################
###########################################  Figures and Analyses production  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/05/2015 
#DATE MODIFIED: 12/11/2015
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - The script does not use bacteria data at the current time.  
#          - Add spacetime object functions for later
#TO DO:
# - Select 2 inches rainfall events and correlates with bacteria data
# - Compute accumulated rain over several days using time series functions
# - Make a movie sequence later on using animation package in R
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

###### Functions used in this script sourced from other files

function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_12112015.R" #PARAM 1
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

CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 #param 2
CRS_reg <- CRS_WGS84 # PARAM 3

file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_flag_val <- NA_value #PARAM6
out_suffix <-"NEST_prism_12082015" #output suffix for the files and ouptu folder #PARAM 7
create_out_dir_param=TRUE #PARAM8
num_cores <- 11 #PARAM 9

rainfall_dir <- "prism_rain/prismrain2012" #PARAM 10
station_data_fname <- file.path(in_dir, "WQ_TECS_Q.txt") #PARAM 11

start_date <- "2012-01-01" #PARAM 12
end_date <- "2012-12-31" #PARAM 13
var_name <- "COL_SCORE" #PARAM 14
var_ID <- "LOCATION_ID" #PARAM 15
year_processed <- "2012" #PARAM 16
threshold_val <- 2*25.4 #PARAM 17, in inches or mm
convert_to_inches <- FALSE #PARAM 18
units_val <- "mm"

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

r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=file.path(in_dir,rainfall_dir),full.names=T))) #rainfall time series stack

#### Read in data and add time flag

plot(r_rainfall,y=1)
data <- read.table(station_data_fname,sep=",",header=T) #this is T-mode using cor matrix

dates_TRIP_START <- gsub(" 0:00:00","",data$TRIP_START_DATE)
data$TRIP_START_DATE_f <- as.Date(strptime(dates_TRIP_START,"%m/%d/%Y"))
data$TRIP_START_DATE_month <- strftime(data$TRIP_START_DATE_f , "%m") # current month of the date being processed
data$TRIP_START_DATE_year <- strftime(data$TRIP_START_DATE_f , "%Y")
data$TRIP_START_DATE_day <- strftime(data$TRIP_START_DATE_f , "%d")
dim(data)
class(data$LOCATION_ID)
data$LOCATION_ID <- as.character(data$LOCATION_ID)
length(unique(data$LOCATION_ID)) #There are 2851 stations

if (convert_to_inches==TRUE){
  r_rainfall <- r_rainfall/25.4 #improve efficiency later? YES!!
}


#### NOW SELECT RELEVANT DATES

idx <- seq(as.Date(start_date), as.Date(end_date), 'day')
#date_l <- strptime(idx[1], "%Y%m%d") # 
dates_l <- format(idx, "%Y%m%d") #  date being processed

data_subset <- data[data$TRIP_START_DATE_f >= as.Date(start_date) & data$TRIP_START_DATE_f <= as.Date(end_date), ]
data_subset$LOCATION_ID <- as.character(data_subset$LOCATION_ID)

#Remote stations without coordinates and make a SPDF
dat_stat <- subset(data_subset, !is.na(LONGITUDE_DECIMAL) & !is.na(LATITUDE_DECIMAL))
coords <- dat_stat[,c('LONGITUDE_DECIMAL','LATITUDE_DECIMAL')]
coordinates(dat_stat) <- coords
proj4string(dat_stat) <- projection(r_rainfall) #this is the NAD83 latitude-longitude

## Remove duplicates rows from stations to identify uniques sations
dat_stat <- remove.duplicates(dat_stat)
#dat_stat$LOCATION_ID <- as.character(dat_stat$LOCATION_ID)
nrow(dat_stat)==length(unique(dat_stat$LOCATION_ID)) #Checking that we have a unique identifier for each station

r_rainfall <- setZ(r_rainfall, idx) #for now, this can also be made into a spacetime object

x <- zApply(r_rainfall, by="day",fun=mean,name="overall mean") #overall mean, takes about a minute
raster_name <- paste("day","_","overall_mean",file_format,sep="")
writeRaster(x, file=raster_name,overwrite=T)
plot_to_file(raster_name) #quick plot of raster to disk

#x <- zApply(r_rainfall, by=c(1,24),fun=mean,name="overall mean") #overall mean
r_date <- getZ(r_rainfall)

## Plot mosaics for Maine for daily predictions in 2014
## Get pixel time series at centroids of tiles used in the predictions

df_ts_pixel <- extract(r_rainfall,dat_stat,df=T,sp=T)
test<-merge(df_ts_pixel,data_subset,by="LOCATION_ID")

#df_ts_pixel <- cbind(summary_metrics_v,df_ts_pixel)
r_ts_name <- names(r_rainfall)
#d_z <- zoo(df_ts_pixel,idx) #make a time series .

res_pix <- 480
col_mfrow <- 1
row_mfrow <- 1

png(filename=paste("Figure1_","histogram","_station_coliform_measurements_frequency_",year_processed,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(table(data_subset$LOCATION_ID),type="h", main="Number of measurements",
     ylab="Frequency of coliform measurements",xlab="Station ID")

dev.off()

#hist(table(data_subset$COL_SCORE))
#Could output in a textfile
range(table(dat_stat$LOCATION_ID)) #1 to 60
mean(table(dat_stat$LOCATION_ID)) #8
median(table(data_subset$LOCATION_ID)) #8
hist(table(data_subset$LOCATION_ID), main="Frequency of Number of measurements by stations")

ref_pol <- create_polygon_from_extent(dat_stat,outDir=out_dir,outSuffix=out_suffix)
  
res_pix <- 480
col_mfrow <- 1.5
row_mfrow <- 1

png(filename=paste("Figure2a_","rainfall_map_",dates_l[1],"_and_stations_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_rainfall,y=1)
plot(dat_stat,add=T)
text(dat_stat,dat_stat$tID,cex=1.4)
plot(ref_pol,border="red",add=T)
legend("topright",legend=c("stations"), 
       cex=1.2, col="black",pch =3,bty="n")

dev.off()

col_mfrow <- 2
row_mfrow <- 1

png(filename=paste("Figure2b_","rainfall_map_",dates_l[1],"_and_stations_zoom_window_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_rainfall,y=1,ext=extent(dat_stat))
plot(dat_stat,add=T,pch=3)
text(dat_stat,dat_stat$tID,cex=1.4)
legend("topright",legend=c("stations"), 
       cex=1.2, col="black",pch =3,bty="n")

dev.off()

#
#make a function later on?

#ADD selection function by ID
#list_selected_pix: by LOCATION ID
#

freq_station <- sort(table(data_subset$LOCATION_ID),decreasing=T) # select top 2 stations in term of availability
list_selected_ID <- names(freq_station)[1:25] #select top 25
list_selected_ID <- names(freq_station) #select top 25

#View(freq_station)

##This will be a function later on...
df_ts_pix <- df_ts_pixel#this contains the pixels with extracted pixels
#list_selected_pix <- 11:14
list_pix <- vector("list",length=length(list_selected_ID))

#debug(plotting_coliform_and_rainfall)
#i <- 1

#test <- lapply(1:length(list_selected_ID),FUN=plotting_coliform_and_rainfall,
#               df_ts_pix=df_ts_pix,data_var=data_var,list_selected_ID=list_selected_ID,plot_fig=T)

#Takes 5mintues or less on bpy50 laptop
num_cores <- 4
list_df_combined <- mclapply(1:length(list_selected_ID),FUN=plotting_coliform_and_rainfall,
                    df_ts_pix=df_ts_pix,data_var=data_var,list_selected_ID=list_selected_ID,plot_fig=T,
                    mc.preschedule=FALSE,mc.cores= num_cores)

save(list_df_combined,file= file.path(out_dir,paste("list_df_combined_obj",out_suffix,".RData",sep="")))

#l_png_files <- mclapply(1:length(unlist(lf_mean_mosaic)),FUN=plot_mosaic,
#                        list_param= list_param_plot_mosaic,
#                        mc.preschedule=FALSE,mc.cores = num_cores)

list_cleaning_df <- lapply(1:length(list_df_combined),FUN=function(i,x){x[[i]]$df_combined},x=list_df_combined)
data_df <- do.call(rbind,list_cleaning_df)
#colnames(data_df) <- list_selected_ID
#data_dz <- zoo(data_dz,idx)

plot(log(data_df$rainfall),log(data_df$COL_SCORE))
plot(data_df$rainfall,data_df$COL_SCORE)

## IDENTIFY 2 inches events?


############################ END OF SCRIPT #######################
