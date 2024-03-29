##############################################  NEST Beach closure project  #######################################
#################################  Analyes for exploration of station measurements  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#It uses time series processed earlier in R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 06/14/2016 
#DATE MODIFIED: 07/09/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - Adding figure option and title of station in plot  

#
#TO DO:

#################################################################################################

###Loading R library and packages                                                      

library(raster)                 # loading the raster package
library(gtools)                 # loading R helper programming tools/functions
library(sp)                     # spatial objects in R
#library(gplots)                 # plotting functions such as plotCI
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
#library(hydrostats)
library(shiny)
library(plyr)                                # Various tools including rbind.fill

###### Functions used in this script sourced from other files

#function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_03272016.R" #PARAM 1
function_analyses_NEST <- "analyses_NEST_functions_07072016.R" #PARAM 1

#script_path <- "." #path to script #PARAM 
script_path <- "/home/parmentier/Data/NEST/R_NEST"
#script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#script_path <- "/home/bparmentier/Google Drive/NEST/NEST_stations_s06/" #path to script #PARAM 
#script_path <- "/home/benoit/data/NEST_stations_s06" #on SSI server
#setwd(script_path)

#script_path <- "/home/parmentier/Data/rainfall/NEST"
#source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses_NEST)) #source all functions used in this script 1.

##### Functions used in this script 


#####  Parameters and argument set up ###########

#in_dir <- "/home/bparmentier/Google Drive/NEST/NEST_stations_s06" #local bpy50 , param 1
#in_dir <- "./NEST_stations_s05" #NCEAS, param 
#in_dir <- "/home/benoit/data/NEST_stations_s06" #U. Maine SSI server, param 
#in_dir <- "." #Use current directory, run locally
#in_dir <- "/home/parmentier/Data/NEST/NEST_stations_s08"
in_dir <- "/home/parmentier/Data/NEST/data_analyses_NEST"
#in_dir <- "/home/bparmentier/Google Drive/NEST/NEST_stations_s08"

CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 #param 2
CRS_reg <- CRS_WGS84 # PARAM 3

file_format <- ".rst" #PARAM 4
NA_value <- -9999 #PARAM5
NA_flag_val <- NA_value #PARAM6
#out_suffix <-"NEST_prism_03282016" #output suffix for the files and ouptu folder #PARAM 7
#out_suffix <- paste0(as.Date(Sys.Date(),format="%Y%m%d"))
out_suffix <- paste0("NEST_analyses_",format(Sys.Date(),"%Y%m%d"))

create_out_dir_param=T #PARAM8
num_cores <- 11 #PARAM 9

#rainfall_dir <- "/home/bparmentier/Google Drive/NEST_Data" #PARAM 10
rainfall_dir <- "./data" #PARAM 10
#reg_ref_rast_name <- "./data/PRISM_ppt_stable_4kmD2_20111222_crop_proj_reg.tif"
reg_ref_rast_name <- "PRISM_ppt_stable_4kmD2_20111222_crop_proj_reg.tif"
#station_data_fname <- file.path("/home/bparmentier/Google Drive/NEST_Data/", "WQ_TECS_Q.txt") #PARAM 11
#station_data_fname <- file.path("data", "MHB_data_2006-2015.csv") #PARAM 11


#station_measurements_MHB_data_fname <- file.path("data", "data_df_rainfall_and_measurements_MHB.txt") #PARAM 11 
#This will change
#station_measurements_DMR_data_fname <- file.path("data", "data_df_rainfall_and_measurements_DMR.txt") #PARAM 11 
#station_measurements_DMR_data_fname <- list.files(path=file.path(in_dir,"data"),
#                                                  pattern="data_df_combined_.*._DMR.txt",
#                                                  full.names=T) #PARAM 11 
#station_measurements_MHB_data_fname <- list.files(path=file.path(in_dir,"data"),
#                                                  pattern="data_df_.*._NEST_prism_MHB.txt",
#                                                  full.names=T) #PARAM 11 

station_measurements_DMR_data_fname <- list.files(path=in_dir,
                                                  pattern="data_df_combined_.*._DMR.txt",
                                                  full.names=T) #PARAM 11 
station_measurements_MHB_data_fname <- list.files(in_dir,
                                                  pattern="data_df_.*._NEST_prism_MHB.txt",
                                                  full.names=T) #PARAM 11 
#data_df_2015_NEST_prism_MHB.txt
start_date <- "2012-01-01" #PARAM 12,  this is the default value, use user define otherwise
end_date <- "2012-12-31" #PARAM 13,  this is the default rer define otherwise
var_name_DMR <- "COL_SCORE" #PARAM 14
var_name_MHB <- "CONCENTRATION" #PARAM 14, MHB data, need to add DMR
#var_name <- var_name_MHB #default for shiny app
var_name <- var_name_DMR #default for anaylses

##Raster start date
start_date_r <- start_date #set the same start date for now
end_date_r <- end_date #set the same start date for now

var_ID <- "LOCATION_ID" #PARAM 15
year_processed <- "2012" #PARAM 16
threshold_val <- 2*25.4 #PARAM 17, in inches or mm
convert_to_inches <- FALSE #PARAM 18
units_val <- "mm"
#data_type <- "MHB" #for Maine beach health used as default, it is user defined otherwise
data_type <- "DMR" #for Maine beach health used as default, it is user defined otherwise

## Change coordinates to x and y and lat long!!!
coord_names_MHB <- c("SITE.LONGITUDE..UTM.","SITE.LATITUDE..UTM.") #MH beach bacteria dataset
coord_names_DMR <- c("LONGITUDE_DECIMAL","LATITUDE_DECIMAL") #cloroforms beach bacteria dataset
#data_df_fname <- "./data/df_ts_pix_2012.txt"
#/home/bparmentier/Google Drive/NEST/NEST_stations_s02/data

SMAZones_fname <- "SMAZoneDissolve.shp"

#dat_stat_location_DMR_fname <- "dat_stat_location_DMR.shp" #Use this for now
dat_stat_location_DMR_fname <- "dat_stat_location_DMR.txt" 
#dat_stat_location_DMR.shp
#dat_stat_location_MHB_fname <- "dat_stat_location_MHB.shp" #Use this for now
dat_stat_location_MHB_fname <- "dat_MHB_location_DMR.txt" 

transform_fun <- "log1p"
plot_fig <- F #if false no figures are created

################# START SCRIPT ###############################

#set up the working directory
#Create output directory

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
setwd(in_dir)
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}


###############

#data_df <- data_df_MHB #default dataset!!
#data_df_DMR
#dat_stat_location_MHB <- readOGR(file.path(in_dir,"/data"),sub(".shp","",dat_stat_location_MHB_fname))
#dat_stat_location_DMR <- readOGR(file.path(in_dir,"/data"),sub(".shp","",dat_stat_location_DMR_fname))
#dat_stat_location_MHB <- readOGR(in_dir,sub(".shp","",dat_stat_location_MHB_fname))
#dat_stat_location_DMR <- readOGR(in_dir,sub(".shp","",dat_stat_location_DMR_fname))

#dat_stat_location_MHB <- readOGR(in_dir,sub(".shp","",dat_stat_location_MHB_fname))
#dat_stat_location_DMR <- readOGR(in_dir,sub(".shp","",dat_stat_location_DMR_fname))

dat_stat_location_DMR <- read.table(file.path(in_dir,dat_stat_location_DMR_fname),sep=",")
coordinates(dat_stat_location_DMR) <- dat_stat_location_DMR[,coord_names_DMR]
proj4string(dat_stat_location_DMR) <- CRS_WGS84
#dat_stat_location_DMR <- readOGR(in_dir,sub(".shp","",dat_stat_location_DMR_fname))

list_ID_tmp <- unique(dat_stat_location_DMR$ID_stat)
list_ID_tmp <- unique(as.character(dat_stat_location_DMR[[var_ID]]))

tb <- read.table(station_measurements_DMR_data_fname[10],sep=",")
#list_ID <- unique(tb$ID_stat)
list_ID <- as.character(unique(tb[[var_ID]]))

selected_ID <- list_ID[1]
#selected_ID <- list_ID[5]  
## Should query the multiple !!!
#tb$ID_stat

selected_val <- selected_ID
#selected_col <- "ID_stat"
selected_col <- var_ID #ID_stat or LOCATION_ID

#undebug(read_select_station)
#df <- read_select_station(station_measurements_DMR_data_fname[10],selected_val,selected_col)

lf <- station_measurements_DMR_data_fname
#var_name <- var_name_MHB
var_name <- var_name_DMR
y_var_name <- var_name
x_var_name <- "rainfall"
#selected_ID <- 1
#
#undebug(run_lm_by_station)
#test <- run_lm_by_station(selected_ID,selected_col,x_var_name,y_var_name,lf,log_val=T,out_dir,out_suffix,num_cores)
out_suffix_str <- paste0(data_type,"_",out_suffix)
#started on 8:09

#debug(run_lm_by_station)
#test <- run_lm_by_station(list_ID[1],                     
#                    selected_col=selected_col,
#                     x_var_name=x_var_name,
#                     y_var_name=y_var_name,
#                     lf=lf,
#                      plot_fig=plot_fig,
#                    log_val=T,
#                     out_dir=out_dir,
#                     out_suffix=out_suffix,
#                     num_cores=1)

#transform_fun <- "log1p" #this will need to be switched for log_val
# run_lm_obj <- mclapply(list_ID[1:11],
#                      FUN=run_lm_by_station,
#                      selected_col=selected_col,
#                      x_var_name=x_var_name,
#                      y_var_name=y_var_name,
#                      lf=lf,
#                      plot_fig=plot_fig,
#                      log_val=T,
#                      out_dir=out_dir,
#                      out_suffix=out_suffix,
#                      num_cores=1,
#                      mc.preschedule=FALSE,
#                      mc.cores = num_cores)

run_lm_obj <- mclapply(list_ID,
                     FUN=run_lm_by_station,
                     selected_col=selected_col,
                     x_var_name=x_var_name,
                     y_var_name=y_var_name,
                     lf=lf,
                     log_val=T,
                     out_dir=out_dir,
                     out_suffix=out_suffix,
                     num_cores=1,
                     mc.preschedule=FALSE,
                     mc.cores = num_cores)
save(run_lm_obj,
     file=file.path(out_dir,paste0("run_lm_obj_",data_type,"_",out_suffix,".RData")))
save(run_lm_obj,
     file=file.path(out_dir,paste0("test_",data_type,"_",out_suffix,".RData")))

### TO DO CHANGE THE TRANSFORMATION TO LOG10 OR SOME OTHER!!!

run_lm_obj[[10]]$plot
run_lm_obj[[10]]$tb_coefficients
run_lm_obj[[10]]

run_lm_obj[[11]]$plot
run_lm_obj[[11]]$tb_coefficients
run_lm_obj[[11]]

#l_tb_coef <- try(lapply(run_lm_obj,FUN=function(x){x[["tb_coefficients"]]}))
#test <- mclapply(run_lm_obj[1701:1711],
#         FUN=function(x){x[["tb_coefficients"]]},
#                            mc.preschedule=FALSE,
#                   mc.cores = 11)
test <- mclapply(run_lm_obj,
                 FUN=function(x){x[["tb_coefficients"]]},
                 mc.preschedule=FALSE,
                 mc.cores = 11)
l_tb_coef <- test
names(l_tb_coef) <- list_ID
l_tb_coef_NA <- remove_from_list_fun(l_tb_coef)$list
tb_coef_combined <- do.call(rbind.fill,l_tb_coef_NA) #create a df for NA tiles with all accuracy metrics

#add the station ID to table by extracting name and replicating rows
list_ID_col_tmp <- unlist(lapply(1:length(l_tb_coef_NA),
                    FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},
                    x=l_tb_coef_NA,y=names(l_tb_coef_NA)))
#adding tile id summary data.frame
tb_coef_combined$ID_stat <- list_ID_col_tmp
tb_coef_combined$data_type <- data_type
tb_coef_combined$transform_fun <- transform_fun


##test<- merge(tb_coef_combined,dat_stat_location_DMR,by.x=c("ID_stat"),by.y=var_ID,all=T,suffixes=c("","_y"))
#> dim(test)
#[1] 3789   17
pos_col <- match("ID_stat",names(dat_stat_location_DMR))

df_coef_combined<- merge(tb_coef_combined,dat_stat_location_DMR[,-pos_col],by.x=c("ID_stat"),by.y=var_ID,suffixes=c("_x","_y"))
#> dim(test2)
#[1] 3203   17

write.table(tb_coef_combined,paste("df_coef_combined_station_",out_suffix,".txt",sep=""),sep=",")
#out_suffix_str <- paste(selected_ID,"_",out_suffix,sep="")

##
#tb_slope <- subset(tb_coef_combined,coef_type="slope",select())
tb_slope <- df_coef_combined[df_coef_combined$coef_type=="slope",]

histogram(tb_slope$p)
tb_slope$p_005 <- as.numeric(tb_slope$p <= 0.05)
tb_slope$p_001 <- as.numeric(tb_slope$p <= 0.01)
histogram(tb_slope$p_001,main="p <= 0.01")
histogram(tb_slope$p_005,main="p <= 0.05")
sum(tb_slope$p_001,na.rm=T)
sum(tb_slope$p_005,na.rm=T)
length(tb_slope$p_005)
sum(is.na(tb_slope$p))

table(tb_slope$p_005)

histogram(tb_slope$n)
mean(tb_slope$n)
range(tb_slope$n)
breaks_val <- seq(-1000,1000,by=100)
range_val <- range(tb_slope$estimate)

breaks_val <- c(range_val[1],breaks_val,range_val[2])
histogram(tb_slope$estimate)
histogram(tb_slope$estimate,breaks=breaks_val)
#histogram(tb_slope$estimate,xlim=c(-1000,1000))
#hist(tb_slope$estimate,breaks=breaks_val)

hist(tb_slope$estimate, xlim=c(-1000,1000),breaks=breaks_val)
barplot(table(tb_slope$estimate))
max(tb_slope$estimate)
min(tb_slope$estimate)

########################### End of script #####################################

#### Now run by year and station as a test???


# ##run by year
# 
# #test_mod <- run_simple_lm(test[[10]],y_var_name,x_var_name,log_val=T)
# l_tb_coef <- mclapply(l_df,
#                       FUN=run_simple_lm,
#                       y_var_name=y_var_name,
#                       x_var_name=x_var_name,
#                       log_val=T,
#                       mc.preschedule=FALSE,
#                       mc.cores = num_cores)
# l_dates <- 2003:2016
# names(l_tb_coef) <- l_dates
# l_tb_coef_NA <- remove_from_list_fun(l_tb_coef)$list
# 
# tb_coef_combined <- do.call(rbind.fill,l_tb_coef_NA) #create a df for NA tiles with all accuracy metrics
# 
# #add the dates to table
# dates_tmp <- lapply(1:length(l_tb_coef_NA),
#                     FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},
#                     x=l_tb_coef_NA,y=names(l_tb_coef_NA))
# #adding tile id summary data.frame
# tb_coef_combined$year <- dates_tmp


