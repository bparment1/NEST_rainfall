
##############################################  NEST Beach closure project  #######################################
#################################  Analyses for exploration of station measurements  #######################################
#This contains functions used to explore the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#It uses time series processed earlier in R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 06/15/2016 
#DATE MODIFIED: 06/22/2016
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

run_simple_lm <- function(df,y_var_name,x_var_name,log_val=T,plot_fig=T,out_suffix="",out_dir="."){
  #
  #
  #http://stackoverflow.com/questions/9923722/how-to-make-lm-display-in-its-output-a-formula-passed-to-it-as-a-variable
  
  #test <- lm(log1p(df$COL_SCORE) ~ log1p(df$rainfall),df)
  
  if(log_val==T){ #change the log transform later!!
    df_tmp <- data.frame(df[[y_var_name]],df[[x_var_name]])
    df_tmp <- na.omit(df_tmp)
    names(df_tmp)<- c(paste("log1p",y_var_name,sep="_"),paste("log1p",x_var_name,sep="_"))
    y_var <- names(df_tmp[1])
    x_var <- names(df_tmp[2])

    formula_str <- formula(paste(y_var,"~",x_var,sep=" "))
    #formula_str <- as.character(paste(y_var,"~",x_var,sep=" "))
    
    #mod <- lm(log1p(df[[y_var_name]]) ~ log1p(df[[x_var_name]]),df_tmp)
    mod <- lm(eval(formula_str),df_tmp)
    #formula(lm(formula_str))
  }
  if(log_val!=T){
    df_tmp <- data.frame(df[[y_var_name]],df[[x_var_name]])
    df_tmp <- na.omit(df_tmp)
    names(df_tmp) <- c(y_var_name,x_var_name)
    y_var <- names(df_tmp[1])
    x_var <- names(df_tmp[2])
    formula_str <- formula(paste(y_var,"~",x_var,sep=" "))
    
    mod <- lm(eval(formula_str),df_tmp)
    #mod <- lm(df[[y_var_name]] ~ df[[x_var_name]],df_tmp)
  }

  ## NUmber of inputs??
  no_obs <- nobs(mod)
  #nrow(mod) #same as nobs output
  
  # Create the character string that you want to print
  tp <- sprintf("%s=%.1f + %.2f %s", all.vars(formula(mod))[1],
                coef(mod)[1], coef(mod)[2], all.vars(formula(mod))[2])
  
  summary_mod_tb <- (summary(mod))
  tb_coefficients <- as.data.frame(summary_mod_tb$coefficients)
  names(tb_coefficients) <- c("estimate","std_error","t_value","p")
  #rownames(tb_coefficients)
  if(nrow(tb_coefficients)>1){
    coef_type_val <- c("intercept","slope")
    p_intercept <- tb_coefficients$p[1]
    p_slope <- tb_coefficients$p[2]
    p_vals <- c(p_intercept,p_slope)
  }else{
    coef_type_val <- c("intercept")
    p_intercept <- tb_coefficients$p[1]
    p_vals <- p_intercept
  }
  
  tb_coefficients$coef_type <- coef_type_val
  tb_coefficients$n <- no_obs
  rownames(tb_coefficients) <- NULL
  p_vals <- paste(format(p_vals, digits=3),collpapse="")
  plot_obj <- xyplot(log1p(df[[y_var_name]]) ~ log1p(df[[x_var_name]]),
              xlab=list(label=x_var_name, cex=1.2),
              ylab=list(label=y_var_name, cex=1.2),
              auto.key=list(x=0.05,y=0.95,text=c(tp,p_vals,no_obs),
                            points=FALSE, lines=FALSE,col=c(1,1,1))
              )
  #format(tb, digits=3)
  layout_m <- c(1.5,1)
  if(plot_fig==T){
    png_filename <- file.path(out_dir,paste("Figure","_","scatter_plot_regression",out_suffix,".png", sep=""))
    png(png_filename,height=480*layout_m[2],width=480*layout_m[1])
    print(plot_obj)
    dev.off()
  }

  #change name of rows and add n columns for the number of inputs in the model!
  run_lm_obj <- list(tb_coefficients,mod,plot_obj,no_obs)
  names(run_lm_obj) <- c("tb_coefficients","mod","plot","nobs")  
  save(run_lm_obj,
       file=file.path(out_dir,paste0("run_lm_obj_station_",data_type,"_",out_suffix,".RData")))
  
  return(run_lm_obj)
}

run_lm_by_station <- function(selected_ID,selected_col,x_var_name,y_var_name,lf,log_val=T,out_dir,out_suffix,num_cores){
  ##Function to run lm model by station with or without log transform
  
  l_df <- mclapply(lf,
                   FUN=read_select_station,
                   selected_val=selected_ID,
                   selected_col=selected_col,
                   mc.preschedule=FALSE,
                   mc.cores = 1)
  
  #test_mod <- run_simple_lm(test[[10]],y_var_name,x_var_name,log_val=T)
  df_combined <- do.call(rbind,l_df) 
  write.table(df_combined,paste("df_combined_station_",selected_ID,"_",out_suffix,".txt",sep=""),sep=",")
  out_suffix_str <- paste(selected_ID,"_",out_suffix,sep="")
  
  if(nrow(df_combined)>0){
    #debug(run_simple_lm)
    run_mod_obj <- run_simple_lm(df=df_combined,
                                 y_var_name=y_var_name,
                                 x_var_name=x_var_name,
                                 log_val=T,
                                 plot_fig=T,
                                 out_suffix=out_suffix_str,
                                 out_dir=out_dir)
    
    
    #run_mod_obj$tb_coefficients
    #run_mod_obj$p
  }else{
    run_mod_obj <- NULL
  }
  save(run_mod_obj,
       file=file.path(out_dir,paste0("run_mod_obj_station_",data_type,"_",out_suffix_str,".RData")))
  return(run_mod_obj)
}

remove_from_list_fun <- function(l_x,condition_class ="try-error"){
  index <- vector("list",length(l_x))
  for (i in 1:length(l_x)){
    if (inherits(l_x[[i]],condition_class)){
      index[[i]] <- FALSE #remove from list
    }else{
      index[[i]] <- TRUE
    }
  }
  l_x<-l_x[unlist(index)] #remove from list all elements using subset
  
  obj <- list(l_x,index)
  names(obj) <- c("list","valid")
  return(obj)
}

###Modifying function lm to print formula if it is given as string or formula
#this is taken from:
lm <- function(...) {
  mf <- match.call()
  mf[[1]] <- quote(stats::lm)
  env <- parent.frame()
  fm <- eval(mf, env)
  fm$call$formula <- formula(fm)
  fm
}

########################### End of script #####################################
