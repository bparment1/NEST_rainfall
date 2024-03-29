##############################################  NEST Beach closure project  #######################################
###########################################  SHINY APP #######################################
#This script is the server side for the Shiny app to explore rainfall and bacteria measurements.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/12/2016 
#DATE MODIFIED: 04/26/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: - Adding DMR data   
#          - 
#TO DO:
#
#################################################################################################

library(shiny)
library(datasets)

# Define server logic required to summarize and view the 
# selected dataset
shinyServer(function(input, output,session) {
  
  
  # Return the requested dataset
  #datasetInput <- reactive({
  #  switch(input$dataset,
  #         "rock" = rock,
  #         "pressure" = pressure,
  #         "cars" = cars)
  #})
  #var_ID <- "LOCATION_ID"
  #data_sets <- c("mtcars", "morley", "rock")
  data_sets <- c("MHB","DMR")
  output$choose_dataset <- renderUI({
    selectInput("dataset", "Dataset", as.list(data_sets),selected="MHB")
  })
  
  datasetInput <- reactive({
    #inputs:
    #data_type <- input$dataset
    if(is.null(input$dataset)){
      data_type <- "MHB"
    }else{
      data_type <- input$dataset
    }
    #data_type <- input$choose_dataset
    start_date <-input$dates[1]
    end_date <-input$dates[2]
    
    #get years processed
    year_processed_start <- year(start_date) #make this work for end dates too!!
    year_processed_end <- year(end_date) #make this work for end dates too!!
    #browser()
    #For faster reading of the data...multicore and reduce data in memory by on the fly loading of data
    list_year_processed <- year_processed_start:year_processed_end
    list_tb <- mclapply(1:length(list_year_processed),function(i){
      filename_str <- file.path("./data",paste0("data_df_combined_",list_year_processed[i],"_" ,data_type,".txt"));
      tb <- read.table(filename_str,sep=",",header=T,fill=T,stringsAsFactors = F) #bacteria measurements
    },mc.preschedule=FALSE,mc.cores = num_cores)
    data_df<- do.call(rbind,list_tb)
    data_df
    #browser()
  })

  # Check boxes
  output$choose_columns <- renderUI({
    # If missing input, return to avoid error later in function
    #if(is.null(input$dataset))
    #  #return()
    #  data_df <- get(data_df)
    #browser()
    # Get the data set with the appropriate name
    #dat <- get(input$dataset)
    data_df <- datasetInput()
    #colnames <- names(dat)
    #Warning: Error in [[: subscript out of bounds
    # Create the checkboxes and select them all by default
    #checkboxGroupInput("columns", "Choose columns", 
    #                   choices  = colnames,
    #                   selected = colnames)
    #
    #list_station <- c("17","")
    list_location_ID <- unique(data_df$LOCATION_ID)
    selectInput("station", "Choose a station:", 
             choices = list_location_ID,selected= "4" ) 
  })

  #datasetInput <- reactive({
  #  data_type <- input$dataset
  #  if(data_type=="MHB"){
  #   data_df <- data_df_MHB
  #   #rm(data_df_MHB)
  #   #rm(data_df_DMR)
  #  }
  #  if(data_type=="DMR"){
  #   data_df <- data_df_DMR
  #   #rm(data_df_MHB)
  #  #rm(data_df_DMR)
  # }
  #  data_df
  #})
  #data_type <- input$dataset
  ##try
  #data_type[1]
  
  
  #if(data_type=="MHB"){
  #  data_df <- data_df_MHB
  #  #rm(data_df_MHB)
  #  #rm(data_df_DMR)
  #}
  #if(data_type=="DMR"){
  #  data_df <- data_df_DMR
  #  #rm(data_df_MHB)
  # #rm(data_df_DMR)
  #}
  
  ##Reactive function to return desired dataset
  #datasetInput <- reactive({
  #  switch(input$dataset,
  #         "MHB" = data_df_MHB,
  #         "DMR" = data_df_DMR)
  #})
  data_df <- reactive({
    data_df<- datasetInput()
  })
  #output$data_df <- reactive({
  #  data_df<- datasetInput()
  #})
  #output$station_ID <- reactive({
  #  data_df <- datasetInput()
  #  station_ID <- unique(data_df$LOCATION_ID)
  #  station_ID
  #})
  
  #0utput$dataset <- reactive({
  #  dataset<- datasetInput()
  #})
  
  ## First prepare plot for the station profile!!
  output$plot_ts <- renderPlot({
    #browser()
    #data_type <- input$dataset #this is null in the first run of the app
    if(is.null(input$dataset)){
      data_type <- "MHB"
    }else{
      data_type <- input$dataset
    }
    
    #Get the relevant data
    data_df <- datasetInput()
    
    if(data_type=="MHB"){
      var_name <- var_name_MHB
    }else{
      var_name <- var_name_DMR
    }

    #id_name<-input$station #is null in the first run
    if(is.null(input$station)){
      id_name <- "4" # station 4 in default data "MHB"
    }else{
      id_name <- input$station
    }
    id_selected <- data_df[[var_ID]]==id_name
    
    ### Not get the data from the time series
    data_pixel <- data_df[id_selected,]
    data_pixel$rainfall <- as.numeric(data_pixel$rainfall)
    d_z_tmp <-zoo(data_pixel$rainfall,as.Date(data_pixel$date))
    #names(d_z_tmp)<- "rainfall"
    #data_pixel <- as.data.frame(data_pixel)
    d_z_tmp2 <- zoo(data_pixel[[var_name]],as.Date(data_pixel$date))
    
    start_date <-input$dates[1]
    end_date <-input$dates[2]
    
    d_z <- window(d_z_tmp,start=start_date,end=end_date)
    d_z2 <- window(d_z_tmp2,start=start_date,end=end_date)
    df2 <- as.data.frame(d_z2)
    names(df2)<- var_name
    #d_z <- subset(d_z,select=c("rainfall"))
    #d_z$rainfall <- as.numeric(as.character(d_z$rainfall))
    #r_ts_name <- names(r_rainfall)
    
    #var_pix_ts <- t(as.data.frame(subset(data_pixel,select=var_name)))
    #pix <- t(data_pixel[1,24:388])#can subset to range later
    
    #pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
    #pix_ts <- subset(as.data.frame(pix_ts),select=r_ts_name)
    #pix_ts <- (as.data.frame(pix_ts))
    ## Process the coliform data
    #d_z <- zoo(pix_ts,idx) #make a time series ...
    #names(d_z_tmp)<- "rainfall"
    #test<-as.data.frame(d_z)
    #plot(test$rainfall,lty=2,ylab="rainfall",xlab="Time",main=paste0("Rainfall for station ",id_name))
    #(start_date)
    plot(d_z,lty=2,ylab="rainfall",xlab="Time",main=paste0("Rainfall from ",start_date," to ", end_date,
                                                    " for station ",id_name))
    
    #abline(h=threshold_val,col="green")
    ##Added to deal with bugs ylim, in some cases, tere are no bacteria measurements...then skip
    nb_na <- sum(is.na(df2[[var_name]]))
    if(nb_na<length(df2[[var_name]])){
      par(new=TRUE)              # key: ask for new plot without erasing old
      #plot(x,y,type="l",col=t_col[k],xlab="",ylab="",lty="dotted",axes=F) #plotting fusion profile
      #plot(log(df2$COL_SCORE),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
      
      plot(log10(df2[[var_name]]),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
      
      #points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
      legend("topleft",legend=c("stations"), 
             cex=1.2,col="red",pch =10,bty="n")
      
      axis(4,cex=1.2)
      mtext(4, text = "bacteria scores", line = 3)
      #title(paste("Station time series",id_name,sep=" "))
    }else{
      mtext("No bacteria measurements for the time period selected")
    }
  })
  
  output$summary <- renderPrint({
    data_df <- datasetInput()
    summary(data_df)
  })
  
  output$raster_map <- renderPlot({

    #browser()
    r_start_date <-input$r_dates[1]
    r_end_date <-input$r_dates[2]
    #browser()
    
    year_processed_start <- year(r_start_date) #make this work for end dates too!!
    year_processed_end <- year(r_end_date) #make this work for end dates too!!
    
    #For faster reading of the data...multicore and reduce data in memory by on the fly loading of data
    list_year_processed <- year_processed_start:year_processed_end
    #year_processed <- year(start_date) #make this work for end dates too!!
    #in_dir_rst <- grep(paste0("prism_ppt_",year_processed), list_dir_rainfall,value=T)
    l_in_dir_rst <- unlist(lapply(list_year_processed, FUN=function(x){grep(paste0("prism_ppt_",x), list_dir_rainfall,value=T)}))
    r_rainfall_filename <- (mixedsort(list.files(pattern="*.tif",path=l_in_dir_rst,full.names=T))) #rainfall time series stack
    
    #r_year_processed <- year(r_start_date) #make this work for end dates too!!
    #in_dir_rst <- grep(paste0("prism_ppt_",r_year_processed), list_dir_rainfall,value=T)
    #This maybe fast in a brick??
    #maybe use glob from Sys instead
    #  lf_mosaic <- lapply(1:length(day_to_mosaic),FUN=function(i){
    #searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
    #print(searchStr)
    #Sys.glob(searchStr)
    #use option quick to make the stack creation faster!! there is no check on the extent and projection
    r_rainfall <- stack(r_rainfall_filename,quick=T,filename="r_rainfall.tif") #rainfall time series stack, maybe use a temp file??
    #writeRaster(r_rainfall, filename="r_rainfall.tif", overwrite=TRUE) #Takes more than 3 minues to write to file so stopped it
    #r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=in_dir_rst,full.names=T))) #rainfall time series stack
    
    if (convert_to_inches==TRUE){
      r_rainfall <- r_rainfall/25.4 #improve efficiency later? YES!!
    }
    ### Deal with bug, length zero argument below
    if(is.null(input$dataset)){
      data_type <- "MHB"
    }else{
      data_type <- input$dataset
    }
    
    if(data_type=="MHB"){
      dat_stat <- dat_stat_location_MHB
      #rm(dat_stat_location_MHB)
      #rm(dat_stat_location_DMR)
    }
    if(data_type=="DMR"){
      dat_stat <- dat_stat_location_DMR
      #rm(dat_stat_location_MHB)
      #rm(dat_stat_location_DMR)
    }
    #browser()
    #date_str <- paste0(year_processed,"-01-01") #change this later!!
    date_str <- input$date_range
    diff_date <- as.numeric((date_str - r_start_date) +1)
    #plot(r_rainfall,y=1,ext=extent(dat_stat),main=paste0("Raster image for ",date_str))
    title_str <- paste0("Rainfall for ",date_str)
    plot(r_rainfall,y=diff_date,ext=extent(dat_stat),main=title_str)
    #Update the plot as well here!!
    plot(dat_stat,add=T,pch=3)
    #text(dat_stat,dat_stat$LOCATION_ID,cex=1.4)
    legend("topright",legend=paste0(data_type," stations"), 
           cex=1.2, col="black",pch =3,bty="n")
  })
  
  # Generate a summary of the dataset
  #output$summary <- renderPrint({
  #  dataset <- datasetInput()
  #  summary(dataset)
  #})
  
  # Show the first "n" observations
  #output$view <- renderTable({
  #  head(datasetInput(), n = input$obs)
  #})
})

###################### END OF SCRIPT ##########################
#Find out about animations
#http://moc.environmentalinformatics-marburg.de/doku.php?id=courses:msc:data-management:code-examples:dm-ce-12-02
#http://search.r-project.org/library/rglwidget/html/playwidget.html
#http://stackoverflow.com/questions/27194893/reset-animation-in-shiny-r-studio
#http://shiny.rstudio.com/articles/action-buttons.html
#http://www.r-bloggers.com/shiny-module-design-patterns-pass-module-inputs-to-other-modules/?utm_source=feedburner&utm_medium=email&utm_campaign=Feed%3A+RBloggers+%28R+bloggers%29
#http://www.r-bloggers.com/animated-plots-with-r/
#http://stackoverflow.com/questions/30492537/shiny-app-embedded-video-not-working-browser-dependent


#sliderInput("date_range", 
#            "Choose Date Range:", 
#            min = as.Date("2016-02-01"), max = Sys.Date(), 
#            value = c(as.Date("2016-02-25"), Sys.Date())
#)
