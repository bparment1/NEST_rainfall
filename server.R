##############################################  NEST Beach closure project  #######################################
###########################################  SHINY APP #######################################
#This script is the server side for the Shiny app to explore rainfall and bacteria measurements.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/12/2016 
#DATE MODIFIED: 04/17/2016
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
shinyServer(function(input, output) {
  
  
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
    selectInput("dataset", "Data set", as.list(data_sets),selected="MHB")
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
  #dataset <- reactive({
  #  filename <- paste0("data_", input$date, ".Rdata")
  #  load(filename)
  #})
  
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
  
  output$summary <- renderPrint({
    data_df <- datasetInput()
    summary(data_df)
  })
  
  ## First prepare plot for the station profile!!
  output$plot_ts <- renderPlot({
    #input$newplot
    # Add a little noise to the cars data
    #cars2 <- cars + rnorm(nrow(cars))
    #plot(cars2)
    #plot_ts <- 
    #browser()
    #data_type <- input$dataset #this is null in the first run of the app
    if(is.null(input$dataset)){
      data_type <- "MHB"
    }else{
      data_type <- input$dataset
    }
    data_df <- datasetInput()
    #if(data_type=="MHB"){
    #  data_df <- data_df_MHB
    #  #rm(data_df_MHB)
    #  #rm(data_df_DMR)
    #}
    #if(data_type=="DMR"){
    #  data_df <- data_df_DMR
    #  #rm(data_df_MHB)
    #  #rm(data_df_DMR)
    #}
    
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
      
      plot(log(df2[[var_name]]),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
      
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
  
  output$raster_map <- renderPlot({
    #if (input$date_in == "Map1")
    #  data <- raster("file1.asc")
    #else if (input$variable == "Map2")
    #  data <- raster("file2.asc")
    #levelplot(data, margin=FALSE, par.settings=GrTheme)
    
    start_date <-input$dates[1]
    end_date <-input$dates[2]
    #browser()
    year_processed <- year(start_date) #make this work for end dates too!!
    in_dir_rst <- grep(paste0("prism_ppt_",year_processed), list_dir_rainfall,value=T)
    r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=in_dir_rst,full.names=T))) #rainfall time series stack
    
    
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
    
    date_str <- paste0(year_processed,"-01-01") #change this later!!
    plot(r_rainfall,y=1,ext=extent(dat_stat),main=paste0("Raster image for ",date_str))
    plot(dat_stat,add=T,pch=3)
    #text(dat_stat,dat_stat$LOCATION_ID,cex=1.4)
    legend("topright",legend=paste0(data_type," stations"), 
           cex=1.2, col="black",pch =3,bty="n")
  })
  
  #dataInput <- reactive({
  #  getSymbols(input$symb, src = "yahoo", 
  #             from = input$dates[1],
  #             to = input$dates[2],
  #             auto.assign = FALSE)
  #})

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