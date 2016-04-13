##############################################  NEST Beach closure project  #######################################
###########################################  SHINY APP #######################################
#This script is the server side for the Shiny app to explore rainfall and bacteria measurements.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/12/2016 
#DATE MODIFIED: 04/12/2016
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
  datasetInput <- reactive({
    switch(input$dataset,
           "MHB" = data_df_MHB,
           "DMR" = data_df_DMR)
  })
  #output$data_df <- reactive({
  #  data_df<- datasetInput()
  #})
  output$dataset <- reactive({
    dataset<- datasetInput()
  })
  
  output$summary <- renderPrint({
    dataset <- datasetInput()
    summary(dataset)
  })
  
  ## First prepare plot for the station profile!!
  output$plot_ts <- renderPlot({
    #input$newplot
    # Add a little noise to the cars data
    #cars2 <- cars + rnorm(nrow(cars))
    #plot(cars2)
    #plot_ts <- 
    data_type <- input$dataset
    
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
    
    id_name<-input$station
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
    par(new=TRUE)              # key: ask for new plot without erasing old
    #plot(x,y,type="l",col=t_col[k],xlab="",ylab="",lty="dotted",axes=F) #plotting fusion profile
    #plot(log(df2$COL_SCORE),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
    plot(log(df2[[var_name]]),pch=10,cex=2.5,col="red", axes=F,ylab="",xlab="")
    
    #points(d_z2$COL_SCORE,col="red",pch=10,cex=2)
    legend("topleft",legend=c("stations"), 
           cex=1.2,col="red",pch =10,bty="n")
    
    axis(4,cex=1.2)
    mtext(4, text = "coliform scores", line = 3)
    #title(paste("Station time series",id_name,sep=" "))
    
    
  })
  
  output$raster_map <- renderPlot({
    #if (input$date_in == "Map1")
    #  data <- raster("file1.asc")
    #else if (input$variable == "Map2")
    #  data <- raster("file2.asc")
    #levelplot(data, margin=FALSE, par.settings=GrTheme)
    
    start_date <-input$dates[1]
    end_date <-input$dates[2]
    
    year_processed <- year(start_date) #make this work for end dates too!!
    in_dir_rst <- grep(paste0("prism_ppt_",year_processed), list_dir_rainfall,value=T)
    r_rainfall <- stack(mixedsort(list.files(pattern="*.tif",path=in_dir_rst,full.names=T))) #rainfall time series stack
    
    
    if (convert_to_inches==TRUE){
      r_rainfall <- r_rainfall/25.4 #improve efficiency later? YES!!
    }
    
    if(input$dataset=="MHB"){
      dat_stat <- dat_stat_location_MHB
      #rm(dat_stat_location_MHB)
      #rm(dat_stat_location_DMR)
    }
    if(input$dataset=="DMR"){
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