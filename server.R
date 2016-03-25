##############################################  NEST Beach closure project  #######################################
###########################################  SHINY APP #######################################
#This script is the server side for the Shiny app to explore rainfall and bacteria measurements.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/12/2016 
#DATE MODIFIED: 03/25/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -   
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

  ## First prepare plot for the station profile!!
  output$plot_ts <- renderPlot({
    input$newplot
    # Add a little noise to the cars data
    #cars2 <- cars + rnorm(nrow(cars))
    #plot(cars2)
    #plot_ts <- 
    id_name<-input$station
    id_selected <- data_df[[var_ID]]==id_name
    
    ### Not get the data from the time series
    data_pixel <- data_df[id_selected,]
    data_pixel <- as.data.frame(data_pixel)
    r_ts_name <- names(r_rainfall)
    
    var_pix_ts <- t(as.data.frame(subset(data_pixel,select=var_name)))
    #pix <- t(data_pixel[1,24:388])#can subset to range later
    
    pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
    #pix_ts <- subset(as.data.frame(pix_ts),select=r_ts_name)
    pix_ts <- (as.data.frame(pix_ts))
    ## Process the coliform data
    d_z <- zoo(pix_ts,idx) #make a time series ...
    names(d_z)<- "rainfall"
    plot(d_z,lty=2,ylab="rainfall",xlab="Time",main=paste0("Rainfall for station ",id_name))
  })
  
  output$raster_map <- renderPlot({
    #if (input$date_in == "Map1")
    #  data <- raster("file1.asc")
    #else if (input$variable == "Map2")
    #  data <- raster("file2.asc")
    #levelplot(data, margin=FALSE, par.settings=GrTheme)
    plot(r_rainfall,y=1,ext=extent(dat_stat),main=paste0("Raster image for ",idx[1]))
    plot(dat_stat,add=T,pch=3)
    text(dat_stat,dat_stat$LOCATION_ID,cex=1.4)
    legend("topright",legend=c("stations"), 
           cex=1.2, col="black",pch =3,bty="n")
  })
  
  dataInput <- reactive({
    getSymbols(input$symb, src = "yahoo", 
               from = input$dates[1],
               to = input$dates[2],
               auto.assign = FALSE)
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