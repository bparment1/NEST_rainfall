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
  output$plot <- renderPlot({
    data <- getSymbols(input$symb, src = "yahoo", 
                       from = input$dates[1],
                       to = input$dates[2],
                       auto.assign = FALSE)
  plot_ts <- plot(d_z,lty=2,ylab="rainfall",xlab="Time",main="")
  # Generate a summary of the dataset
  output$summary <- renderPrint({
    dataset <- datasetInput()
    summary(dataset)
  })
  
  # Show the first "n" observations
  output$view <- renderTable({
    head(datasetInput(), n = input$obs)
  })
})