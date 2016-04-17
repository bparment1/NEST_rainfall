##############################################  NEST Beach closure project  #######################################
#################################  Shiny Application for exploration of station measurements  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#This script is the client, user interface side of the shiny app.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/10/2016 
#DATE MODIFIED: 04/13/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -  Adding DMR data
#          - 
#TO DO:

#################################################################################################

###Loading R library and packages                                                      


library(shiny)

# Define UI for dataset viewer application
shinyUI(fluidPage(
  
  # Application title.
  titlePanel("Station Measurements"),
  
  # Sidebar with controls to select a dataset and specify the
  # number of observations to view. The helpText function is
  # also used to include clarifying text. Most notably, the
  # inclusion of a submitButton defers the rendering of output
  # until the user explicitly clicks the button (rather than
  # doing it immediately when inputs change). This is useful if
  # the computations required to render output are inordinately
  # time-consuming.
  sidebarLayout(
    
    sidebarPanel(
      selectInput("dataset", "Choose a dataset:", 
                  choices = c("MHB", "DMR")),
      #if(dataset=="MHB"){
      #  station_ID <- data_df_MHB$LOCATION_ID
      #}
      #if(dataset=="DMR"){
      #  station_ID <- data_df_DMR$LOCATION_ID
      #}
      #selectInput("station", "Choose a station:", 
      #            choices = station_ID),  
      #selectInput("station", "Choose a station:", 
      #            choices = unique(data_df_MHB$LOCATION_ID)),  
      #selectInput("station", "Choose a station:", 
      #            choices = unique(data_df$LOCATION_ID)),  
      #selectInput("station", "Choose a station:", 
      #         choices = unique(data_df$LOCATION_ID)),  
      uiOutput("choose_dataset"),
      
      uiOutput("choose_columns"),
      #selectInput("station", "Choose a station:", 
      #            choices = output$station_ID),  
      #selectInput("station", "Choose a station:", 
      #            choices = unique(input$dataset$LOCATION_ID)),        
      #numericInput("obs", "Number of observations to view:", 10),
      
      helpText("Note: the plot will show the specified",
               "station profile based on rainfall", 
               "on the full dataset."),
      #helpText("Note: while the data view will show only the specified",
      #         "number of observations, the summary will still be based",
      #         "on the full dataset."),
      
      #dateRangeInput("dates", label = h3("Date range")),
      dateRangeInput("dates", start=start_date,end=end_date,label = h3("Date range")),
      submitButton("Update View")
    ),
    
    mainPanel(
      plotOutput("plot_ts"),
      plotOutput("raster_map")
      #verbatimTextOutput("summary")
    )#,
    
    # Show a summary of the dataset and an HTML table with the
    # requested number of observations. Note the use of the h4
    # function to provide an additional header above each output
    # section.
    #mainPanel(
    #  h4("Summary"),
    #  verbatimTextOutput("summary"),
    #  
    #  h4("Observations"),
    # tableOutput("view")
    #)
  )
))

