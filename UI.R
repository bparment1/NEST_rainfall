##############################################  NEST Beach closure project  #######################################
#################################  Shiny Application for exploration of station measurements  #######################################
#This script explores the correlation between rainfall events and beach closures due to bacteria outbreaks in Maine.
#The script uses time series analyes from R. 

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/10/2016 
#DATE MODIFIED: 03/11/2016
#Version: 1
#PROJECT: NEST beach closures            

#
#COMMENTS: -   
#          - 
#TO DO:

#################################################################################################

###Loading R library and packages                                                      


library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("NEST beach stations"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      #droplist of stations??
      #sliderInput("bins",
      #            "Number of bins:",
      #            min = 1,
      #            max = 50,
      #           value = 30)
      selectInput("dataset", "Choose a station:", 
                  choices = c("rock", "pressure", "cars")),
      
      #numericInput("obs", "Number of observations to view:", 10),
      
      helpText("Note: the plot will show the specified",
               "station profile based on rainfall, ,
               "on the full dataset."),
      
      submitButton("Update View")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("rainfallPlot")
    )
  )
))