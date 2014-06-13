
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Spatial and attraction effects"),
  
  # Sidebar with a slider input for amount of data per subject
  sidebarPanel(
    sliderInput("nsubjects", 
                "Number of subjects:", 
                min = 0, 
                max = 100, 
                value = 50),
    sliderInput("standevWithin", 
                "Standard deviation within subjects:", 
                min = 0, 
                max = 200, 
                value = 100),
    sliderInput("nobs", 
                "Amount of observations per subject and condition:", 
                min = 0, 
                max = 50, 
                value = 25),
    sliderInput("standevSubject", 
                "Standard deviation for subject:", 
                min = 0, 
                max = 200, 
                value = 100),
    sliderInput("numberObsSubject", 
                "Amount of observations for subject and condition:", 
                min = 0, 
                max = 50, 
                value = 25)
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
   plotOutput("effectPlot")
  )
))
