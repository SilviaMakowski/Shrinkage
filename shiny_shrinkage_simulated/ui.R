
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
    sliderInput("standevAcross", 
                "Standard deviation across subjects:", 
                min = 0, 
                max = 100, 
                value = 50),
    sliderInput("standevWithin", 
                "Standard deviation within subjects:", 
                min = 0, 
                max = 100, 
                value = 50),
    sliderInput("obs", 
                "Amount of data per subject (percentage of 30):", 
                min = 0, 
                max = 1, 
                value = 0.5)
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
   plotOutput("effectPlot")
  )
))
