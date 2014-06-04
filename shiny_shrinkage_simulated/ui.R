
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
    sliderInput("obs", 
                "Amount of data per subject (percentage):", 
                min = 0, 
                max = 1, 
                value = 0.5),
    sliderInput("standev", 
                "Standard deviation:", 
                min = 0, 
                max = 200, 
                value = 100),
    sliderInput("nsubjects", 
                "Number of subjects:", 
                min = 0, 
                max = 100, 
                value = 50)
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
   plotOutput("effectPlot")
  )
))