
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel(HTML("<h1>Shrinkage effects</h1> 
              <p>Illustration of shrinkage towards population mean using the two-rectangle cueing paradigm (Egly, Driver & Rafal, 1994).</p>"), 
              windowTitle="Shiny shrinkage"),
  
  # Sidebar with slider input
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
   plotOutput("effectPlot",width = "800px", height = "600px")
  )
))
