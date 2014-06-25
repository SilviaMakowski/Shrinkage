
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel(HTML("<h3>Shrinkage effects</h3> 
              <p>This application can be used for illustrating shrinkage of observed data toward the estimated population mean in linear mixed models. LMM analyses are based on simulated data using mixedDesign() function (Hohenstein & Kliegl, 2013). Means, standard deviations and correlations used in mixedDesign() are calculated from experimental data by Kliegl, Wei, Dambacher, Yan and Zhou (2011) using the two-rectangle cueing paradigm (Egly, Driver & Rafal, 1994). </p>"), 
              windowTitle="Shiny Shrinkage"),
  
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
                "Amount of observations per within subject factor:", 
                min = 0, 
                max = 100, 
                value = 50),
    sliderInput("standevSubject", 
                "Standard deviation for selected subject (red):", 
                min = 0, 
                max = 200, 
                value = 100),
    sliderInput("numberObsSubject", 
                "Amount of observations for selected subject (red):", 
                min = 0, 
                max = 100, 
                value = 50)
  ),
  
  # Show caterpillar plot and scatterplots
  mainPanel(
   plotOutput("effectPlot",width = "800px", height = "600px")
  )
))
