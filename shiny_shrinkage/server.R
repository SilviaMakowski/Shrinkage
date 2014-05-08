
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
library(MASS)
library(ggplot2)
library(lme4)
library(gtools)
library(reshape)
library(scales)

library(plyr)
library(lattice)

load("../inst/KWDYZ.FQPM.rda")

shinyServer(function(input, output) {
   
  output$effectPlot <- renderPlot({
    dat <- c[c$rt>150, 1:5]
    
    ids <- unique(c$id)
    conds <- unique(c$cond1)
    
    dat_half <- NULL
    
    for (i in ids) {
      for (j in conds) {
        new <- dat
        new <- new[new$id==i, 1:5]
        new <- new[new$cond1==j, 1:5]
        new <- new[1:round(0.25*nrow(new)), 1:5]
        dat_half <- rbind(dat_half, new)
      }
    }
    dat <- dat_half
    
    names(dat)[3:4] <- c("tar", "dir")
    dat$id <- factor(dat$id)
    dat$item <- factor(dat$item) 
    
    # Factors  - reorder levels for conditions
    dat$tar <- factor(dat$tar, levels=c("val", "sod", "dos", "dod")) 
    contrasts(dat$tar)
    
    # Contrasts
    # ... target factor
    cmat <- contr.sdif(4) 
    cmat[,3] <- -1*cmat[,3]  # invert gravity contrast
    rownames(cmat) <- c("val", "sod", "dos", "dod")
    colnames(cmat) <- c(".spt", ".obj", ".att")
    contrasts(dat$tar) <- cmat
    
    # ... also as vectors for testing individual random effects in LMM;  alternative use of model.matrix() below
    dat$c1 <- ifelse(dat$tar=="val", -0.75, 0.25)
    dat$c2 <- ifelse(dat$tar=="val" | dat$tar=="sod", -0.5, +0.5)
    dat$c3 <- ifelse(dat$tar=="dod", -0.75, 0.25)
    
    #First plot
    # Transformations
    # ... determine lambda (i.e., power coefficient); boxcox() is from MASS
    # ... ... for exact lambda
    lambdaList <- boxcox(rt ~ id*tar, data=dat)
    
  })
  
})
