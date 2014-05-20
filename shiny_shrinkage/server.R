
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
source("../inst/dotplot.RK.R")

dat <- c[c$rt>150, 1:5]

# Linear Mixed Model
build_model <- function(dat){
  model <- lmer(rt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat)
  return(model)
}

build_data <- function(dat){
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
  return(dat)
}



shinyServer(function(input, output) {
   
  output$effectPlot <- renderPlot({ 
    
    ids <- unique(dat$id)
    conds <- unique(dat$cond1)
    
    dat_part <- NULL
    
    for (i in ids) {
      sub_i <- dat[dat$id==i, 1:5]
      for (j in conds) {
        cond_j <- sub_i[sub_i$cond1==j, 1:5]
        rnd_subset <- cond_j[sample(c(1:nrow(cond_j)),round(input$obs*nrow(cond_j))),1:5]
        dat_part <- rbind(dat_part, rnd_subset)
      }
    }
    
    ## Build model m2 with partial data
    dat_part <- build_data(dat_part)
    print(m2<- build_model(dat_part))  
    m2.coef <- unlist(coef(m2))
    dim(m2.coef) <- c(61, 4)
    m2.coef <- data.frame(m2.coef)
    names(m2.coef) <- c("Mean", "Spatial", "Object", "Attraction")
    ## px4
    m2.ranef <- ranef(m2, postVar = TRUE)
    names(m2.ranef[[1]])[1:4] <- c("Mean", "Spatial", "Object", "Attraction")
    
    
    # Within-subject OLS estimates
    df <- coef(lmList(rt ~ tar | id, data=dat_part))
    
    m2.coef$id <- factor(1:61)
    
    df <- cbind(df, m2.coef)
    
    # Spatial and attraction effects
    px1 <- with(df,
                print(xyplot(tar.spt ~ tar.att, aspect = 1,
                             x1 = Attraction, y1 = Spatial, xlab="Attraction Effect", ylab="Spatial Effect",
                             panel = function(x, y, x1, y1, subscripts, ...) {
                               panel.grid(h = -1, v = -1)
                               x1 <- x1[subscripts]
                               y1 <- y1[subscripts]
                               panel.arrows(x, y, x1, y1, type = "open", length = 0.1, col="gray50", lty=1, angle = 15, ...)
                               panel.points(x, y, col="black", cex=0.9, pch=1)
                               panel.points(x1, y1, col="black", cex=0.7, pch=19)
                             }
                ))     )
    # Mean RT and spatial effects
    px2 <- with(df,
                print(xyplot(tar.spt ~ `(Intercept)`, aspect = 1,
                             x1 = Mean, y1 = Spatial, xlab="Mean RT", ylab="Spatial Effect",
                             panel = function(x, y, x1, y1, subscripts, ...) {
                               panel.grid(h = -1, v = -1)
                               x1 <- x1[subscripts]
                               y1 <- y1[subscripts]
                               panel.arrows(x, y, x1, y1, type = "open", length = 0.1, col="gray50", lty=1, angle = 15, ...)
                               panel.points(x, y, col="black", cex=0.9, pch=1)
                               panel.points(x1, y1, col="black", cex=0.7, pch=19)
                             }
                )) )  
  
    # Figure used in paper: order on spatial effect
    px4 <- print(dotplot.RK(m2.ranef, refvar = 2, layout=c(1,1), scales = list(x = list(relation = 'free', rot=0), y = list(draw=FALSE)), strip = TRUE)[[1]][c(1)])
    
    print(px4, position=c(0, 0, 0.2, 1), more=TRUE)
    print(px1, position=c(0.25, 0, 0.6, 1), more=TRUE)
    print(px2, position=c(0.65, 0, 1, 1))
  })
  
  output$meanPlot <- renderPlot({ 
    ## Build model m1 with full data
    dat_full <- build_data(dat)
    m1<- build_model(dat_full)
    ## px3
    m1.ranef <- ranef(m1, postVar = TRUE)
    names(m1.ranef[[1]])[1:4] <- c("Mean", "Spatial", "Object", "Attraction")
    px3 <- print(dotplot.RK(m1.ranef, refvar = 2, layout=c(1,1), scales = list(x = list(relation = 'free', rot=0), y = list(draw=FALSE)), strip = TRUE)[[1]][c(1)])
    print(px3, position=c(0, 0, 0.2, 1))
    })
})
