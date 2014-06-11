
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
source("../inst/mixedDesign.R")


# Linear Mixed Model
build_model <- function(dat){
  model <- lmer(rt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat)
  return(model)
}

rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

means <- matrix(c(358,392,406,403), ncol=4)

simulate_data <- function(means, standevAcross, standevWithin, number) {
  set.seed(1)
  
  dat <- mixedDesign(W = 4, M = means, SD = standevAcross, n = number, R=.7, long = TRUE)
 
  names(dat)<- c("id","tar","M")
  levels(dat$tar) <- c("val","sod","dos","dod")
  
  dat_sim <-NULL
  
  for (i in 1:nrow(dat)){
    sub_i <- dat[dat$id==i, ]
    for (j in unique(sub_i[,2])){
      tar_j <- sub_i[sub_i$tar==j,]
      r <- rnorm2(30,tar_j[,3],standevWithin)
      rmat <- matrix(c(rep(i,times=30),rep(j,times=30)),ncol=2,nrow=30)
      rmati <- cbind(rmat, r)
      dat_sim <- rbind(dat_sim, rmati)
    }
  } 
  
  dat <- data.frame(dat_sim)
  
  names(dat)<- c("id","tar","rt")
  dat$rt <- as.numeric(dat$rt)
  dat$tar <- factor(dat$tar, levels=c("val", "sod", "dos", "dod")) 
  contrasts(dat$tar)
  
  
  # ... also as vectors for testing individual random effects in LMM;  alternative use of model.matrix() below
  dat$c1 <- ifelse(dat$tar=="val", -0.75, 0.25)
  dat$c2 <- ifelse(dat$tar=="val" | dat$tar=="sod", -0.5, +0.5)
  dat$c3 <- ifelse(dat$tar=="dod", -0.75, 0.25)
  return(dat)
}


shinyServer(function(input, output) {
   
  output$effectPlot <- renderPlot({
    dat <- simulate_data(means, input$standevAcross, input$standevWithin, input$nsubjects)
    ids <- unique(dat$id)
    conds <- unique(dat$tar)
    
    dat_part <- NULL
    
    set.seed(1)
    
    for (i in ids) {
      sub_i <- dat[dat$id==i,]
      for (j in conds) {
        cond_j <- sub_i[sub_i$tar==j,]
        rnd_subset <- cond_j[sample(c(1:nrow(cond_j)),round(input$obs*nrow(cond_j))),]
        dat_part <- rbind(dat_part, rnd_subset)
      }
    }
    
    # Contrasts
    # ... target factor
    cmat <- contr.sdif(4) 
    cmat[,3] <- -1*cmat[,3]  # invert gravity contrast
    rownames(cmat) <- c("val", "sod", "dos", "dod")
    colnames(cmat) <- c(".spt", ".obj", ".att")
    contrasts(dat_part$tar) <- cmat
    
    print(m2<- build_model(dat_part))  
    
    m2.coef <- unlist(coef(m2))
    dim(m2.coef) <- c(input$nsubjects, 4)

    m2.coef <- data.frame(m2.coef)

    names(m2.coef) <- c("Mean", "Spatial", "Object", "Attraction")
   
    m2.ranef <- ranef(m2, postVar = TRUE)

    names(m2.ranef[[1]])[1:4] <- c("Mean", "Spatial", "Object", "Attraction")
    

    # Within-subject OLS estimates
    df <- coef(lmList(rt ~ tar | id, data=dat_part))

    m2.coef$id <- factor(1:input$nsubjects)

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
    px4 <- print(dotplot.RK(m2.ranef, refvar = 2, layout=c(4,1), scales = list(x = list(relation = 'free', rot=0), y = list(draw=FALSE)), strip = TRUE)[[1]][c(1,2,3,4)])
    
    print(px4, position=c(0, 0, 0.2, 1), more=TRUE)
    print(px1, position=c(0.25, 0, 0.6, 1), more=TRUE)
    print(px2, position=c(0.65, 0, 1, 1))
  })
    
  })
  
