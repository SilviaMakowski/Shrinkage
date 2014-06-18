
# This is the server logic for Shiny shrinkage.

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

means <- matrix(c(358,392,406,403), ncol=4)
standevAcrossSubjects <- 60

# correlations for mixedDesign based on experimental data
correlations <- matrix(c(1, .91, .88, .84, 
                         .91, 1, .95, .96, 
                         .88, .95, 1, .96,
                         .84, .96, .96, 1), nrow=4, ncol=4, byrow=TRUE)

# Linear Mixed Model
build_model <- function(dat){
  model <- lmer(rt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat)
  return(model)
}

# Function for generating n data based on mean and standard deviation
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

simulate_data <- function(dat_mixed, standevWithin, numberObs) {
  dat_sim <-NULL
  
  sub_ids <- unique(dat_mixed[,1])
  
  # generate data for each condition based on means from mixedDesign 
  for (i in sub_ids){
    sub_i <- dat_mixed[dat_mixed$id==i, ]
    for (j in unique(sub_i[,2])){
      tar_j <- sub_i[sub_i$tar==j,]
      r <- rnorm2(numberObs,tar_j[,3],standevWithin)
      rmat <- matrix(c(rep(i,times=numberObs),rep(j,times=numberObs)),ncol=2,nrow=numberObs)
      rmati <- cbind(rmat, r)
      dat_sim <- rbind(dat_sim, rmati)
    }
  } 
  
  dat_sim <- data.frame(dat_sim)
  
  names(dat_sim)<- c("id","tar","rt")
  dat_sim$rt <- as.numeric(as.character(dat_sim$rt))
  dat_sim$tar <- factor(dat_sim$tar, levels=c("val", "sod", "dos", "dod")) 
  contrasts(dat_sim$tar)
  
  
  # vectors for testing individual random effects in LMM
  dat_sim$c1 <- ifelse(dat_sim$tar=="val", -0.75, 0.25)
  dat_sim$c2 <- ifelse(dat_sim$tar=="val" | dat_sim$tar=="sod", -0.5, +0.5)
  dat_sim$c3 <- ifelse(dat_sim$tar=="dod", -0.75, 0.25)
  return(dat_sim)
}


shinyServer(function(input, output) {
   
  output$effectPlot <- renderPlot({
    if (input$nobs<2) {
      stop('Not enough observations: random-effects parameters and residual variance unidentifiable. Minimum is 2.')
    }
    if (input$nsubjects<2) {
      stop('Not enough subjects for the given design. Minimum is 2.')
    }
    if (input$numberObsSubject<2) {
      stop('Not enough observations. Minimum is 2.')
    }
    
    # simulate means for each subject and condition
    set.seed(1)
    dat_mixed <- mixedDesign(W = 4, M = means, SD = standevAcrossSubjects, n = input$nsubjects, R=correlations, long = TRUE)
    
    names(dat_mixed)<- c("id","tar","M")
    levels(dat_mixed$tar) <- c("val","sod","dos","dod")
    
    # simulate observations within each subject and condition
    dat <- simulate_data(dat_mixed, input$standevWithin, input$nobs)
    
    # simulate data for selected subject according to user input from sliders 4 and 5
    dat_subject <- dat_mixed[dat_mixed$id==1,]
    dat_sim_Subject <- simulate_data(dat_subject, input$standevSubject, input$numberObsSubject)
    
    # update selected subject's data to user input
    dat <- dat[which(as.numeric(dat$id) > 1),]
    dat <- rbind(dat_sim_Subject, dat)
    
    # Contrasts
    # ... target factor
    cmat <- contr.sdif(4) 
    cmat[,3] <- -1*cmat[,3]  # invert gravity contrast
    rownames(cmat) <- c("val", "sod", "dos", "dod")
    colnames(cmat) <- c(".spt", ".obj", ".att")
    contrasts(dat$tar) <- cmat
    
    # build model on data of all subjects
    print(m2<- build_model(dat))  
    
    m2.coef <- unlist(coef(m2))
    dim(m2.coef) <- c(input$nsubjects, 4)

    m2.coef <- data.frame(m2.coef)
    names(m2.coef) <- c("Mean", "Spatial", "Object", "Attraction")
   
    m2.ranef <- ranef(m2, postVar = TRUE)
    names(m2.ranef[[1]])[1:4] <- c("Mean", "Spatial", "Object", "Attraction")
    

    # Within-subject OLS estimates
    df <- coef(lmList(rt ~ tar | id, data=dat))

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
                               # print selected subject in red
                               panel.points(df[1,4],df[1,2], col="red", cex=0.9, pch=1)
                               panel.points(df[1,8],df[1,6], col="red", cex=0.7, pch=19)
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
                               # print selected subject in red
                               panel.points(df[1,1],df[1,2], col="red", cex=0.9, pch=1)
                               panel.points(df[1,5],df[1,6], col="red", cex=0.7, pch=19)
                             }
                )) )  
    
    # Caterpillar plots for conditional modes: order on spatial effect
    px4 <- print(dotplot.RK(m2.ranef, refvar = 2, layout=c(4,1), scales = list(x = list(relation = 'free', rot=0), y = list(draw=FALSE)), strip = TRUE)[[1]][c(1,2,3,4)])
    
    print(px4, position=c(0, 0, 0.35, 1), more=TRUE)
    print(px1, position=c(0.1, 0.47, 1, 1), more=TRUE)
    print(px2, position=c(0.1, 0, 1, 0.53))
  })
    
})
  
