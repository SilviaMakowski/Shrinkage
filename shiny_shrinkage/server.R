
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
    
    dat_part <- NULL
    
    for (i in ids) {
      sub_i <- dat[dat$id==i, 1:5]
      for (j in conds) {
        cond_j <- sub_i[sub_i$cond1==j, 1:5]
        rnd_subset <- cond_j[sample(c(1:nrow(cond_j)),round(input$obs*nrow(cond_j))),1:5]
        dat_part <- rbind(dat_part, rnd_subset)
      }
    }
 
    dat <- dat_part
    
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
    
    ###########################
    ### Linear Mixed Models ###
    ###########################
    
    #### DV: rt
    # ... fully parameterized model; add correlation parameters (see Tables 1 and 2, top left)
    print(m2 <- lmer(rt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat))  
    
    
    # ... splom of coefficiencts (i.e., fixed effects plus conditional modes)
    m2.coef <- unlist(coef(m2))
    dim(m2.coef) <- c(61, 4)
    m2.coef <- data.frame(m2.coef)
    names(m2.coef) <- c("Mean", "Spatial", "Object", "Attraction")
    
    ###########################################################################
    ### Figure 3: Compare conditional modes and within-subject OLS effects  ###
    ###########################################################################
    
    # Within-subject OLS estimates
    df <- coef(lmList(rt ~ tar | id, data=dat))
    
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
    print(px1, position=c(0, 0, 0.5, 1), more=TRUE)
    print(px2, position=c(0.5, 0, 1, 1))
  })
  
})
