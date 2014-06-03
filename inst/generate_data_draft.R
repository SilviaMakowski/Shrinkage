# generate data...


library(MASS)
library(plyr)
library(ggplot2)

rm(list=ls())
source("inst/mixedDesign.R")

mean.mat <- matrix(c(358,392,406,403), ncol=4)
                   
# Call function 


simulate_data <- function(means, number, standev) {
  set.seed(1)
  dat <- mixedDesign(W = 4, M = means, SD = standev, n = number, R=.7, long = TRUE)
  names(dat)<- c("id","tar","M")
  levels(dat$tar) <- c("val","sod","dos","dod")
  
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
