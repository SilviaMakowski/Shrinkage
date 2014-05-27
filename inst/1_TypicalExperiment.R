## title:  *** mixed-model ANOVA of Factorial Design ***
## author: Reinhold Kliegl

#---------------------------------------------------------

# Example for mixedDesign.R
#
# MIXED-MODEL DESIGN 2 B x 3 W x 2 W x 12 subjects in each group
# between subject factor, e.g., type of nonwords (2 levels)
# within subject factor, e.g., eccentricity (2 levels)
# within subject factor, e.g., primetype    (3 levels)
#
#----------------------------------------------------------

library(ggplot2)
library(plyr)

rm(list=ls())
source("functions/mixedDesign.v0.6.2.R")  

# Generate a data set 
Means <- matrix(c(610, 625, 620, 640, 650, 670, 620, 635, 655, 670, 660, 680), 2, 6)

Means
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]  610  620  650  620  655  660
# [2,]  625  640  670  635  670  680

# One can also fill a matrix by row ...
Means.2 <- matrix(c(610, 625, 620, 640, 650, 670, 620, 635, 655, 670, 660, 680), 2, 6, byrow=TRUE)

Means.2
#       [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]  610  625  620  640  650  670
# [2,]  620  635  655  670  660  680


StDev <- matrix(rep(10, 12), 2, 6)


# ... long format
set.seed(1)
data <- mixedDesign(B=c(2), W=c(2, 3), M=Means, SD=10, R=.4, n=12, long=TRUE)
data


# Not needed here: Check need for transformation 
# (only preliminary; should be based on model residuals)
# Determine lambda (i.e., power coefficient); e.g., boxcox()  from MASS


# Figure of means
# removing between-subject variance for error bars  
GM <- mean(data$DV)
data <- ddply(data, .(id), transform, DV.w = DV - mean(DV) + GM)  

# compute the Morey factor for the adjustment of within-subject SEs (Morey, 2008)
nl <- nlevels(data$W_a)*nlevels(data$W_b)
mf <- sqrt( nl/(nl-1) )  

# compute nobs, means and sd's, normed means and normed sd's, se's, and ci's
table1 <- ddply(data, .(B_A, W_a, W_b), summarise, N=length(DV), M=mean(DV), M_norm=mean(DV.w),
	                                   SD=sd(DV), SD_normed=sd(DV.w), 
	                                   SE=SD_normed/sqrt(N)*mf,
	                                   CI=SE*qt(.975, N-1))

# Morey-based within-subject SEs and CIs using Winston Chang's functions
#source("functions/normDataWithin.R")
#source("functions/summarySE.R")
#source("functions/summarySEwithin.R")
#summarySEwithin(data=data, idvar="id", measurevar="DV", betweenvars="B_A", withinvars=c("W_a", "W_b"))

levels(table1$B_A) <- c("Coltheart", "Wuggy")
levels(table1$W_a) <- c("0", "4")
levels(table1$W_b) <- c("Morph.", "Ortho.", "Unrelat.")

# ... plot highest-order interaction
(
plot1 <- qplot(data=table1, x=W_b, y=M, linetype=W_a, group=W_a, facets= . ~ B_A,
      xlab="Stimulus condition", ylab="Response time [ms]", geom=c("line", "point") ) +
      scale_linetype_discrete("Eccentrícity") +
      geom_errorbar(aes(max = M + 2*SE, min = M - 2*SE), width=0) +  geom_point(size=2) +
      labs(title=" Morphological Priming", legend.position = "right") + theme_bw()
)

# Analysis of Variance
summary(aov(DV ~  B_A*W_a*W_b + Error(id/(W_a*W_b)), data=data))

# There is a significant W_a x W_b interaction 
# ... compute nobs, means and sd's, normed means and normed sd's, se's, and ci's
table2 <- ddply(data, .(W_a, W_b), summarise, N=length(DV), M=mean(DV), M_norm=mean(DV.w),
			  SD=sd(DV), SD_normed=sd(DV.w), 
			  SE=SD_normed/sqrt(N)*mf,
			  CI=SE*qt(.975, N-1))
levels(table2$W_a) <- c("0", "4")
levels(table2$W_b) <- c("Morph.", "Ortho.", "Unrelat.")

# ... plot significant interaction
(
plot2 <- qplot(data=table2, x=W_b, y=M, linetype=W_a, group=W_a, 
	      xlab="Stimulus condition", ylab="Response time [ms]", geom=c("line", "point") ) +
		scale_linetype_discrete("Eccentrícity") +
		geom_errorbar(aes(max = M + 2*SE, min = M - 2*SE), width=0) +  geom_point(size=2) +
		labs(title=" Morphological Priming", legend.position = "right") + theme_bw()
)
