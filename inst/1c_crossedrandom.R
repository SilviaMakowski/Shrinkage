# Example for illustration of mixedDesign() function for crossed random effects
# Hohenstein & Kliegl (2013): Simulation of Factorial Mixed-Model Designs in R: The mixedDesign() function
# 27 June 2013

library(plyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(MASS)

rm(list=ls())
source("functions/mixedDesign.v0.6.2.R")
source("functions/remef.v0.6.10.R")
source("functions/ggCaterpillar.R")

nsubj <- 25
nitem <- 30

# Cell means (between-subject factor levels across rows; within-subject factor levels across columns)
mean.mat <- matrix(rep(c(250, 260,  290, 280, 310, 340), nsubj), nrow=nsubj, ncol=6, byrow=TRUE)

set.seed(1)
dat <- mixedDesign(B = nsubj, W = c(2, 3), M = mean.mat, SD = 150, n = nitem, R=.42, long = TRUE, empirical=FALSE)

# Rename variables and levels
# ... convert Between-Factor to Subj, and id to Item Factor
dcast(dat, Subj ~ Item, length)
names(dat) <- c("Subj", "Item", "Distance", "Preview", "GD")
dat$Item <- factor((as.numeric(as.character(dat$Item))-1) %% nitem + 1)

# ... the two within-subject factors are also within-item factors
levels(dat$Distance) <- c("near", "far")
levels(dat$Preview) <- c("identical", "related", "unrelated")

dcast(dat, Subj ~ Item, length)
dcast(dat, Subj ~ Distance + Preview, mean)
dcast(dat, Item ~ Distance + Preview, mean)

# F1-ANOVA
d.subj <- ddply(dat, .(Subj, Distance, Preview), summarise, gd = mean(GD), n=length(GD))
summary(aov(gd ~ Distance*Preview + Error(Subj/(Distance*Preview)), data=d.subj))

# Compute table of means for interaction
table.subj <- ddply(d.subj, .(Distance, Preview), summarise, N=length(gd), M=mean(gd),
		   SD=sd(gd), SE=SD/sqrt(N)*sqrt(4/3) )   # Morey-factor 4 levels/3

# Plot table
qplot(data=table.subj, x=Preview, y=M, ylab = "Gaze duration (ms)", group=Distance, colour=Distance, 
      geom=c("point", "line")) +
	geom_errorbar(aes(ymax=M+2*SE, ymin=M-2*SE), width=0) + theme_bw()

# F2-ANOVA
d.item <- ddply(dat, .(Item, Distance, Preview), summarise, gd = mean(GD), n=length(GD))
summary(aov(gd ~ Distance*Preview + Error(Subj/(Distance*Preview)), data=d.subj))

# Compute table of means for interaction
table.item <- ddply(d.item, .(Distance, Preview), summarise, N=length(gd), M=mean(gd),
		   SD=sd(gd), SE=SD/sqrt(N)*sqrt(4/3) )   # Morey-factor 4 levels/3

# Plot table
qplot(data=table.item, x=Preview, y=M, ylab = "Gaze duration (ms)", group=Distance, colour=Distance, 
	geom=c("point", "line")) +
	geom_errorbar(aes(ymax=M+2*SE, ymin=M-2*SE), width=0) + theme_bw()

# LMM
contrasts(dat$Distance) <- contr.sdif(2)
contrasts(dat$Preview) <- contr.sdif(3)

summary(m1 <- lmer(GD ~ Distance*Preview + (1 | Subj) + (1 | Item), data=dat ), cor=F)

ggCaterpillar(ranef(m1, condVar=TRUE))

summary(m1a <- lmer(GD ~ Distance*Preview + (1 | Subj), data=dat ))
anova(m1a, m1, refit=FALSE)

summary(m1b <- lmer(GD ~ Distance*Preview + (1 | Item), data=dat ))
anova(m1b, m1, refit=FALSE)

summary(m2a <- lmer(GD ~ Distance*Preview + (1 | Subj) + (1 | Subj : Distance) + (1 | Item), data=dat))
anova(m1, m2a, refit=FALSE)
summary(m2b <- lmer(GD ~ Distance*Preview + (1 | Subj) + (1 | Item) + (1 | Item : Distance) , data=dat))
anova(m1, m2b, refit=FALSE)

summary(m3a <- lmer(GD ~ Distance*Preview + (1 | Subj) + (1 | Subj : Preview) + (1 | Item), data=dat ))
anova(m1, m3a, refit=FALSE)
summary(m3b <- lmer(GD ~ Distance*Preview + (1 | Subj) + (1 | Item) + (1 | Item : Preview) , data=dat ))
anova(m1, m3b, refit=FALSE)
