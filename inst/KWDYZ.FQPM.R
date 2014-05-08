# R script for

# Reinhold Kliegl, Ping Wei, Michael Dambacher, Ming Yan, & Xiaolin Zhou (2011)
# Experimental Effects and Individual Differences in Linear Mixed Models:
# Estimating the Relationship between Spatial, Object, and Attraction Effects in Visual Attention.
# Frontiers in Quantitative Psychology and Measurement

# doi: 10.3389/fpsyg.2010.00238

# Scripts for visualizing conditional modes are adopted from earlier versions of Bates (2010)
# 5 Dec 2010, R. Kliegl

# Corrections. Titus von der Malsburg pointed out two errors in the publication relating to AIC and BIC values reported on page 7:
# (1)  The AIC-value for the model m2 was reported as: 328540; the correct value is: 325840. 
#      This was a transposition typo ("85" instead of "58").
# (2)  The BIC-value for model m1 (325941) is actually smaller than the BIC-value for model m2 (325964).
#      Thus, for BIC model m2 the fit is not better than the one for model m1.
# 5 Jan 2011, R. Kliegl


library(MASS)
library(gtools)
library(ggplot2)
library(lme4)

rm(list=ls())

load("KWDYZ.FQPM.rda")

############################
### Preliminary analyses ###
############################

# Data screening
dat <- c[c$rt>150, 1:5]
names(dat)[3:4] <- c("tar", "dir")
dat$id <- factor(dat$id)
dat$item <- factor(dat$item)

# Factors  - reorder levels for conditions
dat$tar <- factor(dat$tar, levels=c("val", "sod", "dos", "dod")) 
contrasts(dat$tar)

# Contrasts
# ... target factor
cmat.t <- contr.sdif(4) 
ct <- contr.treatment(4)     # contrast matrix (contr.sdif from MASS)
cmat.t[,3] <- -1*cmat.t[,3]  # invert gravity contrast
rownames(cmat.t) <- c("val", "sod", "dos", "dod")
colnames(cmat.t) <- c(".spt", ".obj", ".att")
dat$tar <- C(dat$tar, cmat.t, 3)
contrasts(dat$tar)

# ... also as vectors for testing individual random effects
dat$c1 <- ifelse(dat$tar=="val", -0.75, 0.25)
dat$c2 <- ifelse(dat$tar=="val" | dat$tar=="sod", -0.5, +0.5)
dat$c3 <- ifelse(dat$tar=="dod", -0.75, 0.25)

######################################################
### Repeated-measures Multiple Regression Analyses ###
######################################################

# rmMRA for RT (see Table 1, top right)
within.subject.coeff <- coef(lmList(rt ~ tar | id, data=dat))
within.subject.coeff <- coef(lmList(rt ~ c1 + c2 + c3 | id, data=dat))

(M.rmMRA <- mean(within.subject.coeff))
(SD.rmMRA <- sd(within.subject.coeff))
(SE.rmMRA <- SD.rmMRA/sqrt(61))
(t.rmMRA <- M.rmMRA / SE.rmMRA)

# rmMRA for log RT (see Table 1, bottom right)
lrt <- log(dat$rt)
within.subject.coeff.log <- coef(lmList(lrt ~ tar | id, data=dat))
(M.rmMRA.log <- mean(within.subject.coeff.log))
SD.rmMRA.log <- sd(within.subject.coeff.log)
(SE.rmMRA.log <- SD.rmMRA.log/sqrt(61))
(t.rmMRA.log <- M.rmMRA.log / SE.rmMRA.log)

within.subject.coeff$id <- factor(1:61)

##############################
### Figure 0--not in paper ###
##############################

# ... aggregate to subject x condition level
dat.rs <- melt(dat, id=c("id", "tar"), measure=c("rt"), na.rm=TRUE)
dat.id <- data.frame(cast(dat.rs, id + tar ~ ., function(x) c(RT=mean(x), N=length(x) ) ))

# ... remove between-subject variance for error bars  
GM <- mean(tapply(dat.id$RT, dat.id$id, mean))
dat.id <- ddply(dat.id, .(id), transform, RT.w = RT - mean(RT) + GM)  

# ...  aggregate to condition level
(M0.id.w <- cast(melt(dat.id, id.var=c("id","tar"), measure.var="RT.w"), tar  ~ ., 
                function(x) c(M=mean(x), SE=sd(x)/sqrt(length(x)), N=length(x) ) ) ) 
				# ... x labels 
levels(M0.id.w$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")

# ... plot
(Figure0 <- qplot(x=tar, y=M, group=1, data=M0.id.w, , ylim=c(320, 420), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
                 geom_errorbar(aes(max=M+1.96*SE, min=M-1.96*SE, width=0.1)) + geom_point(size=3) +
                 theme_bw(base_size=11) )
#dev.copy2pdf(file="/Users/kliegl/documents/manuscripts/ZDK/Figures/Figure0.Means.pdf", width=5, height=5)


###########################
### Linear Mixed Models ###
###########################

#### DV: rt
# Baseline model (random-intercept)
print(m0 <- lmer(rt ~ tar + (1 | id), data=dat))
m0@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters

# Test significance of variances
print(m1 <- lmer(rt ~ c1 + c2 + c3 + (1 | id) + (0 + c1 | id) + (0 + c2 | id) + (0 + c3 | id), data=dat))
m1@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters

# Fully parameterized model (see Tables 1 and 2, bottom left)
print(m2 <- lmer(rt ~ c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat))  
m2@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters

anova(m0, m1, m2)

### DV: log(rt)
# ... fully parameterized model (see Tables 1 and 2, bottom left)
print(m2.log <- lmer(log(rt) ~ c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat))  # equivalent to m1, m1a
m2.log@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters


########################################################
### Figure 2: Caterpillar plots of conditional modes ###
########################################################

#  Extra plots (not in the paper)
# ... splom of coefficiencts (i.e., fixed effects plus conditional modes)
m2.coef <- unlist(coef(m2))
dim(m2.coef) <- c(61, 4)
m2.coef <- data.frame(m2.coef)
names(m2.coef) <- c("Mean", "Spatial", "Object", "Attraction")

splom(~ m2.coef[,1:4])

# ... histograms of effects only
m2.coef.long <- melt(m2.coef, measure=names(m2.coef)[2:4])
histogram( ~ value | variable, data = m2.coef.long,  
          xlab="Conditional modes", type = "density",
          panel = function(x, ...) {
              panel.histogram(x, ...)
              panel.mathdensity(dmath = dnorm, col = "black",
                                args = list(mean=mean(x),sd=sd(x)))
          } )

# ... scatterplot of conditional modes
m2.ranef <- ranef(m2, postVar = TRUE)
names(m2.ranef[[1]])[1:4] <- c("Mean", "Spatial", "Object", "Attraction")
plot(m2.ranef, aspect = 1, type = c("g", "p"))[[1]]   # plot.mer

# ... conditional modes sorted by intercept
trellis.device(color = FALSE) # set to black and white

dotplot(m2.ranef, layout=c(4,1), scales = list(x = list(relation = 'same', rot=0), y=list(draw=FALSE)), strip = TRUE)[[1]][c(1,2,3,4)]
# ... ... default: relation = "same" 

source("dotplot.RK.R") 
# Figure used in paper: order on spatial effect
dotplot.RK(m2.ranef, refvar = 2, layout=c(4,1), scales = list(x = list(relation = 'free', rot=0), y = list(draw=FALSE)), strip = TRUE)[[1]][c(1,2,3,4)]

# Alternative: use same x-scale for all panels
# dotplot.RK(m2.ranef, refvar = 2, layout=c(4,1), scales = list(x = list(relation = 'same', rot=0), y = list(draw=FALSE)), strip = TRUE)[[1]][c(1,2,3,4)]

###########################################################################
### Figure 3: Compare conditional modes and within-subject OLS effects  ###
###########################################################################

# Within-subject OLS estimates
df <- coef(lmList(rt ~ tar | id, data=dat))
dat$lrt <-log(dat$rt)
df.log <- coef(lmList(lrt ~ tar | id, data=dat))
mean(df)
cor(df)
mean(df.log)
cor(df.log)

op <- options(digits=2)
cor(m2.coef[, 1:4])
sd(m2.coef)
options(op)

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

# Arrange the two plots vertically
print(px1, position=c(0, 0.475, 1, 1), more=TRUE)
print(px2, position=c(0, 0, 1, 0.525))


# Other shrinkage plots, not included in paper
# ... mean RT and object effects
with(df,
    print(xyplot(tar.obj ~ `(Intercept)`, aspect = 1,
           x1 = Mean, y1 = Object, xlab="Mean RT", ylab="Object Effect",
           panel = function(x, y, x1, y1, subscripts, ...) {
               panel.grid(h = -1, v = -1)
               x1 <- x1[subscripts]
               y1 <- y1[subscripts]
               panel.arrows(x, y, x1, y1, type = "open", length = 0.1, col="gray50", lty=1, angle = 15, ...)
               panel.points(x, y, col="black", cex=0.5, pch=1)
               panel.points(x1, y1, col="black", cex=0.5, pch=19)
           },
                  key = list(space = "top", columns = 2,
                  text = list(c("Linear mixed model", "Within-group")),
                  points = list(col = trellis.par.get("superpose.symbol")$col[1:2],
                  pch = trellis.par.get("superpose.symbol")$pch[1:2]))
              )))

# ... mean RT and attraction effects
with(df,
    print(xyplot(tar.att ~ `(Intercept)`, aspect = 1,
           x1 = Mean, y1 = Attraction, xlab="Mean RT", ylab="Attraction Effect",
           panel = function(x, y, x1, y1, subscripts, ...) {
               panel.grid(h = -1, v = -1)
               x1 <- x1[subscripts]
               y1 <- y1[subscripts]
               panel.arrows(x, y, x1, y1, type = "closed", length = 0.1, col="gray50", lty=1, angle = 15, ...)
               panel.points(x, y, col="black", cex=0.5, pch=1)
               panel.points(x1, y1, col="black", cex=0.5, pch=19)
           },
                  key = list(space = "top", columns = 2,
                  text = list(c("Linear mixed model", "Within-group")),
                  points = list(col = trellis.par.get("superpose.symbol")$col[1:2],
                  pch = trellis.par.get("superpose.symbol")$pch[1:2]))
              )))


######################################################
### Figure 4: Speed-group plot and post-hoc ANOVAs ###
######################################################

# Assign subjects to speed groups based on median of condition means
within.subject.coeff$Speed <- quantcut(within.subject.coeff$'(Intercept)', q=seq(0, 1, 1/2))  
levels(within.subject.coeff$Speed) <- c("fast", "slow")

dat.id <- merge(dat.id, within.subject.coeff, by.x="id", by.y="id")

# Post-hoc ANOVA and follow-ups 
summary(m1     <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id))
summary(m1.spt <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id, subset=tar=="val" | tar=="sod"))
summary(m1.obj <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id, subset=tar=="dos" | tar=="sod"))
summary(m1.att <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id, subset=tar=="dos" | tar=="dod"))

# Aggregate to condition x speed-group level--between-subject variance is in the data
(M1.id <- cast(melt(dat.id, id.var=c("id","tar", "Speed"), measure.var="RT"), tar + Speed ~ ., 
                function(x) c(M=mean(x), SE=sd(x)/sqrt(length(x)), N=length(x) ) ) ) 

levels(M1.id$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")


# Plot 
(Figure4a <- qplot(x=tar, y=M, shape=Speed, group=Speed, data=M1.id, ylim=c(300, 500), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]")
             geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3) +
             scale_shape_discrete("Speed group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81)) )

# Aggregate to condition x speed-group level--between-subject variance is removed from the data
(M1.id.w <- cast(melt(dat.id, id.var=c("id","tar", "Speed"), measure.var="RT.w"), tar + Speed ~ ., 
                function(x) c(M=mean(x), SE=sd(x)/sqrt(length(x)), N=length(x) ) ) ) 

levels(M1.id.w$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")

# Plot--included as Figure 4 in paper
(Figure4b <- qplot(x=tar, y=M, shape=Speed, group=Speed, data=M1.id.w, ylim=c(320, 420), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
            geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3)  +
            scale_shape_discrete("Speed group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81)) )

# Put them together
M1 <- rbind(M1.id, M1.id.w)
M1$Panel <- factor(rep(c("With between-subject variance in RT", "W/o between-subject variance in RT" ), each=8))
M1$Panel <- relevel(M1$Panel, 2)

# Plot
(Figure4 <- qplot(x=tar, y=M, shape=Speed, group=Speed, data=M1, facets=Panel ~ ., 
	             ylim=c(300, 500), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
           geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3)  + 
           scale_shape_discrete("Speed group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81)) )

#dev.copy2pdf(file="/Users/kliegl/documents/manuscripts/ZDK/Figures/Figure4.Speed.pdf", width=5, height=5)

###########################################################
### Figure 5: Attraction-group plot and post-hoc ANOVAs ###
###########################################################

# Assign subjects to attraction groups based on tar.att >= 0 (i.e., DOD < DOS)
within.subject.coeff$Attraction <- ifelse(within.subject.coeff$tar.att >= 0, "yes", "no") 
within.subject.coeff$Attraction <- factor(within.subject.coeff$Attraction)

dat.id <- merge(dat.id, within.subject.coeff[ , c("Attraction", "id")], by.x="id", by.y="id")

# Post-hoc ANOVA and follow ups; leave out DOD rts (too redundant w/ Attraction factor)
summary(m2     <- aov(RT ~ Attraction*tar + Error(id/tar), data=dat.id, subset=tar != "dod"))
summary(m2.spt <- aov(RT ~ Attraction*tar + Error(id/tar), data=dat.id, subset=tar=="val" | tar=="sod"))
summary(m2.obj <- aov(RT ~ Attraction*tar + Error(id/tar), data=dat.id, subset=tar=="dos" | tar=="sod"))

# Aggregate to condition x attraction-group level--between-subject variance is in the data
(M2.id <- cast(melt(dat.id, id.var=c("id","tar", "Attraction"), measure.var="RT"), tar + Attraction ~ ., 
                function(x) c(M=mean(x), SE=sd(x)/sqrt(length(x)), N=length(x) ) ) ) 

levels(M2.id$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")

# Plot 
(Figure5a <- qplot(x=tar, y=M, shape=Attraction, group=Attraction, data=M2.id, ylim=c(300, 500), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
            geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3)  + 
            scale_shape_discrete("Attraction group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81)) )

# Aggregate to condition x attraction-group level--between-subject variance is removed from the data
(M2.id.w <- cast(melt(dat.id, id.var=c("id","tar", "Attraction"), measure.var="RT.w"), tar + Attraction ~ ., 
                function(x) c(M=mean(x), SE=sd(x)/sqrt(length(x)), N=length(x) ) ) ) 

levels(M2.id.w$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")

# Plot--included as Figure 5 in paper
(Figure5b <- qplot(x=tar, y=M, shape=Attraction, group=Attraction, data=M2.id.w, ylim=c(320, 420), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
            geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3)  + 
            scale_shape_discrete("Attraction group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81))  )

# Put them together
M2 <- rbind(M2.id, M2.id.w)
M2$Panel <- factor(rep(c("With between-subject variance in mean RT", "W/o between-subject variance in mean RT" ), each=8))
M2$Panel <- relevel(M2$Panel, 2)

# Plot
(Figure5 <- qplot(x=tar, y=M, shape=Attraction, group=Attraction, data=M2, facets=Panel ~ ., 
	             ylim=c(300, 500), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
           geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3)  + 
           scale_shape_discrete("Attraction group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81)) )

#dev.copy2pdf(file="/Users/kliegl/documents/manuscripts/ZDK/Figures/Figure5.Attraction.pdf", width=5, height=5)
