# R script for
# Reinhold Kliegl, Ping Wei, Michael Dambacher, Ming Yan, & Xiaolin Zhou (in press)
# Experimental Effects and Individual Differences in Linear Mixed Models:
# Estimating the Relation between Spatial, Object, and Attraction Effects in Visual Attention.
# Frontiers in Quantitative Psychology and Measurement

# Scripts for visualizing conditional modes are adopted from earlier versionds of Bates (2010)
# 5 Dec 2010, R. Kliegl


# Titus von der Malsburg pointed out two errors in the publication relating to AIC and BIC values
# reported on page 7:
#  (1)  The AIC-value for the model m2 was reported as: 328540; the correct value is: 325840.
#       This was a transposition error (85 instead of 58). 
#  (2)  The BIC-value for model m1 (325941) is actually smaller than the BIC-value for model m2 (325964).
#       Thus, for BIC model m2 fits worse than model m1.
# 5 Jan 2011, R. Kliegl

# Compute various rt transformations using boxcox(); new script for computing within-subject SEs
# 7 April 2011, R. Kliegl

# A few additions:
# ... model tests for power transformed RTs, using lambda from boxcox()
# ... use model matrix for tests of var components (instead of vectors c1, c2, c3)
# ... plot residuals 
# 9 April 2011, R. Kliegl

# A few revisions
# ... added package requests (reshape, scales) to accommodate change in ggplot2)
# ... replace mean(x) and sd(x) with sapply(x, mean) and sapply(x, sd)
# 21 April 2012, R. Kliegl

library(MASS)
library(ggplot2)
library(lme4)
library(gtools)
library(reshape)
library(scales)
library(plyr)
library(lattice)

rm(list=ls())

load("inst/KWDYZ.FQPM.rda")

############################
### Preliminary analyses ###
############################

# Data screening
# ... practice, catch, and no-response trials have already been removed in an earlier step

# ... remove extreme scores (best identified in LMM residuals)
#ibad <- which(d$rt < MIN | d$rt > MAX)
#d <- d[-ibad, ]

dat <- c[c$rt>150, 1:5]

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

# Transformations
# ... determine lambda (i.e., power coefficient); boxcox() is from MASS
# ... ... for exact lambda
lambdaList <- boxcox(rt ~ id*tar, data=dat)
lambda <- lambdaList$x[which.max(lambdaList$y)]  # 0.4242424

# ... for graph of profile of likelihood only 
boxcox(rt ~ id*tar, data=dat)

# ... alternative DVs
dat$srt <- dat$rt/1000
dat$lrt <- log(dat$rt)
dat$qrt <- sqrt(dat$rt)         # close lambda, from boxcox: 0.424 ~ 0.5 ~ sqrt(rt)
dat$prt <- dat$rt^(lambda)      # exact lambda 

cor(cbind(dat$srt, dat$lrt, dat$qrt, dat$prt))
str(dat)

# ... visualize distribution 
# ... ... using lattice package
boxplot(srt ~ tar, data=dat)
boxplot(lrt ~ tar, data=dat)
boxplot(qrt ~ tar, data=dat)
boxplot(prt ~ tar, data=dat)

# ... ... using reshape/ggplot2 package; transform from wide to long format, including type of transformation as factor "DV"
dat.rs <- melt(dat, id.vars=c( "id", "tar"), measure.vars=c("rt"), variable="DV")
qplot(x=tar, y=value, data=dat.rs, geom="boxplot", xlab="Cue-Target Relation", ylab="Reaction Time") + 
  facet_grid(DV ~ ., scales="free_y")

# Remove between-subject variance
# ...aggregate to subject-defined cell means
dat.id <- data.frame(cast(dat.rs, DV + id + tar ~ ., function(x) c(rt=mean(x), SD=sd(x), N=length(x))))

# ... normalize to GM = 0
dat.id <- ddply(dat.id, .(id, DV), transform, rt.w = rt - mean(rt))  # for rt.w everybody has mean of zero

# ... compute cell means
M  <- cast(dat.id, value="rt", DV+tar ~ ., 
           function (x) c(rt=mean(x), SD=sd(x), N=length(x)) )

# ... compute standard error for RTs AFTER removal of between-subject variance (within-subject SE) and add to M
M$SE <- cast(dat.id, value="rt.w", DV+tar ~ ., function(x) sd(x)/sqrt(length(x)) )[ , 3]  # only third column

# ... nicer labels for factor levels                  
levels(M$DV)  <- c("untransformed [s]", "logarithmic [ms]", "square root [ms]", "lambda [-s^(0.424)]")
levels(M$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")

# ... generate figure for profile of condition means (not in paper)
(plot1 <- qplot(x=tar, y=rt, data=M, xlab = "Cue-Target Relation", ylab = "Reaction Time", 
            geom=c("line", "point"), group=1) +
facet_grid(DV ~ ., scales="free_y") + 
geom_errorbar(aes(max=rt+2*SE, min=rt-2*SE, width=0.1)) + theme_bw()  )

######################################################
### Repeated-measures Multiple Regression Analyses ###
######################################################
N.id <- length(unique(dat$id))

# rmMRA for RT (see Table 1, top right)
(within.subject.coeff <- coef(lmList(rt ~ tar | id, data=dat)))
# ... same as
(within.subject.coeff.vectors <- coef(lmList(rt ~ c1 + c2 + c3 | id, data=dat)))

(M.rmMRA <- sapply(within.subject.coeff, mean))
(SD.rmMRA <- sapply(within.subject.coeff, sd))
(SE.rmMRA <- SD.rmMRA/sqrt(N.id))
(t.rmMRA <- M.rmMRA / SE.rmMRA)

# rmMRA for log RT (see Table 1, bottom right)
within.subject.coeff.lrt <- coef(lmList(lrt ~ tar | id, data=dat))
(M.rmMRA.lrt <- sapply(within.subject.coeff.lrt, mean))
SD.rmMRA.lrt <- sapply(within.subject.coeff.lrt, sd)
(SE.rmMRA.lrt <- SD.rmMRA.lrt/sqrt(N.id))
(t.rmMRA.lrt <- M.rmMRA.lrt / SE.rmMRA.lrt)

# rmMRA for powertransformed RT (not in paper)
within.subject.coeff.prt <- coef(lmList(prt ~ tar | id, data=dat))
(M.rmMRA.prt <- sapply(within.subject.coeff.prt, mean))
SD.rmMRA.prt <- sapply(within.subject.coeff.prt, sd)
(SE.rmMRA.prt <- SD.rmMRA.prt/sqrt(N.id))
(t.rmMRA.prt <- M.rmMRA.prt / SE.rmMRA.prt)

###########################
### Linear Mixed Models ###
###########################

#### DV: rt
# ... baseline model (only random intercept)
print(m0 <- lmer(rt ~ 1 + tar + (1 | id), data=dat))
m0@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters

# ... add variance components for effects
print(m1 <- lmer(rt ~ 1 + c1 + c2 + c3 + (1 | id) + (0 + c1 | id) + (0 + c2 | id) + (0 + c3 | id), data=dat))
m1@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters

# ... fully parameterized model; add correlation parameters (see Tables 1 and 2, top left)
print(m2 <- lmer(rt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat))  
m2@dims   # p = n of fixed effects; q = n of random effects; np = n of variance component parameters

# ... ... an alternative is to extract vectors c1, c2, and c3 from model matrix of m0
m2.x <- lmer(rt ~ tar + (model.matrix(m0)[ ,2:4] | id), data=dat)  # m2.xtra and m2 yield same estimates

# ... ... this allows very efficient tests of individual variance components, for example:
for (i in 1:3) {
  cat("\n\n\n", "Testing:", colnames(cmat)[i], "\n")
  .mx <- lmer(rt ~ tar + (1 | id) + (0 + model.matrix(m0)[,i+1] | id), data=dat)
  print(anova(m0, .mx))
}

# ... check significance of variance components (m0 vs m1) and significance of correlation parameters (m1 vs. m2)
anova(m0, m1, m2)

# ... check residuals
qqmath(resid(m2))
plot(fitted(m2), resid(m2))
abline(h=0)

### DV: log(rt)
# ... baseline model (only random intercept)
print(m0.lrt <- lmer(lrt ~ 1 + tar + (1 | id), data=dat))

# ... add variance components for effects
print(m1.lrt <- lmer(lrt ~ 1 + c1 + c2 + c3 + (1 | id) + (0 + c1 | id) + (0 + c2 | id) + (0 + c3 | id), data=dat))

# ... fully parameterized model; add correlation parameters (see Tables 1 and 2, bottom left)
print(m2.lrt <- lmer(lrt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat))  # equivalent to m1, m1a

# ... check significance of variance components (m0.lrt vs m1.lrt) and significance of correlation parameters (m1.lrt vs. m2.lrt)
anova(m0.lrt, m1.lrt, m2.lrt)

# ... check residuals
qqmath(resid(m2.lrt))
plot(fitted(m2.lrt), resid(m2.lrt))
abline(h=0)

### DV: power transformed (rt) (not in paper)
# ... baseline model (only random intercept)
print(m0.prt <- lmer(prt ~ 1 + tar + (1 | id), data=dat))

# ... add variance components for effects
print(m1.prt <- lmer(prt ~ 1 + c1 + c2 + c3 + (1 | id) + (0 + c1 | id) + (0 + c2 | id) + (0 + c3 | id), data=dat))

# ... fully parameterized model; add correlation parameters 
print(m2.prt <- lmer(prt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | id), data=dat), cor=FALSE)  # equivalent to m1, m1a

# ... check significance of variance components (m0.prt vs m1.prt) and significance of correlation parameters (m1.prt vs. m2.prt)
anova(m0.prt, m1.prt, m2.prt)

# ... check residuals
qqmath(resid(m2.prt))
plot(fitted(m2.prt), resid(m2.prt), pch=".")
abline(h=0)

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
# Note: The following was fixed to run again; revisions may have broken the original code (21/04/2012; rk)
dat.id.srt <- subset(dat.id, DV=="srt")
dat.id.srt$RT <- 1000*dat.id.srt$rt   # Convert to ms

GM <- mean(tapply(dat.id.srt$RT, dat.id.srt$id, mean))
dat.id.srt$RT.w <- 1000*dat.id.srt$rt.w + GM  # Convert to ms

# Assign subjects to speed groups based on median of condition means
within.subject.coeff$Speed <- quantcut(within.subject.coeff$'(Intercept)', q=seq(0, 1, 1/2))  
levels(within.subject.coeff$Speed) <- c("fast", "slow")

# ... add id as factor
within.subject.coeff$id <- factor(1:61)

# ... combine
dat.id.srt <- merge(dat.id.srt, within.subject.coeff, by.x="id", by.y="id")


# Post-hoc ANOVA and follow-ups 
summary(m1     <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id.srt))
summary(m1.spt <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id.srt, subset=tar=="val" | tar=="sod"))
summary(m1.obj <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id.srt, subset=tar=="dos" | tar=="sod"))
summary(m1.att <- aov(RT ~ Speed*tar + Error(id/tar), data=dat.id.srt, subset=tar=="dos" | tar=="dod"))

# Aggregate to condition x speed-group level--between-subject variance is in the data
(M1.id <- cast(melt(dat.id.srt, id.var=c("id","tar", "Speed"), measure.var="RT"), tar + Speed ~ ., 
                function(x) c(M=mean(x), SE=sd(x)/sqrt(length(x)), N=length(x) ) ) ) 

levels(M1.id$tar) <- c("valid",
                      "same object -\ndifferent location", 
                      "different object -\nsame location",
                      "different object -\ndiagonal location")


# Plot 
(Figure4a <- qplot(x=tar, y=M, shape=Speed, group=Speed, data=M1.id, ylim=c(300, 500), geom=c("point", "line"),
                 xlab="Cue-target relation", ylab="Response time [ms]") +
             geom_errorbar(aes(max=M+2*SE, min=M-2*SE, width=0.1)) + geom_point(size=3) +
             scale_shape_discrete("Speed group") + theme_bw(base_size=12) + opts(legend.position = c(0.3, 0.81)) )

# Aggregate to condition x speed-group level--between-subject variance is removed from the data
(M1.id.w <- cast(melt(dat.id.srt, id.var=c("id","tar", "Speed"), measure.var="RT.w"), tar + Speed ~ ., 
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

dat.id.srt <- merge(dat.id.srt, within.subject.coeff[ , c("Attraction", "id")], by.x="id", by.y="id")

# Post-hoc ANOVA and follow ups; leave out DOD rts (too redundant w/ Attraction factor)
summary(m2     <- aov(RT ~ Attraction*tar + Error(id/tar), data=dat.id.srt, subset=tar != "dod"))
summary(m2.spt <- aov(RT ~ Attraction*tar + Error(id/tar), data=dat.id.srt, subset=tar=="val" | tar=="sod"))
summary(m2.obj <- aov(RT ~ Attraction*tar + Error(id/tar), data=dat.id.srt, subset=tar=="dos" | tar=="sod"))

# Aggregate to condition x attraction-group level--between-subject variance is in the data
(M2.id <- cast(melt(dat.id.srt, id.var=c("id","tar", "Attraction"), measure.var="RT"), tar + Attraction ~ ., 
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
(M2.id.w <- cast(melt(dat.id.srt, id.var=c("id","tar", "Attraction"), measure.var="RT.w"), tar + Attraction ~ ., 
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
