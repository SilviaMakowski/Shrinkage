# Shrinkage

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

# Data screening
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

# Transformations
# ... determine lambda (i.e., power coefficient); boxcox() is from MASS
# ... ... for exact lambda
lambdaList <- boxcox(rt ~ id*tar, data=dat)
lambda <- lambdaList$x[which.max(lambdaList$y)]  # 0.4242424

# ... alternative DVs ##### maybe later.... 
dat$srt <- dat$rt/1000
dat$lrt <- log(dat$rt)
dat$qrt <- sqrt(dat$rt)         # close lambda, from boxcox: 0.424 ~ 0.5 ~ sqrt(rt)
dat$prt <- dat$rt^(lambda)      # exact lambda

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

## new data dat2

