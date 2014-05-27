## title: Simulation and computation of within-subject errors
## author: Reinhold Kliegl

# Script: withinSubjErrExample.R
# Date: 11. 12.2012

library(ggplot2)
library(reshape2)
library(plyr)

rm(list=ls())

source("functions/mixedDesign.v0.6.2.R")  

# Another example of mixedDesign(); this time involving means and correlations of IQ scores
# Generate table of means for Age (2) x Gender (2) x IQ-Scale (2) x Time (3) 
# - "Age" is varied between subjects (young, old),
# - "Gender" is varied between subjects (female, male) 
# - "Scale"  is varied within subjets (verbal IQ -VIQ, performance IQ - PIQ); "Scale" is a repeated measure
# - "Time"  is varied within subjets (1, 2,  3); "Time" is a repeated measure. 

# We have several hypotheses about means, partly reproducing some prejudices
#  (1) Men score higher than women on PIQ
#  (2) Women score higher than men on VIQ
#  (3) PIQ decreases more with age than VIQ
#  (4) Perfomance improves over assessments more on 

#             VIQ_1 VIQ_2 VIQ_3 PIQ_1 PIQ_2 PIQ_3
# yng-female   118   121   122   120   123   122
# yng-male     112   115   114   118   120   121
# old-female   114   117   117   104   106   106
# old-male     108   112   112   112   113   112

IQ_means <- matrix(c(118,   121,   122,   120,   123,   122,
			   112,   115,   114,   118,   120,   121,
			   114,   117,   117,   104,   106,   106,
			   108,   112,   112,   112,   113,   112), nrow=4, ncol=6, byrow=TRUE)

# We have two hypotheses about correlations:
#  (1) Correlation are high between 3 VIQ measures and between 3 PIQ measures, but decrease with distance (simplex structure)
#  (2) Correlations are intermediate between VIQ and PIQ, but higher for same time points, again simplex structure 

IQ_corrs <- matrix(c(  1, .8, .7, .5, .4, .3, 
	                .8,  1, .8, .4, .5, .4,
	   	          .7, .8,  1, .3, .4, .5,
	                .5, .4, .3,  1, .8, .7,
                      .4, .5, .4, .8,  1, .8,  
			    .3, .4, .5, .7, .8,  1), nrow=6, ncol=6, byrow=TRUE)

set.seed(1)   # Let's all generate the same data
IQdata.w <- mixedDesign(B=c(2,2), W=c(2,3), M=IQ_means, SD=10, R=IQ_corrs, n=10, long=FALSE)

names(IQdata.w) <- c("Age", "Gender", "Subj",  "VIQ_1", "VIQ_2", "VIQ_3", "PIQ_1", "PIQ_2", "PIQ_3")
levels(IQdata.w$Age) <- c("young", "old")
levels(IQdata.w$Gender) <- c("female", "male")
str(IQdata.w)

# A quick check that the raw data agree with our specification
ddply(IQdata.w, .(Age, Gender), summarise, 
	viq_1=mean(VIQ_1), viq_2=mean(VIQ_2), viq_3=mean(VIQ_3), 
	piq_1=mean(PIQ_1), piq_2=mean(PIQ_2), piq_3=mean(PIQ_3), 
	N=length(VIQ_1), r12=cor(PIQ_1, PIQ_2), r13=cor(PIQ_1, PIQ_3),
	r14=cor(PIQ_1, VIQ_1), r15=cor(PIQ_1, VIQ_2), r16=cor(PIQ_1, VIQ_3) )

# Long format
set.seed(1)   # Let's generate the same data again, but in long format
IQdata <- mixedDesign(B=c(2,2), W=c(2,3), M=IQ_means, SD=10, R=IQ_corrs, n=10, long=TRUE)

names(IQdata) <- c("Age", "Gender", "Subj", "Scale", "Time", "Score")
levels(IQdata$Age) <- c("young", "old")
levels(IQdata$Gender) <- c("female", "male")
levels(IQdata$Scale) <- c("VIQ", "PIQ")
levels(IQdata$Time) <- c("1", "2", "3")

# Check
str(IQdata)
dcast(IQdata, value.var="Score", Age+Gender ~ Scale + Time, mean)

# To go back from long to wide format use dcast() 
IQdata.w.2 <- dcast(IQdata, value.var="Score", Age + Gender + Subj ~ Scale + Time, mean)

# Check
identical(IQdata.w, IQdata.w.2)

table <- ddply(IQdata, .(Age, Gender, Scale, Time), summarise, 
		                                        N=length(Score), M=mean(Score),
		                                        SD=sd(Score), SE=SD/sqrt(N))
	
# Mixed-model ANOVA (use long format!)
summary(aov(Score ~ Age*Gender*Scale*Time + Error(Subj/(Scale*Time)), data=IQdata))
 
# NOT SIGNIFICANT: Age x Gender interaction  -- both factors vary between subjects; no correction necessary
(table.AxG <- ddply(IQdata, .(Age, Gender), summarise, N=length(Score), M=mean(Score), 
 			                                     SE = sd(Score)/sqrt(N))  )
 
 # ... check
source("functions/summarySE.R") # by Winston Chang
(temp <- summarySE(IQdata, measurevar="Score", groupvars=c("Age", "Gender") ))
 
qplot(data = table.AxG, x = Age, y = M, group = Gender, shape = Gender, geom=c("point", "line")) + 
 	geom_point(size=3) +
 	geom_errorbar(aes(ymax = M + 2*SE, ymin = M - 2*SE), width=.05) + theme_bw()

# For computation of standard errors associated with within-subject factors, 
# ... we must remove between-subject differences (norming)
GM <- mean(IQdata$Score)

# ... remove differences between subjects, but keep effects between levels of within-subject factors
IQdata <- ddply(IQdata, .(Subj), transform, Score.w = Score - mean(Score) + GM)  

# ... check 
source("functions/normDataWithin.R")   # by Winston Chang 
IQdata <- normDataWithin(IQdata, idvar="Subj", measurevar="Score", betweenvars=c("Age", "Gender"))


# Main effect of Time and three types of SEs (raw, Cousineau, Morey)
# Morey factor:  sqrt ( n of measures / (n of measures - 1) )
mf <- sqrt(3/2)  # Three levels of Time
(table.T <- ddply(IQdata, .(Time), summarise, 
	          N = length(Score), M=mean(Score), 
	          SE = sd(Score)/sqrt(N),
	          SE_C = sd(Score_norm)/sqrt(N),
	          SE_M = sd(Score_norm)/sqrt(N)*mf,   # 3 levels
	          CI_M = sd(Score_norm)/sqrt(N)*mf*qt(.975, N-1)
))

# ... check
source("functions/summarySEwithin.R")
(temp <- summarySEwithin(IQdata, idvar="Subj",  measurevar="Score", betweenvars=c(), withinvars="Time"))

qplot(data = table.T, x = Time, y = M, group = 1, geom=c("point", "line")) + 
	geom_errorbar(aes(ymax = M + SE, ymin = M - SE), width=.1) +                      # Bad choice
	geom_errorbar(aes(ymax = M + SE_C, ymin = M - SE_C), width=.1, colour="blue") +
	geom_errorbar(aes(ymax = M + SE_M, ymin = M - SE_M), width=.1, colour="red") +
	theme_bw()

# Gender x Scale interaction 
mf <- sqrt(2/1)  # Scale has two levels
(table.GxS <- ddply(IQdata, .(Gender, Scale), summarise, 
			  N = length(Score), M=mean(Score), 
			  SE = sd(Score)/sqrt(N),
			  SE_C = sd(Score_norm)/sqrt(N),
			  SE_M = sd(Score_norm)/sqrt(N)*mf,   # 3 levels
			  CI_M = sd(Score_norm)/sqrt(N)*mf*qt(.975, length(Score)-1)
))

# ... check
(temp <- summarySEwithin(IQdata, idvar="Subj",  measurevar="Score", betweenvars="Gender", withinvars="Scale"))


qplot(data=table.GxS, x = Gender, y = M, group = Scale, shape = Scale, geom=c("point", "line")) + 
	geom_point(size=3) +
	geom_errorbar(aes(ymax = M + SE, ymin = M - SE), width=.1) +                      # Bad choice??
	geom_errorbar(aes(ymax = M + SE_C, ymin = M - SE_C), width=.1, colour="blue") +
	geom_errorbar(aes(ymax = M + SE_M, ymin = M - SE_M), width=.1, colour="red") +
	theme_bw()

# Age x Scale interaction 
mf <- sqrt(2/1)  # Scale has two levels
(table.AxS <- ddply(IQdata, .(Age, Scale), summarise, 
			  N = length(Score), M=mean(Score), 
			  SE = sd(Score)/sqrt(N),
			  SE_C = sd(Score_norm)/sqrt(N),
			  SE_M = sd(Score_norm)/sqrt(N)*mf,   # 3 levels
			  CI_M = sd(Score_norm)/sqrt(N)*mf*qt(.975, length(Score)-1)
))

# ... check
(temp <- summarySEwithin(IQdata, idvar="Subj",  measurevar="Score", betweenvars="Age", withinvars="Scale"))


qplot(data=table.AxS, x = Age, y = M, group = Scale, shape = Scale, geom=c("point", "line")) + 
	geom_point(size=3) +
	geom_errorbar(aes(ymax = M + SE, ymin = M - SE), width=.1) +                      # Bad choice??
	geom_errorbar(aes(ymax = M + SE_C, ymin = M - SE_C), width=.1, colour="blue") +
	geom_errorbar(aes(ymax = M + SE_M, ymin = M - SE_M), width=.1, colour="red") +
	theme_bw()

# NOT SIGNIFICANT: Age x Gender x Scale interaction 
mf <- sqrt(2/1)  # Scale has two levels
(table.AxGxS <- ddply(IQdata, .(Age, Gender, Scale), summarise, 
			    N = length(Score), M=mean(Score), 
			    SE = sd(Score)/sqrt(N),
			    SE_C = sd(Score_norm)/sqrt(N),
			    SE_M = sd(Score_norm)/sqrt(N)*mf,   # 3 levels
			    CI_M = sd(Score_norm)/sqrt(N)*mf*qt(.975, length(Score)-1)
))

# ... check
(temp <- summarySEwithin(IQdata, idvar="Subj",  measurevar="Score", betweenvars=c("Age", "Gender"), withinvars="Scale"))

qplot(data = table.AxGxS, x = Age, y = M,  
	group = Scale, shape = Scale, geom=c("point", "line"),
	facets = . ~ Gender) + geom_point(size=3) +
	geom_errorbar(aes(ymax = M + SE, ymin = M - SE), width=.1) +
	geom_errorbar(aes(ymax = M + SE_C, ymin = M - SE_C), width=.1, colour="blue") +
	geom_errorbar(aes(ymax = M + SE_M, ymin = M - SE_M), width=.1, colour="red") +
	theme_bw()

# We could specify  different correlation matrices for old and young. 
# Dedifferentiation hypothesis: Correlations between measures of different abilities increase with age
# We use IQ_corr for young adults

Ryf <- Rym <- IQ_corrs  # We use IQ_corr for young adults

Rof <- Rom <- matrix(c(  1, .8, .7, .7, .6, .5, 
			      .8,  1, .8, .6, .7, .6,
				.7, .8,  1, .5, .6, .7,
				.7, .6, .5,  1, .8, .7,
				.6, .7, .6, .8,  1, .8,  
				.5, .6, .7, .7, .8,  1), nrow=6, ncol=6, byrow=TRUE)

IQ_corrs2 <- list(Ryf, Rym, Rof, Rom)

set.seed(1)   # Let's generate the same data again, but in long format
IQdata2 <- mixedDesign(B=c(2,2), W=c(2,3), M=IQ_means, SD=10, R=IQ_corrs2, n=10, long=TRUE)
names(IQdata2) <- c("Age", "Gender", "Subj", "Scale", "Time", "Score")
levels(IQdata2$Age) <- c("young", "old")
levels(IQdata2$Gender) <- c("female", "male")
levels(IQdata2$Scale) <- c("VIQ", "PIQ")
levels(IQdata2$Time) <- c("1", "2", "3")

# Does the change in correlations change the ANOVA F-statistics? Compare with above
summary(aov(Score ~ Age*Gender*Scale*Time + Error(Subj/(Scale*Time)), data=IQdata2))

summary(aov(Score ~ Age*Gender*Scale*Time + Error(Subj/(Scale*Time)), data=IQdata))
