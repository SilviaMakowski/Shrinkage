## title: Example for illustration of mixedDesign() function
## author: Reinhold Kliegl

# Citation: Hohenstein & Kliegl (2013): Simulation of Factorial Mixed-Model Designs in R: The mixedDesign() function

library(MASS)
library(plyr)
library(ggplot2)

source("functions/mixedDesign.v0.6.2.R")

# Cell means (between-subject factor levels across rows; within-subject factor levels across columns)
mean.mat <- matrix(c(950, 950,  950, 1050, 1050, 1050, 
			   850, 800,  750,  950,  900,  850,
			   800, 700,  600,  900,  800,  700,
			   830, 730,  630,  930,  830,  730), nrow=4, ncol=6, byrow=TRUE)

# Call function 
set.seed(1)
data <- mixedDesign(B = 4, W = c(2, 3), M = mean.mat, SD = 60, n = 20, long = TRUE)

# Rename variables and levels
names(data) <- c("Age", "Subj", "Prime", "Frequency", "RT")
levels(data$Age) <- c("Children", "Teenager", "Young \nadults", "Old \nadults")
levels(data$Prime) <- c("Related", "Unrelated")
levels(data$Frequency) <- c("Low", "Medium", "High")

# Compute table of means for full factorial
table <- ddply(data, .(Age, Frequency, Prime), summarise, N=length(RT), M=mean(RT),
		                                              SD=sd(RT), SE=SD/sqrt(N) )
# Note: SE's are valid only for comparisons between age groups (between-subject factor)

# Plot table
qplot(data=table, x=Age, y=M, ylab = "Response time (ms)", group=Frequency, colour=Frequency, 
	facets = . ~ Prime, geom=c("point", "line")) +
	geom_errorbar(aes(ymax=M+2*SE, ymin=M-2*SE), width=0) + theme_bw()

# Age(4) x Prime (2) x Freq(3) mixed-model ANOVA
summary(aov(RT ~ Age*Frequency*Prime + Error(Subj/(Frequency*Prime)), data=data))


# Significant interaction for Age x Frequency
table.a_f <- ddply(data, .(Age, Frequency), summarise, N=length(RT), M=mean(RT), 
			                                     SD=sd(RT), SE=SD/sqrt(N) )

qplot(data=table.a_f, x=Age, y=M, ylab = "Response time (ms)", 
	group=Frequency, colour=Frequency, shape=Frequency,
	geom=c("point", "line")) +
	geom_errorbar(aes(ymax=M+2*SE, ymin=M-2*SE), width=0) + theme_bw()

