# A few quick demos
# F1-F2-ANOVAs from trial x trial file
# How to plot a main effect with a line in qplot()
# 5 July 2013, r. kliegl

# Updated after comment by Ómar I. Jóhannesson
# 6 July 2013, r. kliegl

library(languageR)  # for lexdec data
library(gtools)     # for quantcut()
library(MASS)       # for boxcox()
library(plyr)       # for ddply(), arrange()
library(ggplot2)    # for qplot()
library(lme4)       # for lmer()

rm(list=ls())

# ------------------------------------------------

# Let's generate a 3-level factor frequency; it will be within-subject, but between-item
lexdec$FRQ <- quantcut(lexdec$Frequency, q=seq(0, 1, 1/3), labels = c("low", "med", "high"))

# RTs are in log format, but should be converted to speed  
boxcox(exp(RT) ~ Subject, data=lexdec)
lexdec <- transform(lexdec, speed=1000/exp(RT) )

# ------------------------------------------------

# Aggregate across words within NativeLanguage and FRQ; Subject is random factor in F1-ANOVA
df.subject <- ddply(lexdec, .(Subject, NativeLanguage, FRQ), summarise, N=length(speed), SPEED = mean(speed))

# full factorial table for plot
(table.F1 <-ddply(df.subject, .(NativeLanguage, FRQ), summarise, 
			     N=length(SPEED), M=mean(SPEED), SD=sd(SPEED),
			     SE=SD/sqrt(N)) )  # valid only for between-subject factor

# F1-ANOVA 
# df.subject, not lexdec, is the data frame we need for the F1-ANOVA
# NativeLanguage (between-subject), FRQ (within-subject), Subject must be of type factor in df.subject
str(df.subject)
summary(aov(SPEED ~ NativeLanguage*FRQ + Error(Subject/FRQ), data=df.subject))

qplot(data=table.F1, x=FRQ, y=M, group=NativeLanguage, colour=NativeLanguage, geom=c("point", "line")) +
	geom_errorbar(aes(ymax=M+SE, ymin=M-SE), width=0.02) + ylab("Response speed") + theme_bw()

# ------------------------------------------------

# Aggregate across subjects within NativeLanguage and FRQ; Word is random factor in F2-ANOVA
df.word <- ddply(lexdec, .(Word, NativeLanguage, FRQ), summarise, N=length(speed), SPEED = mean(speed))

# F2-ANOVA
# df.word, not lexdec, is the data frame we need for the F2-ANOVA
# NativeLanguage (within-item), FRQ (between-item), Item must be of type factor in df.word
summary(aov(SPEED ~ NativeLanguage*FRQ + Error(Word/NativeLanguage), data=df.word))

# table for main effect of FRQ
# full factorial table for plot
(table.F2 <-ddply(df.word, .(FRQ), summarize, 
			N=length(SPEED), M=mean(SPEED), SD=sd(SPEED),
			SE=SD/sqrt(N)) )  # valid only for between-item factor))

qplot(data=table.F2, x=FRQ, y=M, group=1,  geom=c("point", "line")) +
	geom_errorbar(aes(ymax=M+SE, ymin=M-SE), width=0.02) + ylab("Response speed") + theme_bw()


# ------------------------------------------------


# Why you need Subject to be a factor or character, but not numeric, for aov()?
# ... convert factor level to number and letter
subid <- data.frame(id=1:length(levels(df.subject$Subject)), 
			  idc=letters[1:length(levels(df.subject$Subject))],
			  Subject=levels(df.subject$Subject))
str(subid)

df.subject <- merge(df.subject, subid)
df.subject$idc <- as.character(df.subject$idc) 
str(df.subject)

# ... subject is of type factor--THAT IS A CORRECT RESULT
summary(aov(SPEED ~ NativeLanguage*FRQ + Error(Subject/FRQ), data=df.subject))

# ... subject is of type character--THAT IS A CORRECT RESULT, TOO 
summary(aov(SPEED ~ NativeLanguage*FRQ + Error(idc/FRQ), data=df.subject))

# ... subject is of type integer--WRONG RESULT!
summary(aov(SPEED ~ NativeLanguage*FRQ + Error(id/FRQ), data=df.subject))

# Why is this wrong? In mixed-model ANOVAs every within-subject factor and 
# every interaction between within-subject factors have their own error term.
# The output of the last anova shows that there is only one error term for
# all of the sources of variance. 

# Demonstration of the effect of non-unique subject codes
df <- data.frame(expand.grid(id=c("S1", "S2", "S3", "S4"), B_A=c("A1", "A2"), W_a=c("a1", "a2", "a3")))
df <- cbind(df, Subj=interaction(df$B_A, df$id, sep="_"), dv=rnorm(nrow(df), 100, 10))
df <- arrange(df, B_A, Subj)
str(df)

# Correct subject
summary(aov(dv ~ B_A*W_a + Error(Subj/W_a), data=df))

# "Wrong" subject
summary(aov(dv ~ B_A*W_a + Error(id/W_a), data=df))

