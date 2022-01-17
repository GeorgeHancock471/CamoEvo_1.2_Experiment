# Contents

#[0] SETUP AND IMPORT DATA

#[1] STATISTICAL ANALYSES
#   Data Distribution
#   Selection Correlation
#   Fitness Correlation


#[2] FANCY PLOTS
#    Selection Correlation
#    Fitness Correlation

#[3] PCA ANALYSES
#    Data Setup
#    Plots

#[4] SURVIVAL TIME

# -----------------------------------------------------------------------------------------------------------------
#  (0) SETUP AND IMPORT DATA
# -----------------------------------------------------------------------------------------------------------------

# Datasets from CamoEvo have been modified to include calculations of a second measure of contrast difference.
# Contrast_Difference2 = euclidean difference BgLocal_Contast/BgLocal_Luminance - Target_Contast/Target_Luminance.
# Additionally the order, internal order and status (survived/dead) of each individual has been added.

# LAB measures of mean, contrast and GabRat have been re-added posthoc. This was due to at the time of
# The experiment CamoEvo did not output the SD or GabRat of a* or b*, this has since been changed.


# This was done using excel.



#Get Packages
#=============================
library(ggplot2)
library(ggimage)
library(scales)
library(Matrix)
library(lme4)
library(nlme)
library(lmerTest)
library(pbkrtest)
library(car)
library(olsrr)
library(lattice)

sessionInfo()

#Get Data
#=============================
A1<-read.table("Data_Output_treatment_A1.txt", header=TRUE)
A2<-read.table("Data_Output_treatment_A2.txt", header=TRUE)
A3<-read.table("Data_Output_treatment_A3.txt", header=TRUE)

B1<-read.table("Data_Output_treatment_B1.txt", header=TRUE)
B2<-read.table("Data_Output_treatment_B2.txt", header=TRUE)
B3<-read.table("Data_Output_treatment_B3.txt", header=TRUE)

C1<-read.table("Data_Output_treatment_C1.txt", header=TRUE)
C2<-read.table("Data_Output_treatment_C2.txt", header=TRUE)
C3<-read.table("Data_Output_treatment_C3.txt", header=TRUE)

#Add Background
#==========================
A1$background <- "scrub"
A2$background <- "scrub"
A3$background <- "scrub"

B1$background <- "leaf"
B2$background <- "leaf"
B3$background <- "leaf"

C1$background <- "veg"
C2$background <- "veg"
C3$background <- "veg"

#Add Population
#==========================
A1$population <- "1"
A2$population <- "2"
A3$population <- "3"

B1$population <- "1"
B2$population <- "2"
B3$population <- "3"

C1$population <- "1"
C2$population <- "2"
C3$population <- "3"

#Add Treatment
#==========================
A1$treatment <- "Scrub1"
A2$treatment <- "Scrub2"
A3$treatment <- "Scrub3"

B1$treatment <- "Leaf1"
B2$treatment <- "Leaf2"
B3$treatment <- "Leaf3"

C1$treatment <- "Veg1"
C2$treatment <- "Veg2"
C3$treatment <- "Veg3"


#Merge Data
#==========================
Data <- rbind(A1,A2,A3,B1,B2,B3,C1,C2, C3)



# -----------------------------------------------------------------------------------------------------------------
#  (1) STATISTICAL ANALYSES
# -----------------------------------------------------------------------------------------------------------------

# [0] Data Distribution
#========================================================================

#Time
#........
hist(log(AllData$Survival_Time), breaks=10)

shapiro.test(log(AllData$Survival_Time))

#Significant (2.2e-16), 95% confidence interval



#Luminance
#........
hist(log(AllData$Luminance_Difference), breaks=10)

shapiro.test(log(AllData$Luminance_Difference))

#Significant ( 2.2e-16), 95% confidence interval


#Contrast
#........
hist(log(AllData$Contrast2_Difference), breaks=10)

shapiro.test(log(AllData$Contrast2_Difference))

#Significant different from normal (2.2e-16), 95% confidence interval


#Colour
#........
hist(log(AllData$Colour_Difference), breaks=10)
shapiro.test(log(AllData$Colour_Difference))

#Significant (2.2e-16), 95% confidence interval


#GabRat
#........
hist(log(AllData$GabRat_Edge_Disruption), breaks=10)
shapiro.test(log(AllData$GabRat_Edge_Disruption))

#Significant (2.2e-16), 95% confidence interval



# [1] Linear Models of Selection
#========================================================================

AllData$uniquePop <- paste(AllData$treatment,AllData$Generation,sep="_")

AllData$uniquePop

AllData




# Fitness
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,log(Fitness), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

# Positive Correlation


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(log(Fitness) ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significant increase p<0.0001
# significant polynomial for scrub (flattening)
# significantly lower estimate for vegetation (not increasing as fast)


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')






# L Mean Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,LocalDifMeanL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#  Negative Correlation


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(LocalDifMeanL ~ poly(Generation,2) * background   + (1|treatment), data = AllData)


# significant decrease p<0.0001
# significant polynomial for scrub (flattening)
# significantly lower estimate for vegetation (not increasing as fast)


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')




# A Mean Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,LocalDifMeanA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#  Negative Correlation


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(LocalDifMeanA ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significant decrease p<0.0001



#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')



# B Mean Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,LocalDifMeanB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#  Negative Correlation, but not for scrubland


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(LocalDifMeanB ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significant decrease p<0.0001
# both scrubland and vegetation have significantly higher estimates (slower decrease) then leaflitter


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')






# L Contrast Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,LocalDifContrastL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#  No correlation except leaflitter


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(LocalDifContrastL ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significant decrease in leaflitter p<0.0001
# no significant decrease for scrub or leaflitter


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')






# A Contrast Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,LocalDifContrastA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#  No correlation except scrubland


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(LocalDifContrastA ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significant decrease in leaflitter p<0.0001
# significant polynomial for scrub (flattening)
# significantly lower estimate for vegetation (not increasing as fast)


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')




# B Contrast Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,LocalDifContrastB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#  No correlation except scrubland


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(LocalDifContrastB ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significant decrease in leaflitter p<0.0001
# significant polynomial for scrub (flattening)
# significantly lower estimate for vegetation (not increasing as fast)


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')


# GabRatL
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,gabRatL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

# Increases and faster for scrubland


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(gabRatL ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significantly increases p<0.0001
# increases significantly faster for scrub and is flattening


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')


# GabRatA
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,gabRatA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

# Increases for leaflitter and scrub but decreases for vegetation


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(gabRatA ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significantly increases p<0.0001
# does not increase for vegetation

#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')



# GabRatB
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes(Generation,gabRatB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

# Increases for leaflitter but not for the others


#Model
#-------------------

#treatment = Participant
#Background = the the images

m1 <- lmer(gabRatB ~ poly(Generation,2) * background   + (1|treatment), data = AllData)

# Then look at the model summary:
anova(m1)
summary(m1)

# significantly increases for leaflitter only p<0.0001

#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')


# [2] Linear Models of Survival Time
#========================================================================



# Mean Luminance Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( LocalDifMeanL,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#Survival Time appears to decrease with luminance difference.


#LMER Model
#-------------------

#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ Luminance_Difference * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)

# Significant negative correlation p<0.0001
# Vegetation decreases less quickly but not significantly


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')







# Mean A Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( LocalDifMeanA,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#Survival Time appears to decrease with A difference.

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ LocalDifMeanA * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)


# Significant negative correlation p<0.0001

# Vegetation does not have a negative correlation


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')





# Mean B Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( LocalDifMeanB,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#Survival Time appears to decrease with B difference but not for scrubland.

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ LocalDifMeanB * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)


# Significant negative correlation p<0.0001

# Scrubland does not have a significant interaction.


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')



# L Contrast Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( LocalDifContrastL,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#No Obvious relationship.

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ LocalDifContrastL * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)


# Significant negative correlation p<0.01

# Vegetation has a positive correlation p<0.01


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')






# A Contrast Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( LocalDifContrastA,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#No Obvious relationship.

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ LocalDifContrastA * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)


# Significant negative correlation p<0.001

# Vegetation has no correlation p<0.01


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')






# B Contrast Difference
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( LocalDifContrastB,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#Positive correlation for scrub and veg, negative correlation for leaflitter

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ LocalDifContrastB * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)


summary(m1)

ranef(m1)


# No significant correlation, except scurbland which significantly increases

# Vegetation has no correlation p<0.01


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')





# GabRat L
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( gabRatL ,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#Positive correlation 

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ gabRatL  * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)

# Significant positive correlation p<0.001



#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')




# GabRat A
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( gabRatA ,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#Positive correlation for scrub, no obious correlation for the others.

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ gabRatA  * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)

# No significant correlation except for Scrub, p<0.01


#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')


# GabRat B
#==============================

#Quick Plot
#-------------------

p1 <- ggplot(AllData, aes( gabRatB ,log(Survival_Time), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)
p1

#No obvious correlation

head(AllData)


#LMER Model
#-------------------
#treatment = Participant
#Int_Order = Slide Order for each Generation (0-23)
#Background = the Background Images (scrub,vegetation & leaflitter)

m1 <- lmer(log(Survival_Time) ~ gabRatB  * background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)
m2 <- lmer(log(Survival_Time) ~ background + (1|Generation) + (1|treatment) + (1|Int_Order), data = AllData)

# Then look at the model summary:
anova(m1,m2)

summary(m1)
anova(m1)
rand(m1)
ranef(m1)


# Significant but weak increase 



#Check Residuals
#-------------------

#Linearity

plot(resid(m1))

#Normality

qqmath(m1)

library(effects)
library(sjPlot)

plot_model(m1, type='diag')




# -----------------------------------------------------------------------------------------------------------------
#  (3) FANCY PLOTS
# -----------------------------------------------------------------------------------------------------------------


# HEX values=c("#cc9900","#997300","#664d00", "#c43b3b", "#9d2f2f", "#762323",  " #00cc00"," #009900","#006600"   )) +  

cbPalette <- c( "#56B4E9", "#E69F00","#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

greyPalette <- c( "#000000", "#4d4d4d", "#999999")



  
# [0] Selection (Treatment Measures)
#========================================================================================================
  
  
  # Fitness
  #........................................
  
  p1 <- ggplot(AllData, aes(Generation, (Fitness)/1000, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)+
    geom_hline(yintercept=15, 
               color = "red", size=0.85, alpha=0.5)
  
  
  p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
    scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
    scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
    theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
    xlab("Generation") + ylab("Survival Time (seconds)")+
    scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
    scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
    scale_y_log10()+
    coord_cartesian( ylim = c(0.300, 15))+
    theme(
      axis.text = element_text(family = "sans", size = 22),
      axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
      axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
      axis.line = element_line(colour = "black", size=0),
      panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")
  
  
  
# L Mean Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, LocalDifMeanL, col = background, fill = background)) +
    geom_point(shape =21,  alpha = 0.1, size=2)  +
    geom_smooth(method=glm,formula=y~poly(x,2),
                se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("L Mean Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 50, by = 10))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# a* mean Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, LocalDifMeanA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,1),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("A Mean Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 70, by = 10))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# b* mean Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, LocalDifMeanB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,1),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("B Mean Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 70, by = 10))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")


# L* Contrast Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, LocalDifContrastL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("L Contrast Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 30, by = 5))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# A* Contrast Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, LocalDifContrastA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("A Contrast Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 30, by = 5))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# B* Contrast Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, LocalDifContrastB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("B Contrast Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 30, by = 5))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# GabRat L
#........................................

p1 <- ggplot(AllData, aes(Generation, gabRatL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("L Edge Disruption")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# GabRat A
#........................................

p1 <- ggplot(AllData, aes(Generation, gabRatA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("A Edge Disruption")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# GabRat B
#........................................

p1 <- ggplot(AllData, aes(Generation, gabRatB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=2)  +
  geom_smooth(method=glm,formula=y~poly(x,1),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("B Edge Disruption")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))+ 
  theme(
    axis.text = element_text(family = "sans", size = 22),
    axis.title.x = element_text(family = "sans", size = 26, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 26, margin=margin(0,5,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# [1] Fitness Correlation
#========================================================================================================

# L* Mean Difference
#........................................

p1 <- ggplot(AllData, aes( LocalDifMeanL,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")


# A* Mean Difference
#........................................

p1 <- ggplot(AllData, aes( LocalDifMeanA,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# B* Mean Difference
#........................................

p1 <- ggplot(AllData, aes( LocalDifMeanB,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# L* Contrast Difference
#........................................

p1 <- ggplot(AllData, aes( LocalDifContrastL,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# A* Contrast Difference
#........................................

p1 <- ggplot(AllData, aes( LocalDifContrastA,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# B* Contrast Difference
#........................................


p1 <- ggplot(AllData, aes( LocalDifContrastB,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")







# GabRat L
#........................................


p1 <- ggplot(AllData, aes( gabRatL,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# GabRat A
#........................................

p1 <- ggplot(AllData, aes( gabRatA,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# GabRat B
#........................................

p1 <- ggplot(AllData, aes( gabRatB,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")









# [2] Selection (Target Measures)
#========================================================================================================


# TargetL
#........................................

p1 <- ggplot(AllData, aes(Generation, TmeanL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Target L")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")





# TargetA
#........................................

p1 <- ggplot(AllData, aes(Generation, TmeanA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Target A")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-60, 60, by = 10))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")





# TargetB
#........................................

p1 <- ggplot(AllData, aes(Generation, TmeanB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Target B")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-10, 80, by = 10))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Target Contrast L
#........................................

p1 <- ggplot(AllData, aes(Generation, TcontrastL, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Target Contrast b")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-10, 80, by = 10))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# Target Contrast A
#........................................

p1 <- ggplot(AllData, aes(Generation, TcontrastA, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Target Contrast b")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-10, 80, by = 10))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# Target Contrast B
#........................................

p1 <- ggplot(AllData, aes(Generation, TcontrastB, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Target Contrast b")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-10, 80, by = 10))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")















# [3] Fitness Correlation, Bg and Target Measures
#========================================================================================================





# Local Contrast L
#........................................


p1 <- ggplot(AllData, aes( LcontrastL,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Local Contrast A
#........................................


p1 <- ggplot(AllData, aes( LcontrastA,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Local Contrast B
#........................................


p1 <- ggplot(AllData, aes( LcontrastB,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Target Mean L
#........................................


p1 <- ggplot(AllData, aes( TmeanL,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Target Mean A
#........................................


p1 <- ggplot(AllData, aes( TmeanA,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")





# Target Mean B
#........................................


p1 <- ggplot(AllData, aes( TmeanB,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")







# Target Contrast L
#........................................


p1 <- ggplot(AllData, aes( TcontrastL,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")



# Target Contrast A
#........................................


p1 <- ggplot(AllData, aes( TcontrastA,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Target Contrast B
#........................................


p1 <- ggplot(AllData, aes( TcontrastB,Fitness, col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")






# [3] Miscellaneous
#========================================================================================================

#Fitness & Order

p1 <- ggplot(AllData, aes(Int_Order,(Fitness))) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method = "lm", se = TRUE, alpha=0.1)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))    


#Fitness & Status
p1 <- ggplot(AllData, aes(Generation,(Fitness/1000), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE, alpha=0.1,aes(linetype=Status))

p1 + theme_classic()  +
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 15, by = 1.50))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))  


#Luminance & Status
p1 <- ggplot(AllData, aes(Generation,(Luminance_Difference), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE, alpha=0.1,aes(linetype=Status))


p1 + theme_classic()  +
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 70, by = 5))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))  

#Colour & Status
p1 <- ggplot(AllData, aes(Generation,(Colour_Difference), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE, alpha=0.1,aes(linetype=Status))


p1 + theme_classic()  +
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 70, by = 5))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))  




#GabRat & Status
p1 <- ggplot(AllData, aes(Generation,(GabRat_Edge_Disruption), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE, alpha=0.1,aes(linetype=Status))


p1 + theme_classic()  +
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.05))+ 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))  


#-------------------------------------------------------------------------------------------------------------------


head(demoData3)


p1 <- ggplot(demoData2, aes((Survival_Time), Colour_Difference, col = background, fill = background)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  + 
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 


p1 <- ggplot(demoData2, aes((Survival_Time), Luminance_Difference, col = background, fill = background)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  + 
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 


p1 <- ggplot(demoData2, aes((Survival_Time),GabRat_Edge_Disruption, col = background, fill = background)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  + 
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 



p1 <- ggplot(demoData2, aes((Survival_Time), PC1, col = background, fill = background)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  + 
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 





p1 <- ggplot(demoData2, aes(Colour_Difference, Survival_Time, col = background, fill = background)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  + 
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 


p1 <- ggplot(demoData2, aes(x=Generation, y= Colour_Difference, color=background)) + geom_point(shape=1) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE)

p1 + theme_light()  + theme(text = element_text(size = 16)) +  xlab("Generation")  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)   


p1 <- ggplot(demoData2, aes(x=Generation, y= Luminance_Difference, color=background)) + geom_point(shape=1) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE)

p1 + theme_light()  + theme(text = element_text(size = 16)) +  xlab("Generation")  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)   



p1 <- ggplot(demoData2, aes(x=Generation, y= GabRat_Edge_Disruption, color=background)) + geom_point(shape=1) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE)

p1 + theme_light()  + theme(text = element_text(size = 16)) +  xlab("Generation")  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)   



p1 <- ggplot(demoData2, aes(x=GabRat_Edge_Disruption, y= (Capture_Time), color=background)) + geom_point(shape=1) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE)

p1 + theme_light()  + theme(text = element_text(size = 16)) +  xlab("GabRatL")  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)






p1 <- ggplot(demoData2, aes(x=Colour_Difference, y= (Capture_Time), color=background)) + geom_point(shape=1) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE)

p1 + theme_light()  + theme(text = element_text(size = 16))  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)


p1 <- ggplot(demoData2, aes(x=Luminance_Difference, y= (Capture_Time), color=background)) + geom_point(shape=1) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              method.args=list(family="quasipoisson")
              ,  se=TRUE)

p1 + theme_light()  + theme(text = element_text(size = 16))  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)



p1 <- ggplot(C1, aes(x=Target_A, y= Target_B, color=Generation, fill=Generation)) +
  geom_point(shape =21, color="black", alpha = 0.25, size=4)

p1 + theme_light()  + theme(text = element_text(size = 16))  + theme(
  axis.title.x = element_text(family = "sans", size = 15, margin=margin(5,0,0,0)), 
  axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,20,0,0)), 
)










p1 <- ggplot(B3, aes( Target_A, Target_B, col = Generation, fill = Generation)) +
  geom_point(shape =21, col="black",  alpha = 0.25, size=5)

p1 <- p1 + theme_classic()  +  xlim(-60, 60) +  ylim(-60, 60) +
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + theme(legend.position = "none")

p1


p1 <- ggplot(B3, aes( BgLocal_A, BgLocal_B, col = Generation, fill = Generation)) +
  geom_point(shape =21, col="black",  alpha = 0.25, size=5)

p1 <- p1 + theme_classic()  +  xlim(-60, 60) +  ylim(-60, 60) +
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + theme(legend.position = "none")

p1


p1 <- ggplot(C1, aes( BgLocal_A, BgLocal_B, col = "red", fill = "red")) +
  geom_point(shape =21, col="black",  alpha = 0.00, size=5)+
  stat_ellipse(geom="polygon", col = "black", alpha = 1)

p1 <- p1+ theme_classic()  +  xlim(-60, 60) +  ylim(-60, 60) +
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + theme(legend.position = "none")

p1




p1 <- ggplot(C3, aes( Luminance_Difference, Colour_Difference, col = Generation, fill = Generation)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  +  xlim(0, 60) +  ylim(0, 60) +
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 



p1 <- ggplot(demoData2, aes(Generation, Colour_Difference, col = Generation, fill = Generation)) +
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  +  xlim(0, 15) +  ylim(0, 100) +
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)), 
    panel.border = element_rect(colour = "black", fill=NA, size=1.5)) 





#========================================================================================================
# Selection (Background)
#========================================================================================================


# Luminance Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, NorLum, col = background, fill = background)) +
  geom_smooth(method=glm,formula=y~poly(x,2),
               se=TRUE)


p1 + theme_classic() +
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Luminance Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.2))+ 
  coord_cartesian(xlim =c(0, 15), ylim = c(-1.0, 1.0))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")


# Colour Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, NorCol, col = background, fill = background)) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Colour Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.2))+ 
  coord_cartesian(xlim =c(0, 15), ylim = c(-1.0, 1.0))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")


# Contrast Difference
#........................................

p1 <- ggplot(AllData, aes(Generation, NorCon, col = background, fill = background)) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Contrast Difference")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.2))+ 
  coord_cartesian(xlim =c(0, 15), ylim = c(-1.0, 1.0))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")





# GabRat 
#........................................

p1 <- ggplot(AllData, aes(Generation, NorGab, col = background, fill = background)) +
  geom_smooth(method=glm,formula=y~poly(x,2),
              se=TRUE)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  xlab("Generation") + ylab("Edge Disruption (GabRatL)")+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.2))+ 
  coord_cartesian(xlim =c(0, 15), ylim = c(-1.0, 1.0))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,20,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")


# Fitness
#........................................

p1 <- ggplot(AllData, aes(Generation, NorFit, col = background, fill = background)) +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=1.25)+
  geom_hline(yintercept=15000, 
             color = "red", size=0.5, alpha=0.5)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  xlab("Generation") + ylab("Survival Time (milliseconds)")+
  scale_x_continuous(breaks = seq(0, 15, by = 3))+ 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.2))+ 
  coord_cartesian(xlim =c(0, 15), ylim = c(-1.0, 1.0))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")






#========================================================================================================
# Fitness Correlations (Background)
#========================================================================================================

# Luminance Difference
#........................................

p1 <- ggplot(noTimeOutData, aes( Luminance_Difference,(Fitness), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method = "lm", alpha = 0.1, size=0)+
  stat_smooth(method = "lm",se=FALSE,  alpha = 0.5, size=1.25)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# Colour Difference
#........................................

p1 <- ggplot(noTimeOutData, aes( Colour_Difference,(Fitness), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method = "lm",  alpha = 0.1, size=0)+
  stat_smooth(method = "lm",  se=FALSE,  alpha = 0.5, size=1.25)




p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 70, by = 10))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")





# Contrast Difference
#........................................

p1 <- ggplot(noTimeOutData, aes( Contrast2_Difference,(Fitness), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method = "lm", alpha = 0.1, size=0)+
  stat_smooth(method = "lm", se=FALSE,  alpha = 0.5, size=1.25)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 2, by = 0.5))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")


# GabRat 
#........................................

p1 <- ggplot(noTimeOutData, aes( GabRat_Edge_Disruption,(Fitness), col = background, fill = background)) +
  geom_point(shape =21,  alpha = 0.1, size=1)  +
  geom_smooth(method = "lm",  alpha = 0.1, size=0)+
  stat_smooth(method = "lm", se=FALSE,  alpha = 0.5, size=1.25)


p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1))+ 
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 15000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 14, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,20,0,0)),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  + theme(legend.position="none")




# -----------------------------------------------------------------------------------------------------------------
#  (3) PCA ANALYSES
# -----------------------------------------------------------------------------------------------------------------


# [0] Data Setup
#======================================================================

#PC All Generations
#---------------------------------
geneData = Data[35:71]
colData = cbind(Gen15Data[72:78])


pcData = geneData

#prcomp
myPr <- prcomp(pcData, scale = TRUE, center = TRUE) # this is what runs the PCR, using the columns with genes

myPr

head(myPr) 
summary(myPr) 



#Get correlations with original variables for PC1 and PC2


loadings <- myPr$rotation

loadings

PC1 = loadings[,1]
PC1

PC2 = loadings[,2]
PC2

PCs = loadings[,1:2]

PC1Rank = sqrt(PC1^2)
PC1Rank
PC1Ranked <- PC1Rank[,order(PC1)] 
PC1Ranked

PC2Rank = sqrt(PC2^2)
PC2Rank


PCsRank = sqrt(PCs^2)
PCsRank 

PCsRank<- as.data.frame(PCsRank)

PCsRank$PC1R <- as.numeric(PCsRank$PC1)
PCsRank$PC2R <- as.numeric(PCsRank$PC2)

PCsRank = PCsRank[,3:4]
PCsRank <- rbind(PCs, PCsRank)

# Get Most Important for PC1
PCsRanked1 <- PCsRank[order(PCsRank$PC1R, decreasing = TRUE),]  
PCsRanked1

# Get Most Important for PC2
PCsRanked2 <- PCsRank[order(PCsRank$PC2R, decreasing = TRUE),]  
PCsRanked2

AllData = cbind(Data, myPr$x[,1:2])  # This will create the new dataset which

#Calculate Difference in Local Mean
AllData$LocalDifMeanL <- sqrt((AllData$TmeanL - AllData$LmeanL)^2)
AllData$LocalDifMeanA <- sqrt((AllData$TmeanA - AllData$LmeanA)^2)
AllData$LocalDifMeanB <- sqrt((AllData$TmeanB - AllData$LmeanB)^2)

#Calculate Difference in Local Contrast
AllData$LocalDifContrastL <- sqrt((AllData$TcontrastL - AllData$LcontrastL)^2)
AllData$LocalDifContrastA <- sqrt((AllData$TcontrastA - AllData$LcontrastA)^2)
AllData$LocalDifContrastB <- sqrt((AllData$TcontrastB - AllData$LcontrastB)^2)


#Calculate Difference in Global Mean
AllData$GlobalDifMeanL <- sqrt((AllData$TmeanL - AllData$BmeanL)^2)
AllData$GlobalDifMeanA <- sqrt((AllData$TmeanA - AllData$BmeanA)^2)
AllData$GlobalDifMeanB <- sqrt((AllData$TmeanB - AllData$BmeanB)^2)

#Calculate Difference in Global Contrast
AllData$GlobalDifContrastL <- sqrt((AllData$TcontrastL - AllData$BcontrastL)^2)
AllData$GlobalDifContrastA <- sqrt((AllData$TcontrastA - AllData$BcontrastA)^2)
AllData$GlobalDifContrastB <- sqrt((AllData$TcontrastB - AllData$BcontrastB)^2)




#Calculate Difference in CoV
AllData$LocalDifCoVL <- sqrt((AllData$TcontrastL/AllData$TmeanL - AllData$LcontrastL/AllData$LmeanL)^2)
AllData$LocalDifCoVA <- sqrt((AllData$TcontrastA/AllData$TmeanA - AllData$LcontrastA/AllData$LmeanA)^2)
AllData$LocalDifCoVB <- sqrt((AllData$TcontrastB/AllData$TmeanB - AllData$LcontrastB/AllData$LmeanB)^2)

#PC Gen 15
#---------------------------------
Gen15Data <- subset(Data, data$Generation==15)

geneData = Gen15Data[33:69]
colData =   cbind(Gen15Data[70:76])

pcData = colData 




#prcomp
myPr <- prcomp(pcData, scale = TRUE, center = TRUE) # this is what runs the PCR, using the columns with genes

myPr

head(myPr) 
summary(myPr) 



#Get correlations with original variables for PC1 and PC2


loadings <- myPr$rotation

loadings

PC1 = loadings[,1]
PC1

PC2 = loadings[,2]
PC2

PCs = loadings[,1:2]

PC1Rank = sqrt(PC1^2)
PC1Rank
PC1Ranked <- PC1Rank[,order(PC1)] 
PC1Ranked

PC2Rank = sqrt(PC2^2)
PC2Rank


PCsRank = sqrt(PCs^2)
PCsRank 

PCsRank<- as.data.frame(PCsRank)

PCsRank$PC1R <- as.numeric(PCsRank$PC1)
PCsRank$PC2R <- as.numeric(PCsRank$PC2)

PCsRank = PCsRank[,3:4]
PCsRank <- rbind(PCs, PCsRank)

# Get Most Important for PC1
PCsRanked1 <- PCsRank[order(PCsRank$PC1R, decreasing = TRUE),]  
PCsRanked1
#Pc1 ptn_dim_wdt,col_bot_lmv,ptn_grd_cvr,spk_two_lvl,eem_sig_lv,col_mac_rgv


# Get Most Important for PC2
PCsRanked2 <- PCsRank[order(PCsRank$PC2R, decreasing = TRUE),]  
PCsRanked2
#Pc2 spk_one_rto,col_top_rgv, col_bot_rgv

Gen15PC = cbind(Gen15Data, myPr$x[,1:2])  # This will create the new dataset which







# [1] PCA  Analysis Plots
#======================================================================

# All Generations Data
#........................................
AllData0 <- subset(AllData, data$Generation==0)
AllData5 <- subset(AllData, data$Generation==5)
AllData10 <- subset(AllData, data$Generation==10)
AllData15 <- subset(AllData, data$Generation==15)


AllData0_15 <- rbind(AllData0, AllData15)
AllData0_15$Generation <- as.factor(AllData0_15$Generation )


ScrubData <-subset(AllData, data$background=="scrub")
LeafData <-subset(AllData, data$background=="lea")
VegData <-subset(AllData, data$background=="veg")

noTimeOutData <- AllData[!(AllData$Fitness>= 15000),]


AllDataBest <-subset(AllData, AllData$Status=="Survived")
AllDataBestLast <- rbind(subset(AllDataBest, AllDataBest$Generation==15),subset(AllDataBest, AllDataBest$Generation==14),subset(AllDataBest, AllDataBest$Generation==13)  )
head(A)


#Start
#-------------
p1 <- ggplot( AllDataBest , aes(PC1, PC2, col = background, fill =  background)) +
  stat_ellipse(geom="polygon", col = "black", alpha = 0.25,aes(linetype=population)) + 
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)


p1 + theme_classic()  + scale_x_continuous(breaks = seq(-8, 8, by = 2))+ 
  scale_y_continuous(breaks = seq(-8, 8, by = 2))+ 
  ylim(-8, 8)  +  xlim(-8, 8)  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))   



#Final
#-------------
p1 <- ggplot( AllData15 , aes(PC1, PC2, col = background, fill =  background)) +
  stat_ellipse(geom="polygon", col = "black", alpha = 0.25,aes(linetype=population)) + 
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()  +  ylim(-5, 5)  +  xlim(-5, 5)   + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))   





# Generation 15
#........................................

p1 <- ggplot( Gen15PC, aes(PC1, PC2, col = background, fill =  background)) +
  stat_ellipse(geom="polygon", col = "black", alpha = 0.25,aes(linetype=population)) + 
  geom_point(shape =21, col = "black", alpha = 0.25, size=5)

p1 + theme_classic()   + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))   




# [2] Phenotype Comparison
#======================================================================


#Target Colour
#---------------------
p1 <- ggplot( AllData15, aes(TmeanA, TmeanB, col = background, fill =  background)) +
  stat_ellipse(geom="polygon", col = "black", alpha = 0.25,) + 
  geom_point(shape =22,  alpha = 0.15, size=3)+
  geom_hline(yintercept=0, 
             color = "black", size=0.5, alpha=0.5) +
  geom_vline(xintercept=0, 
             color = "black", size=0.5, alpha=0.5)


p1 + theme_classic()   + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  ylim(-20, 70)  +  xlim(-70, 70)+
  coord_cartesian(xlim =c(-70, 70), ylim = c(-10, 70))+
  scale_y_continuous(breaks = seq(-20, 70, by = 10))+ 
  scale_x_continuous(breaks = seq(-70, 70, by = 20))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))   


#Background Colour
#---------------------

p1 <- ggplot( AllData, aes(LcontrastA, LcontrastB, col = background, fill =  background)) +
  stat_ellipse(geom="polygon", col = "black", alpha = 0.25,) + 
  geom_point(shape =22,  alpha = 0.01, size=3)+
  geom_hline(yintercept=0, 
             color = "black", size=0.5, alpha=0.5) +
  geom_vline(xintercept=0, 
             color = "black", size=0.5, alpha=0.5)

p1 + theme_classic()   + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  coord_cartesian(xlim =c(-70, 70), ylim = c(-10, 70))+
  scale_y_continuous(breaks = seq(-20, 70, by = 10))+ 
  scale_x_continuous(breaks = seq(-70, 70, by = 20))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))   


#Target Luminance
#---------------------
p1 <- ggplot( Gen15PC, aes(Target_Luminance, col = background, fill =  background)) +
  geom_histogram(alpha=0.35,binwidth=5)


p1 + theme_classic()   + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  ylim(0, 100)  +  xlim(0, 100)+
  coord_cartesian(xlim =c(0, 100), ylim = c(0, 100))+
  scale_y_continuous(breaks = seq(-0, 100, by = 10))+ 
  scale_x_continuous(breaks = seq(0, 100, by = 10))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))   


#Background Luminance
#---------------------
p1 <- ggplot( Gen15PC, aes(BgLocal_Luminance, col = background, fill =  background)) +
  geom_histogram(alpha=(0.2),binwidth=5)


p1 + theme_classic()   + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+ 
  ylim(0, 100)  +  xlim(0, 100)+
  coord_cartesian(xlim =c(0, 100), ylim = c(0, 100))+
  scale_y_continuous(breaks = seq(-0, 100, by = 10))+ 
  scale_x_continuous(breaks = seq(0, 100, by = 10))+ 
  scale_colour_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  scale_fill_manual(breaks = c("scrub", "leaf", "veg"),values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5))   









# -----------------------------------------------------------------------------------------------------------------
# (4) SURVIVAL TIME
# -----------------------------------------------------------------------------------------------------------------

AllData$Distance <- sqrt((AllData$X_Coordinate - 739)^2) + sqrt((AllData$Y_Coordinate - 565)^2) 
AllData$Time_Difference = AllData$Capture_Time - AllData$Response_Time
AllData$Time_Used1 = sqrt((AllData$Survival_Time - AllData$Response_Time)^2)
AllData$Time_Used2 = sqrt((AllData$Survival_Time - AllData$Capture_Time)^2)

AllData$Difficulty = ifelse(AllData$Capture_Time>1500,"Best","Worst")


p1 <- ggplot(AllData, aes(Distance, Time_Difference)) +
  geom_smooth(aes(y=Response_Time),linetype="solid",  alpha = 0.1, size=0)+
  stat_smooth(aes(y=Response_Time),linetype="solid",   se=FALSE,  alpha = 0.5, size=0.75)+
  geom_smooth(aes(y=Capture_Time),linetype="dashed",  alpha = 0.1, size=0)+
  stat_smooth(aes(y=Capture_Time),linetype="dashed",   se=FALSE,  alpha = 0.5, size=0.75)+
  geom_smooth(aes(y=Survival_Time), linetype="dotdash",  alpha = 0.1, size=0)+
  stat_smooth(aes(y=Survival_Time), linetype="dotdash", se=FALSE,  alpha = 0.5, size=0.75)


p1 + theme_classic()  + scale_linetype_manual( values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  xlab("Distance (px)") + ylab("Time (milliseconds)")+
  scale_y_continuous(breaks = seq(0, 5000, by = 1000))+ 
  coord_cartesian( ylim = c(1000, 5000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  




p1 <- ggplot(AllData, aes(Distance, Time_Difference)) +
  geom_smooth(  alpha = 0.1, size=0)+
  stat_smooth( se=FALSE,  alpha = 0.5, size=0.75)

p1 + theme_classic()  + scale_linetype_manual(breaks = c("1", "2", "3"),  values=c("dotdash", "longdash","solid"))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size = 14), line = element_line(size=1.4))  + 
  xlab("Distance (px)") + ylab("Time Difference (milliseconds)")+
  scale_y_continuous(breaks = seq(0, 10, by = 0.5))+ 
  scale_y_log10()+
  coord_cartesian( ylim = c(300, 1000))+
  theme(
    axis.title.x = element_text(family = "sans", size = 16, margin=margin(5,0,0,0)), 
    axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,5,0,0)),
    axis.text = element_text(family = "sans", size = 16),
    axis.line = element_line(colour = "black", size=0),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2))  



