#load packages
library(ggplot2)
library(moments)
library(MASS)
library(dplyr)
library(emmeans)
library(viridis)


dat <- read.table("Finished Data.txt", header=TRUE)

dat$Group <- as.factor(dat$Group)

summary(dat)
par(mfrow=c(1,1))
hist(dat$MaxThreat)
skewness(dat$MaxThreat)
#Data is negligibly positively skewed
hist(dat$CurrentThreat)
skewness(dat$CurrentThreat)
#Data is highly positively skewed

#BoxCox transformation for Current Threat
cbc<-boxcox(CurrentThreat+1~1, data=dat)
lambda<-cbc$x[which.max(cbc$y)]
dat$CurrentThreatBC <- ((dat$CurrentThreat+1)^lambda-1)/lambda
hist(dat$CurrentThreatBC)
skewness(dat$CurrentThreatBC)
#skew is largely reduced



#Section 1: Current Threat Analysis

summary(dat)
modelx <- lm(CurrentThreatBC ~ Heterogeneity + Population + Group, data = dat)
summary(modelx)
#Overall model including all variables

#Check assumptions of model
par(mfrow=c(2,2))
plot(modelx)
#create histogram of residuals
par(mfrow=c(1,1))
ggplot(data = dat, aes(x = modelx$residuals)) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


#Check whether island grouping improves model
modely <- lm(CurrentThreatBC ~ Heterogeneity + Population, data = dat)
anova(modelx, modely)
#Grouping improves model non-significantly

#Check whether Population improves model
modelz <- lm(CurrentThreatBC ~ Heterogeneity + Group, data = dat)
anova(modelx, modelz)
#Population improves model significantly

#Check whether Heterogeneity improves model
modelw <- lm(CurrentThreatBC ~ Population + Group, data = dat)
anova(modelx, modelw)
#Heterogeneity improves model but not significantly

#Minimum Adequate Model:
modelv <- lm(CurrentThreatBC ~ Population, data = dat)
summary(modelv)

#Chosen final model:
modelx %>% summary()

#Creating 95% CI values for effect of heterogeneity and human population:
#Define lambda value
bc_tran <- make.tran("boxcox", alpha = lambda)

#Fit linear model INSIDE environment
modelx_forCI <- with(bc_tran, 
              lm(linkfun(CurrentThreat) ~ Heterogeneity + Population + Group, 
                 data = dat))

#Get the back-transformed estimate and 95% CI for continuous predictors
slope_x1 <- emtrends(modelx_forCI, ~ Heterogeneity, var = "Heterogeneity", type = "response")
cat("--- Estimates and 95% CI for Heterogeneity ---\n")
print(slope_x1)

slope_x2 <- emtrends(modelx_forCI, ~ Population, var = "Population", type = "response")
cat("\n--- Estimates and 95% CI for Population ---\n")
print(slope_x2)


#Section 1.5: Analysis of difference in current threat between island groups

groupmodel <- aov(CurrentThreatBC ~ Group, data=dat)
summary(groupmodel)
#Group means are not significantly different

TukeyHSD(groupmodel, which = "Group")
#No significant differences between islands



#Section 2: Maximum Threat Analysis


#Fitting a Gamma GLM with log link to account for slight skew:
modelxx <- glm(MaxThreat ~ Heterogeneity + Group + Population, data = dat, family=(Gamma(link="log")))
summary(modelxx)

#Check assumptions of model
par(mfrow=c(2,2))
plot(modelxx)
#create histogram of residuals
ggplot(data = dat, aes(x = modelxx$residuals)) +
  geom_histogram(fill = 'white', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

#Extract the raw coefficients (on the log scale)
raw_coefs <- coef(modelxx)

#Extract the confidence intervals (on the log scale)
raw_cis <- confint(modelxx)

#Back-transformation
back_transformed_coefs <- exp(raw_coefs)
back_transformed_cis   <- exp(raw_cis)

#Combine into a results table
results_table <- cbind(
  Exponentiated_Estimate = back_transformed_coefs, 
  Lower_95_CI = back_transformed_cis[,1], 
  Upper_95_CI = back_transformed_cis[,2])
results_table


#Section 2.5: Analysis of difference in current threat between island groups

groupmodel2 <- aov(MaxThreat ~ Group, data=dat)
summary(groupmodel2)
#Group means are not significantly different

TukeyHSD(groupmodel2, which = "Group")
#No significant differences between islands



#Section 3: Plotting
#3.1: Current Threat ~ Population

#Dataset for prediction line
PopDat1<- data.frame(
  Population = seq(min(dat$Population),
                   max(dat$Population),
                   length.out = 100),
  Heterogeneity = mean(dat$Heterogeneity),
  Group = factor(rep("Highland", 100),
                 levels = levels(dat$Group))
)
#Creating prediction
pred.pop1 <- predict(modelx,
                     newdata = PopDat1,
                     interval = "confidence")
#Backtransforming upper and lower confidence bands
PopDat1$fit <- ((pred.pop1[,"fit"]*lambda)+1)^(1/lambda)-1
PopDat1$lwr <- ((pred.pop1[,"lwr"]*lambda)+1)^(1/lambda)-1
PopDat1$upr <- ((pred.pop1[,"upr"]*lambda)+1)^(1/lambda)-1
#Plotting line plot
ggplot(PopDat1,
       aes(x = Population,
           y = fit)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(x = "Island Population",
       y = "Predicted Current Threat Posed by Non-Native Species") +
  theme_classic()


#3.2: Current Threat ~ Heterogeneity

#Dataset for prediction line
HetDat1<- data.frame(
  Heterogeneity = seq(min(dat$Heterogeneity),
                   max(dat$Heterogeneity),
                   length.out = 100),
  Population = mean(dat$Population),
  Group = factor(rep("Highland", 100),
                 levels = levels(dat$Group))
)
#Creating prediction
pred.het1 <- predict(modelx,
                     newdata = HetDat1,
                     interval = "confidence")
#Backtransforming upper and lower confidence bands
HetDat1$fit <- ((pred.het1[,"fit"]*lambda)+1)^(1/lambda)-1
HetDat1$lwr <- ((pred.het1[,"lwr"]*lambda)+1)^(1/lambda)-1
HetDat1$upr <- ((pred.het1[,"upr"]*lambda)+1)^(1/lambda)-1
#Plotting line plot
ggplot(HetDat1,
       aes(x = Heterogeneity,
           y = fit)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(x = "Island Habitat Heterogeneity",
       y = "Predicted Current Threat Posed by Non-Native Species") +
  theme_classic()



#3.3: Current Threat ~ Island Group

#Dataset for predicted values whilst accounting for 
#Heterogeneity and Population
GroupDat1 <- data.frame(
  Heterogeneity = mean(dat$Heterogeneity),
  Population = mean(dat$Population),
  Group = factor(levels(dat$Group),
                 levels = levels(dat$Group))
)

#Creating prediction
pred.group1 <- predict(modelx,
                newdata = GroupDat1,
                interval = "confidence")

#Backtransforming upper and lower confidence intervals
GroupDat1$fit <- ((pred.group1[,"fit"]*lambda)+1)^(1/lambda)-1
GroupDat1$lwr <- ((pred.group1[,"lwr"]*lambda)+1)^(1/lambda)-1
GroupDat1$upr <- ((pred.group1[,"upr"]*lambda)+1)^(1/lambda)-1

#Plotting boxplot
ggplot(dat, aes(x = Group, y = CurrentThreat)) +
  geom_boxplot(alpha = 0.5) +
  geom_boxplot(staplewidth = 0.5) +
 # geom_point(data = GroupDat1,
          #   aes(x = Group, y = fit),
         #    colour = "red",
          #   size = 3,
         #    inherit.aes = FALSE) +
#  geom_errorbar(data = GroupDat1,
              #  aes(x = Group,
                 #   ymin = lwr,
                 #   ymax = upr),
               # colour = "red",
               # width = 0.2,
                #inherit.aes = FALSE) +
  labs(x = "Island Group", y = "Current Threat Posed by Non-Native Species") +
  theme_classic()


#3.4: Max Threat ~ Population

modelxxpred <- lm(MaxThreat ~ Heterogeneity + Group + Population, data = dat)
summary(modelxxpred)

#Dataset for prediction line
PopDat2<- data.frame(
  Population = seq(min(dat$Population),
                   max(dat$Population),
                   length.out = 100),
  Heterogeneity = mean(dat$Heterogeneity),
  Group = factor(rep("Highland", 100),
                 levels = levels(dat$Group))
)
#Creating prediction
pred.pop2 <- predict(modelxxpred,
                     newdata = PopDat2,
                     interval = "confidence")
#Backtransforming upper and lower confidence bands
PopDat2$fit <- (pred.pop2[,"fit"]^(1/4))
PopDat2$lwr <- (pred.pop2[,"lwr"]^(1/4))
PopDat2$upr <- (pred.pop2[,"upr"]^(1/4))
#Plotting line plot
ggplot(PopDat2,
       aes(x = Population,
           y = fit)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(x = "Island Population",
       y = "Predicted Maximum Threat Posed by Non-Native Species") +
  theme_classic()




#3.5: Max Threat ~ Heterogeneity

#Dataset for prediction line
HetDat2<- data.frame(
  Heterogeneity = seq(min(dat$Heterogeneity),
                      max(dat$Heterogeneity),
                      length.out = 100),
  Population = mean(dat$Population),
  Group = factor(rep("Highland", 100),
                 levels = levels(dat$Group))
)
#Creating prediction
pred.het2 <- predict(modelxxpred,
                     newdata = HetDat2,
                     interval = "confidence")
#Backtransforming upper and lower confidence bands
HetDat2$fit <- (pred.het2[,"fit"]^(1/4))
HetDat2$lwr <- (pred.het2[,"lwr"]^(1/4))
HetDat2$upr <- (pred.het2[,"upr"]^(1/4))
#Plotting line plot
ggplot(HetDat2,
       aes(x = Heterogeneity,
           y = fit)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(x = "Island Habitat Heterogeneity",
       y = "Predicted Maximum Threat Posed by Non-Native Species") +
  theme_classic()



#3.6: Max Threat ~ Island Group

#Dataset for predicted values whilst accounting for 
#Heterogeneity and Population
GroupDat2 <- data.frame(
  Heterogeneity = mean(dat$Heterogeneity),
  Population = mean(dat$Population),
  Group = factor(levels(dat$Group),
                 levels = levels(dat$Group))
)

#Creating prediction
pred.group2 <- predict(modelxxpred,
                       newdata = GroupDat2,
                       interval = "confidence")

#Backtransforming upper and lower confidence intervals
GroupDat2$fit <- pred.group2[,"fit"]^(1/4)
GroupDat2$lwr <- pred.group2[,"lwr"]^(1/4)
GroupDat2$upr <- pred.group2[,"upr"]^(1/4)

#Plotting boxplot
ggplot(dat, aes(x = Group, y = MaxThreat)) +
  geom_boxplot(alpha = 0.5) +
  geom_boxplot(staplewidth = 0.5) +
  geom_point(data = GroupDat2,
             aes(x = Group, y = fit),
             colour = "red",
             size = 3,
             inherit.aes = FALSE) +
  geom_errorbar(data = GroupDat2,
                aes(x = Group,
                    ymin = lwr,
                    ymax = upr),
                colour = "red",
                width = 0.2,
                inherit.aes = FALSE) +
  coord_cartesian(ylim = c(5, 15)) +
  labs(x = "Island Group", y = "Maximum Threat Posed by Non-Native Species") +
  theme_classic()




#Section 4: Graphs of Literature Review
#4.1: Impact Score and Confidence
ICdat <- read.table("ImpactConfidence.txt", header=TRUE)
ICdat$Score<-as.factor(ICdat$Score)
ICdat$Confidence<-as.factor(ICdat$Confidence)
ScoreOrder <- c("MC", "Minor", "Moderate", "Major", "Massive")
ICdat$Score <- factor(ICdat$Score, levels = ScoreOrder)
ICtotals <- ICdat %>% count(Score)
ICdat$Confidence <- factor(ICdat$Confidence, levels = c("High", "Medium", "Low"))


ggplot(ICdat, aes(x = Score, fill=Confidence)) +
  geom_bar(color = "black") +
  labs(
    title = "",
    x = "EICAT Impact Score",
    y = "Number of Publications",
    fill = "Confidence"
  ) + 
  scale_fill_viridis_d(option = "viridis") +
  scale_x_discrete(labels = c(
    "MC" = "Minimal Concern",
    "Minor" = "Minor",
    "Moderate" = "Moderate",
    "Major" = "Major",
    "Massive" = "Massive"
  )) + 
  geom_text(
    data = ICtotals,
    aes(x = Score, y = n, label = n),
    inherit.aes = FALSE,
    vjust = -0.5
  ) +
  theme_classic()

#4.2: Species, Mechanism and Impact Score
SMdat <- read.table("SpeciesMechanism.txt", header=TRUE)
SMdat$Species<-as.factor(SMdat$Species)
SMdat$Confidence<-as.factor(SMdat$Mechanism)
SMtotals <- SMdat %>% count(Species)


ggplot(SMdat, aes(x = Species, fill=Mechanism)) +
  geom_bar(color = "black") +
  labs(
    title = "",
    x = "Species",
    y = "Number of Publications",
    fill = "Impact Mechanism"
  ) + 
  scale_fill_viridis_d(option = "viridis") +
  scale_x_discrete(labels = c(
    "Bcanadensis" = "Branta canadensis",
    "Cnippon" = "Cervus nippon",
    "Cselloana" = "Cortaderia selloana",
    "Elodea" = "Elodea canadensis/nuttallii",
    "Gunnera" = "Gunntera tinctoria/cryptica",
      "Hranunculoides" = "Hydrocotyle ranunculoides",
      "Lamericanus" = "Lysichiton americanus",
      "Nvison" = "Neovison vison",
      "Rponticum" = "Rhododendron ponticum"
  )) +
  geom_text(
    data = SMtotals,
    aes(x = Species, y = n, label = n),
    inherit.aes = FALSE,
    vjust = -0.5
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#4.3: Species and Area of Study

SAdat <- read.table("SpeciesArea.txt", header=TRUE)
SAdat$Climate<-as.factor(SAdat$Climate)
SAdat$Species<-as.factor(SAdat$Species)
SAtotals <- SAdat %>% count(Climate)

ClimateOrder <- c("Scotisland", "Scot", "UKisland", "UK", "EUislandtemp", "EUtemp", "EUisland", "EU", "NonEUisland", "NonEU", "NA")
SAdat$Climate <- factor(SAdat$Climate, levels = ClimateOrder)

ggplot(SAdat, aes(x = Climate, fill=Species)) +
  geom_bar(color = "black") +
  labs(
    title = "",
    x = "Location of Evidence Collection",
    y = "Number of Publications",
    fill = "Species"
  ) + 
  scale_fill_viridis_d(option = "viridis",labels = c(
    "Bcanadensis" = "Branta canadensis",
    "Cnippon" = "Cervus nippon",
    "Cselloana" = "Cortaderia selloana",
    "Elodea" = "Elodea canadensis/nuttallii",
    "Gunnera" = "Gunntera tinctoria/cryptica",
    "Hranunculoides" = "Hydrocotyle ranunculoides",
    "Lamericanus" = "Lysichiton americanus",
    "Nvison" = "Neovison vison",
    "Rponticum" = "Rhododendron ponticum"
  )) +
  scale_x_discrete(labels = c(
    "Scotisland" = "Scottish Island", 
    "Scot" = "Mainland Scotland", 
    "UKisland" = "Island of Ireland, England or Wales", 
    "UK" = "Mainland Ireland, England or Wales", 
    "EUislandtemp" = "Island in Temperate Europe", 
    "EUtemp" = "Temperate Europe", 
    "EUisland" = "Island in non-temperate Europe", 
    "EU" = "Non-temperate Europe", 
    "NonEUisland" = "Island outside of Europe", 
    "NonEU" = "Mainland outside of Europe", 
    "NA" = "Ex-situ"
  )) +
  geom_text(
    data = SAtotals,
    aes(x = Climate, y = n, label = n),
    inherit.aes = FALSE,
    vjust = -0.5
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dat$Diff <- dat$MaxThreat-dat$CurrentThreat


