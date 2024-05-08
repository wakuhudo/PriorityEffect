library(ggplot2)
library(dplyr)
library(plyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(reshape2)
library(sjPlot)
library(see)
library(brms)
library(DHARMa)
library(cmdstanr)

TMBs1 <- read.csv(file = "TMBs1.csv")
TMBs2NA <- read.csv(file = "TMBs2.csv")
TMBs1NA <- read.csv("TMBs1NA.csv")
TMBs1sex <- read.csv("TMBs1sex.csv")
TMBs2sex <- read.csv("TMBs2sex.csv")


###Start brms
print("Stage 1 Host abundance with ar1, by parasitoid * host, ")

BrmS1HAbun <- brm((Habundance)~  firstP * HostID * Gen + (1|cage) + ar(p = 1) ,
                  data = TMBs1, family = poisson())



sumBrmS1HAbun_tbl <-
  posterior_summary(BrmS1HAbun) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

BrmS1HAbun

sumBrmS1HAbun_tbl

plot(BrmS1HAbun)


SumS1HAbun <- summary(GlmS1HAbun0)
SumS1HAbun

AnovaS1HAbun <- Anova(Glm\S1HAbun0)
AnovaS1HAbun

print("Stage 1 Parasitoid abundance with ar1, by parasitoid * host")

BrmS1PAbunAR1 <- (brm((Pabundance)~  (firstP * Gen*Habundance)  + ar(p=1) + (1|HostID) ,
                      data = TMBs1, family = poisson))
sumBrmS1PAbun_tbl <-
  posterior_summary(BrmS1PAbunAR1) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1PAbun_tbl

plot(BrmS1PAbun)

BrmS1PAbunAR2 <- (brm((Pabundance)~  (firstP * Gen*Habundance)  + ar(p=2) + (1|HostID) ,
                      data = TMBs1, family = poisson))

sumBrmS1PAbunAR2_tbl <-
  posterior_summary(BrmS1PAbunAR2) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1PAbunAR2_tbl
plot(BrmS1PAbunAR2)






AnovaS1PAbun <- Anova(GlmS1PAbun)
AnovaS1PAbun

print("stage 1 sex abundance with ar1 by parasitoid *host")

GlmS1SexAbun <- (glmmTMB((cbind(Male.size,Female.size))~  (firstP * Gen)  + ar1(GenID -1|cage) + (1|HostID) ,
                         data = TMBs1sex, family = binomial))
testZeroInflation(GlmS1SexAbun)
##no sigificance 

SumS1SexAbun <- summary(GlmS1SexAbun)
SumS1SexAbun

AnovaS1SexAbun <- Anova(GlmS1SexAbun)
AnovaS1SexAbun

print("plot abundance interaction with sex")
##explanation?

FigS1sexabun <- plot_model(glmmTMB((cbind(Male.size,Female.size))~  (Gen*firstP)  + ar1(GenID -1|cage) + (1|HostID) ,
                                   data = TMBs1sex, family = binomial), type = "int")
FigS1sexabun

print("Stage 1 parasitism rate with Gen, by P*H")
GlmS1Prate <- (glmmTMB((cbind(Pabundance,nonpara)) ~ firstP*Gen+ (1|HostID)  + ar1(GenID-1|Cage.) , data = TMBs1sex, family = binomial))

testZeroInflation(GlmS1Prate)
##no significance 

SumS1Prate <- summary(GlmS1Prate)
SumS1Prate

AnovaS1Prate <- Anova(GlmS1Prate)
AnovaS1Prate

print("Stage 1 tot resistance rate by sex with gen by P*H")
GlmS1Hres <- (glmmTMB((cbind(resist,nonresist)) ~ (Gen + firstP + sex)^2 + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs1sex, family = binomial))

testZeroInflation(GlmS1Hres)
#no significance 

SumS1Hres <- summary(GlmS1Hres)
SumS1Hres

AnovaS1Hres <- Anova(GlmS1Hres)
AnovaS1Hres

print("Stage 1 tot resistance rate by gen * sex, figure")
GlmS1Hres0 <- (glmmTMB((cbind(resist,nonresist)) ~ (Gen * sex *firstP) + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs1sex, family = binomial))

FigS1sexres <- plot_model(GlmS1Hres0, type = "int")
FigS1sexres


print("Stage 2 Host abundance with ar1, by 2parasitoids*host, block with hostid")

GlmS2HAbun <- (glmmTMB((Habundance) ~  (firstP + secondP + control+ Gen)^3+ ar1(GenID -1|cage)  +(1|HostID),
                       data = TMBs2, family = poisson))
testZeroInflation(GlmS2HAbun)
##zero inflated

SumS2HAbun <- summary(GlmS2HAbun)
SumS2HAbun

AnovaS2HAbun <- Anova(GlmS2HAbun)
AnovaS2HAbun


print("Stage 2 Parasitoid abundance with ar1, by 2parasitoids")

GlmS2PAbun <- (glmmTMB((Pabundance) ~  (firstP +secondP+control +Gen)^3 +ar1(GenID-1|cage)  +(1|HostID),
                       data = TMBs2, family = poisson))

testZeroInflation(GlmS2PAbun)
##zero inflated

Sums2PAbun <- summary(GlmS2PAbun)
Sums2PAbun

Anovas2PAbun <- Anova(GlmS2PAbun)
Anovas2PAbun


print("Stage 2 Host Sex rate with ar1, by 2 Ps")

GlmS2SexAbun <- glmmTMB(cbind(Male.size, Female.size) ~ (firstP + secondP + control + Gen)^3 + ar1(GenID -1|cage) + (1|HostID), 
                        data = TMBs2sex, family = binomial)
testZeroInflation(GlmS2SexAbun)
## no significance 

SumS2SexAbun <- summary(GlmS2SexAbun)
SumS2SexAbun

AnovaS2SexAbun <- Anova(GlmS2SexAbun)
AnovaS2SexAbun

print("Stage 2 plot abundance interaction with sex")

FigS2sexabun <- plot_model(glmmTMB((cbind(Male.size,Female.size))~  (Gen*firstP*secondP*control)  + ar1(GenID -1|cage) + (1|HostID) ,
                                   data = TMBs2sex, family = binomial), type = "int")
FigS2sexabun



print("Stage 2 parasitism rate with Gen, by P*H")

GlmS2Prate <- (glmmTMB((cbind(Pabundance,nonpara)) ~ (firstP + secondP + control + Gen)^3 + ar1(GenID-1|Cage.) + (1|HostID), data = TMBs2, family = binomial))

testZeroInflation(GlmS2Prate)
##Zero deficit?  

SumS2Prate <- summary(GlmS2Prate)
SumS2Prate

AnovaS2Prate <- Anova(GlmS2Prate)
AnovaS2Prate

print("Stage 2 tot resistance rate by sex with gen by P*H*control")


GlmS2Hres <- (glmmTMB((cbind(resist,nonresist)) ~ (firstP+secondP+control +sex +Gen)^3 + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs2sex, family = binomial))
testZeroInflation(GlmS2Hres)
##zero inflated 

SumS2Hres <- summary(GlmS2Hres)
SumS2Hres

AnovaS2Hres <- Anova(GlmS2Hres)
AnovaS2Hres

##now compare with sex abun over time 
AnovaS2SexAbun


print("Stage 2 resistance rate by gen * sex, figure")
GlmS2Hres0 <- (glmmTMB((cbind(resist,nonresist)) ~ (Gen * sex *firstP*secondP) + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs2sex, family = binomial))

FigS2sexres <- plot_model(GlmS2Hres0, type = "int")
FigS2sexres

## For Ar2, choose line 29 r[3,1] 
trace(glmmTMB:::formatVC, edit = T)

### Stage 1 Ar2 Pabundance w/ host abun and p

GlmS1PAbunAr2 <- (glmmTMB((Pabundance)~  (firstP * Gen *Habundance)  + ar1(GenID -1|cage) + (1|HostID),
                          data = TMBs1sex, family = poisson))
#testZeroInflation(GlmS1PAbunAr2)

SummS1PAbunAr2 <- summary(GlmS1PAbunAr2)

SummS1PAbunAr2

SumS1PAbun

AnovaS1PAbunAr2 <- Anova(GlmS1PAbunAr2)
AnovaS1PAbunAr2


### Stage 2 Ar2 Pabundance w/ host abun and p

GlmS2PAbunAr2 <- (glmmTMB((Pabundance)~  (firstP +secondP+control +Gen+Habundance)^3  + ar1(GenID -1|cage) + (1|HostID),
                          data = TMBs2, family = poisson))
testZeroInflation(GlmS2PAbunAr2)
##not zero inflated 

sumS2PAbunAr2 <- summary(GlmS2PAbunAr2)

sumS2PAbunAr2

Sums2PAbun

AnovaS2PAbunAr2 <- Anova(GlmS2PAbunAr2)

AnovaS2PAbunAr2

Anovas2PAbun

library(forecast)
fit <- auto.arima(TMBs2$Pabundance, xreg=cbind(TMBs2$Gen) )
summary(fit)
coef(fit)

arima(TMBs2$Pabundance, order = c(2,0,0), xreg = TMBs2$Gen)
library(nlme)
lme(log(Pabundance) ~ Habundance + firstP,  data = TMBs2,random = ~ 1|Cage.,correlation = corARMA(p =2, q = 1) )

###test for PCA GEn 8,9,10

##a PC exploration 
library(factoextra)
Popsizeinter[is.na(Popsizeinter)] <-0
pcadata  <- Popsizeinter[,c(5,6,8,9,28:30)]
pcadata <- na.omit(pcadata)
pca <- prcomp((pcadata),  scale = T)
PC <- as.data.frame(pca$x)

Popsizeinter <- Popsizeinter %>% drop_na(Male.BM)


Popsizeinter$pc1 <- PC$PC1
Popsizeinter$pc2 <- PC$PC2


PCAbyGenID <- ggplot(Popsizeinter, 
                     aes(x = pc1, 
                         y = pc2, 
                         color = as.factor(GenID))) +
  geom_point() +
  stat_ellipse() + ggtitle("PCA w/ Ellipse by GenID") +
  theme_bw()


BiplotbyPredictor <- fviz_pca_biplot(pca, 
                                     repel = TRUE,
                                     col.var = "deepskyblue",
                                     title = "Biplot w/ Predictors", geom="point")

##present from here 

PCAbyGenID

BiplotbyPredictor

manovainter <- manova(cbind(Male.size, Male.BM,Female.size,Fe.bm,Egg, Larvae,Pupae) ~ as.factor(GenID), data = Popsizeinter)
summa <- summary.manova(manovainter)
summa
summa$Eigenvalues



anovasexrateinter <- aov(Malerate ~ GenID, data = Popsizeinter)
summary(anovasexrateinter)

anovafecuncityinter <- aov(TotFecundity ~ GenID, data = Popsizeinter)
anovaegginter <- aov(Egg ~GenID, data = Popsizeinter)

summary(anovaegginter)
summary(anovafecuncityinter)

TMBs2$Malerate <- TMBs2$Male.size/TMBs2$Habundance
ggplot(TMBs2, aes(as.numeric(GenID),Malerate)) + geom_point()+ theme_bw()+
  geom_smooth(method = "loess") +ggtitle("2nd gen male rate w/o is.na <- 0")

TMBs2[is.na(TMBs2)] <-0

ggplot(TMBs2, aes(as.numeric(GenID),Malerate)) + geom_point()+ theme_bw()+
  geom_smooth(method = "loess") +ggtitle("2nd gen male rate w/ is.na <- 0")


##partial Correlation 
Pabun <- PopsizeP$Pabundance[1:1280]
Pabun <- as.data.frame(Pabun)
Pabun$AR1P<- unlist(c(rep(0,160),PopsizeP$Pabundance[1:1120]))
Pabun$AR2P <- unlist(c(rep(0,320),PopsizeP$Pabundance[1:960]))
Pabun$AR3P <- unlist(c(rep(0,480),PopsizeP$Pabundance[1:800]))

Pabun$cage <- PopsizeP$Cage.[1:1280]

#glmcorP <- glmmTMB(Pabun ~ AR1P + AR2P+ (1|cage), data = pcorP)
#R2_BM <- partR2(glmcorP, data = pcorP, R2_type = "marginal", nboot = 10)
library(correlation)
library(lme4)

correlation(Pabun, multilevel = T,partial = T)


Pabun2 <- PopsizeP$Pabundance[1281:2560]
Pabun2 <- as.data.frame(Pabun2)
Pabun2$AR1P<- unlist(c(rep(0,160),PopsizeP$Pabundance[1281:2400]))
Pabun2$AR2P <- unlist(c(rep(0,320),PopsizeP$Pabundance[1281:2240]))
Pabun2$AR3P <- unlist(c(rep(0,480),PopsizeP$Pabundance[1281:2]))
Pabun2$cage <- PopsizeP$Cage.[1281:2560]

correlation(Pabun2, multilevel = T,partial = T)


#### GLMMTMB old cript 
(glmer(log(as.numeric(Habundance)) ~ Gen * firstP + (1|cage), data = TMBs1, family = poisson ))

##extremely zero inflated 


AnovaS1PAbun <- Anova(GlmS1PAbun)
AnovaS1PAbun

print("stage 1 sex abundance with ar1 by parasitoid *host")

GlmS1SexAbun <- (glmmTMB((cbind(Male.size,Female.size))~  (firstP * Gen)  + ar1(GenID -1|cage) + (1|HostID) ,
                         data = TMBs1sex, family = binomial))
testZeroInflation(GlmS1SexAbun)
##no sigificance 

SumS1SexAbun <- summary(GlmS1SexAbun)
SumS1SexAbun

AnovaS1SexAbun <- Anova(GlmS1SexAbun)
AnovaS1SexAbun

print("plot abundance interaction with sex")
##explanation?

FigS1sexabun <- plot_model(glmmTMB((cbind(Male.size,Female.size))~  (Gen*firstP)  + ar1(GenID -1|cage) + (1|HostID) ,
                                   data = TMBs1sex, family = binomial), type = "int")
FigS1sexabun

print("Stage 1 parasitism rate with Gen, by P*H")
GlmS1Prate <- (glmmTMB((cbind(Pabundance,nonpara)) ~ firstP*Gen+ (1|HostID)  + ar1(GenID-1|Cage.) , data = TMBs1sex, family = binomial))

testZeroInflation(GlmS1Prate)
##no significance 

SumS1Prate <- summary(GlmS1Prate)
SumS1Prate

AnovaS1Prate <- Anova(GlmS1Prate)
AnovaS1Prate

print("Stage 1 tot resistance rate by sex with gen by P*H")
GlmS1Hres <- (glmmTMB((cbind(resist,nonresist)) ~ (Gen + firstP + sex)^2 + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs1sex, family = binomial))

testZeroInflation(GlmS1Hres)
#no significance 

SumS1Hres <- summary(GlmS1Hres)
SumS1Hres

AnovaS1Hres <- Anova(GlmS1Hres)
AnovaS1Hres

print("Stage 1 tot resistance rate by gen * sex, figure")
GlmS1Hres0 <- (glmmTMB((cbind(resist,nonresist)) ~ (Gen * sex *firstP) + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs1sex, family = binomial))

FigS1sexres <- plot_model(GlmS1Hres0, type = "int")
FigS1sexres


print("Stage 2 Host abundance with ar1, by 2parasitoids*host, block with hostid")

GlmS2HAbun <- (glmmTMB((Habundance) ~  (firstP + secondP + control+ Gen)^3+ ar1(GenID -1|cage)  +(1|HostID),
                       data = TMBs2, family = poisson))
testZeroInflation(GlmS2HAbun)
##zero inflated

SumS2HAbun <- summary(GlmS2HAbun)
SumS2HAbun

AnovaS2HAbun <- Anova(GlmS2HAbun)
AnovaS2HAbun


print("Stage 2 Parasitoid abundance with ar1, by 2parasitoids")

GlmS2PAbun <- (glmmTMB((Pabundance) ~  (firstP +secondP+control +Gen)^3 +ar1(GenID-1|cage)  +(1|HostID),
                       data = TMBs2, family = poisson))

testZeroInflation(GlmS2PAbun)
##zero inflated

Sums2PAbun <- summary(GlmS2PAbun)
Sums2PAbun

Anovas2PAbun <- Anova(GlmS2PAbun)
Anovas2PAbun


print("Stage 2 Host Sex rate with ar1, by 2 Ps")

GlmS2SexAbun <- glmmTMB(cbind(Male.size, Female.size) ~ (firstP + secondP + control + Gen)^3 + ar1(GenID -1|cage) + (1|HostID), 
                        data = TMBs2sex, family = binomial)
testZeroInflation(GlmS2SexAbun)
## no significance 

SumS2SexAbun <- summary(GlmS2SexAbun)
SumS2SexAbun

AnovaS2SexAbun <- Anova(GlmS2SexAbun)
AnovaS2SexAbun

print("Stage 2 plot abundance interaction with sex")

FigS2sexabun <- plot_model(glmmTMB((cbind(Male.size,Female.size))~  (Gen*firstP*secondP*control)  + ar1(GenID -1|cage) + (1|HostID) ,
                                   data = TMBs2sex, family = binomial), type = "int")
FigS2sexabun



print("Stage 2 parasitism rate with Gen, by P*H")

GlmS2Prate <- (glmmTMB((cbind(Pabundance,nonpara)) ~ (firstP + secondP + control + Gen)^3 + ar1(GenID-1|Cage.) + (1|HostID), data = TMBs2, family = binomial))

testZeroInflation(GlmS2Prate)
##Zero deficit?  

SumS2Prate <- summary(GlmS2Prate)
SumS2Prate

AnovaS2Prate <- Anova(GlmS2Prate)
AnovaS2Prate

print("Stage 2 tot resistance rate by sex with gen by P*H*control")


GlmS2Hres <- (glmmTMB((cbind(resist,nonresist)) ~ (firstP+secondP+control +sex +Gen)^3 + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs2sex, family = binomial))
testZeroInflation(GlmS2Hres)
##zero inflated 

SumS2Hres <- summary(GlmS2Hres)
SumS2Hres

AnovaS2Hres <- Anova(GlmS2Hres)
AnovaS2Hres

##now compare with sex abun over time 
AnovaS2SexAbun


print("Stage 2 resistance rate by gen * sex, figure")
GlmS2Hres0 <- (glmmTMB((cbind(resist,nonresist)) ~ (Gen * sex *firstP*secondP) + ar1(GenID-1|Cage.)+(1|HostID) , data = TMBs2sex, family = binomial))

FigS2sexres <- plot_model(GlmS2Hres0, type = "int")
FigS2sexres

## For Ar2, choose line 29 r[3,1] 
trace(glmmTMB:::formatVC, edit = T)

### Stage 1 Ar2 Pabundance w/ host abun and p

GlmS1PAbunAr2 <- (glmmTMB((Pabundance)~  (firstP * Gen *Habundance)  + ar1(GenID -1|cage) + (1|HostID),
                          data = TMBs1sex, family = poisson))
#testZeroInflation(GlmS1PAbunAr2)

SummS1PAbunAr2 <- summary(GlmS1PAbunAr2)

SummS1PAbunAr2

SumS1PAbun

AnovaS1PAbunAr2 <- Anova(GlmS1PAbunAr2)
AnovaS1PAbunAr2


### Stage 2 Ar2 Pabundance w/ host abun and p

GlmS2PAbunAr2 <- (glmmTMB((Pabundance)~  (firstP +secondP+control +Gen+Habundance)^3  + ar1(GenID -1|cage) + (1|HostID),
                          data = TMBs2, family = poisson))
testZeroInflation(GlmS2PAbunAr2)
##not zero inflated 

sumS2PAbunAr2 <- summary(GlmS2PAbunAr2)

sumS2PAbunAr2

Sums2PAbun

AnovaS2PAbunAr2 <- Anova(GlmS2PAbunAr2)

AnovaS2PAbunAr2

Anovas2PAbun

library(forecast)
fit <- auto.arima(TMBs2$Pabundance, xreg=cbind(TMBs2$Gen) )
summary(fit)
coef(fit)

arima(TMBs2$Pabundance, order = c(2,0,0), xreg = TMBs2$Gen)
library(nlme)
lme(log(Pabundance) ~ Habundance + firstP,  data = TMBs2,random = ~ 1|Cage.,correlation = corARMA(p =2, q = 1) )

###test for PCA GEn 8,9,10

##a PC exploration 
library(factoextra)
Popsizeinter[is.na(Popsizeinter)] <-0
pcadata  <- Popsizeinter[,c(5,6,8,9,28:30)]
pcadata <- na.omit(pcadata)
pca <- prcomp((pcadata),  scale = T)
PC <- as.data.frame(pca$x)

Popsizeinter <- Popsizeinter %>% drop_na(Male.BM)


Popsizeinter$pc1 <- PC$PC1
Popsizeinter$pc2 <- PC$PC2


PCAbyGenID <- ggplot(Popsizeinter, 
                     aes(x = pc1, 
                         y = pc2, 
                         color = as.factor(GenID))) +
  geom_point() +
  stat_ellipse() + ggtitle("PCA w/ Ellipse by GenID") +
  theme_bw()
