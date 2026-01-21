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
library(marginaleffects)
library(bayesplot)
library(reshape2)
library(patchwork)
library(rstan)
library(StanHeaders)
library(tidybayes)
library(projpred)
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
TMBs2sex <- read.csv(file = "/Users/Apple/Desktop/PhD/Code/Data Analysis_Project1/May9/TMBs2sex.csv")
TMBs2NA <- read.csv(file = "/Users/Apple/Desktop/PhD/Code/PriorityEffect/TMBs2.csv")
TMBs1NA <- read.csv("/Users/Apple/Desktop/PhD/Code/PriorityEffect/TMBs1NA.csv")
TMBs1NA <-TMBs1NA[1:1120,]

mean((subset(TMBs1NA, GenID == 2)$Male.size+subset(TMBs1NA, GenID == 2)$Female.size)/5)
sd((subset(TMBs1NA, GenID == 2)$Male.size+subset(TMBs1NA, GenID == 2)$Female.size)/5)
sd((subset(TMBs1NA, GenID == 2)$totavalar)/5)
TMBs2NA$firstP <- factor(TMBs2NA$firstP)
TMBs2NA$secondP <- factor(TMBs2NA$secondP)

###Start brms
print("Stage 1 Host abundance with ar1, by parasitoid * host, ")

BrmS1HAbun <- brm((Habundance)~  firstP * HostID * Gen + (1|cage) + ar(p = 1) ,
                  data = TMBs1, family = poisson())

saveRDS(BrmS1HAbun, file = "Model/BrmS1HAbun.rds")
BrmS1HAbun<- readRDS(file = "Model/BrmS1HAbun.rds")

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

print("Stage 1 Parasitoid abundance with ar1, by parasitoid * host")

BrmS1PAbunAR1 <- (brm((Pabundance)~  (firstP * Gen)  + ar(p=1) + (1|HostID) ,
                      data = TMBs1, family = poisson))

saveRDS(BrmS1PAbunAR1, file = "Model/BrmS1PAbunAR1.rds")
BrmS1PAbunAR1<- readRDS(file = "Model/BrmS1PAbunAR1.rds")


sumBrmS1PAbun_tbl <-
  posterior_summary(BrmS1PAbunAR1) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1PAbun_tbl

plot(BrmS1PAbunAR1)


#print("stage 1 sex abundance with ar1 by parasitoid *host")

#GlmS1SexAbun <- (glmmTMB((cbind(Male.size,Female.size))~  (firstP*Gen)   + (1|HostID)+(1|Cage.) ,
#                         data =na.omit(TMBs1NA), family = binomial))


conditional_effects(BrmS1Prate0, "Gen:HostID:firstP")

print("Stage 1 parasitism rate with Gen, by P*H")
BrmS1Prate0 <- brm(Pabundance | trials(totavalar) ~ firstP * mo(Gen)*HostID + ar(p = 1) + (1|Cage.) ,
                   data = TMBs1NA , cores = 4 , family = binomial)
BrmS1Prate
##back here


draws1 <- as_draws_df(BrmS1Prate0)
draws <-
  draws1 %>% 
  select(starts_with("b")) 
draws[,5:8] <- draws[,5:8] * 7


f <- function(x) {
  r <- quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}



melted_draws <- melt(draws)

# Calculate lower and upper quartiles for each variable

quartiles<- data.frame()

for (i in unique(melted_draws$variable) ) {
  test <- subset(melted_draws, variable == i)
  quartiles0 <- test %>%
    summarize(
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    mutate(
      asterisk = ifelse((lower > 0 | upper < 0), "*", ""),
      variable = i)
  quartiles <- rbind(quartiles,quartiles0)
}


# Check the output of quartiles
print(quartiles)
S1Prate <- ggplot(data = melted_draws, aes(x = value, y = variable)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fill = "black") +stat_summary(fun.data = f, geom="boxplot")+
  scale_y_discrete(limits = rev) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  ggtitle("2nd Stage Para") +
  geom_text(data = quartiles, aes(x = 4, y = variable, label = asterisk), 
            size = 5, color = "red")


ce <- conditional_effects(BrmS1Prate0, effects = "Gen", conditions = data.frame(firstP = rep(c("T", "H"), each = 2),
                                                                               HostID = rep(c("M", "S"), times = 2)))

# Convert the list to a data frame
ce_df <- ce[[1]]

# Create the facet plot
S1CorrPrate <- ggplot(ce_df, aes(x = Gen, y = estimate__)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  facet_grid(firstP ~ HostID) +
  labs(x = "Gen", y = "Estimated Prate", title = "Conditional Effects of Gen on S1 Prate") +
  theme_minimal()

library(ggpubr)




saveRDS(BrmS1Prate0, file = "Model/BrmS1Prate0.rds")
BrmS1PrateAR1<- readRDS(file = "Model/BrmS1PrateAR1.rds")

library(tidybayes)
fixgen <-   as_draws_df(BrmS1Prate0) %>% 
  transmute(Estimate = (bsp_moGen)*7) %>% 
  median_qi(.width = .95) %>% 
  mutate_if(is.double, round, digits = 2)




sumBrmS1PrateAR1_tbl <-
  posterior_summary(BrmS1PrateAR1) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1PrateAR1_tbl
plot(BrmS1PrateAR1)


print("Stage 1 Parasitoid abundance with ar2, by parasitoid * host")
BrmS1PAbunAR2 <- (brm((Pabundance)~  (firstP * Gen )  + ar(p=2) +(1|Cage.)+ (1|HostID) ,
                      data = TMBs1, family = poisson))

saveRDS(BrmS1PAbunAR2, file = "Model/BrmS1PAbunAR2.rds")
BrmS1PAbunAR2<- readRDS(file = "Model/BrmS1PAbunAR2.rds")


sumBrmS1PAbunAR2_tbl <-
  posterior_summary(BrmS1PAbunAR2) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1PAbunAR2_tbl
plot(BrmS1PAbunAR2)

BrmS1PAbunAR2TestHabun
plot(BrmS1PAbunAR2TestHabun)

###left from here: 


print("Stage 1 parasitism rate with Gen, ar2, by P*H")
BrmS1PrateAR2 <- brm(Pabundance | trials(totavalar) ~ firstP * Gen +ar(p = 2) + (1|Cage.) + (1|HostID)  , 
                     data = TMBs1NA, family = binomial)


saveRDS(BrmS1PrateAR2, file = "Model/BrmS1PrateAR2.rds")
BrmS1PrateAR2<- readRDS(file = "Model/BrmS1PrateAR2.rds")


sumBrmS1PrateAR2_tbl <-
  posterior_summary(BrmS1PrateAR2) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1PrateAR2_tbl
plot(BrmS1PrateAR2)



TMBs1NA$totdot <- is.numeric(TMBs1NA$totdot)
TMBs1NA$totdot[is.na(TMBs1NA$totdot)] <- 0

print("Stage 1 tot resistance rate by sex with gen by P*H")
BrmS1Hres0 <- brm((totdot) | trials(Pabundance+totdot)  ~ mo(Gen) * firstP *HostID+ ar(p=1)+(1|Cage.) , data = TMBs1NA,cores = 4, family = binomial)
BrmS1Hres

saveRDS(BrmS1Hres0, file = "Model/BrmS1Hres0.rds")
BrmS1Hres<- readRDS(file = "Model/BrmS1Hres.rds")


sumBrmS1Hres_tbl <-
  posterior_summary(BrmS1Hres) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS1Hres_tbl
plot(BrmS1Hres)




draws1 <- as_draws_df(BrmS1Hres0)
draws <-
  draws1 %>% 
  select(starts_with("b")) 
draws[,5:8] <- draws[,5:8] * 7


f <- function(x) {
  r <- quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}


melted_draws <- melt(draws)

# Calculate lower and upper quartiles for each variable

quartiles<- data.frame()

for (i in unique(melted_draws$variable) ) {
  test <- subset(melted_draws, variable == i)
  quartiles0 <- test %>%
    summarize(
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    mutate(
      asterisk = ifelse((lower > 0 | upper < 0), "*", ""),
      variable = i)
  quartiles <- rbind(quartiles,quartiles0)
}


# Check the output of quartiles
print(quartiles)
S1Hres <- ggplot(data = melted_draws, aes(x = value, y = variable)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fill = "black") +stat_summary(fun.data = f, geom="boxplot")+
  scale_y_discrete(limits = rev) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  ggtitle("1st Stage Hres") +
  geom_text(data = quartiles, aes(x = 4, y = variable, label = asterisk), 
            size = 5, color = "red")




S1Res <- ggplot(data = melt(draws), aes(x=value, y=variable)) + geom_boxplot(width = 0.01, outlier.shape = NA,fill = "black")+     scale_y_discrete(limits=rev)+
  stat_summary(fun.data = f, geom="boxplot")+theme_bw() + geom_vline(xintercept = 0) + ggtitle("1st Stage Res") 

ce <- conditional_effects(BrmS1Hres0, effects = "Gen", conditions = data.frame(firstP = rep(c("T", "H"), each = 2),
                                                                      HostID = rep(c("M", "S"), times = 2)))

# Convert the list to a data frame
ce_df <- ce[[1]]

# Create the facet plot
S1CorrHres <- ggplot(ce_df, aes(x = Gen, y = estimate__)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  facet_grid(firstP ~ HostID) +
  labs(x = "Gen", y = "Estimated Res", title = "Conditional Effects of Gen on S1 Res") +
  theme_minimal()




####left2

print("Stage 2 Host abundance with ar1, by 2parasitoids*host, block with hostid")

BrmS2HAbun <- (brm((Habundance) ~  (firstP + secondP + control+ Gen)^3+ ar(p=1)+(1|Cage.)  +(1|HostID),
                       data = TMBs2NA, family = poisson))

saveRDS(BrmS2HAbun, file = "Model/BrmS2HAbun.rds")
BrmS2HAbun<- readRDS(file = "Model/BrmS2HAbun.rds")


sumBrmS2HAbun_tbl <-
  posterior_summary(BrmS2HAbun) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS2HAbun_tbl
plot(BrmS2HAbun)






print("Stage 2 Parasitoid abundance with ar1, by 2parasitoids")

BrmS2PAbunAR1 <- (brm((Pabundance) ~  (firstP +secondP+control +Gen)^3 +ar(p=1)+  (1|Cage.)+(1|HostID),
                       data = TMBs2NA, family = poisson))


saveRDS(BrmS2PAbunAR1, file = "Model/BrmS2PAbunAR1.rds")
BrmS2PAbunAR1<- readRDS(file = "Model/BrmS2PAbunAR1.rds")


sumBrmS2PAbunAR1_tbl <-
  posterior_summary(BrmS2PAbunAR1) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS2PAbunAR1_tbl
plot(BrmS2PAbunAR1)

print("Stage 2 Parasitoid abundance with ar2, by 2parasitoids")

BrmS2PAbunAR2 <- (brm((Pabundance) ~  (firstP +secondP+control +Gen)^3 +ar(p=2)+  (1|Cage.)+(1|HostID),
                   data = TMBs2NA, family = poisson))


saveRDS(BrmS2PAbunAR2, file = "Model/BrmS2PAbunAR2.rds")
BrmS2PAbunAR2<- readRDS(file = "Model/BrmS2PAbunAR2.rds")


sumBrmS2PAbunAR2_tbl <-
  posterior_summary(BrmS2PAbunAR2) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS2PAbunAR2_tbl
plot(BrmS2PAbunAR2)


#print("Stage 2 Host Sex rate with ar1, by 2 Ps")

#GlmS2SexAbun <- glmmTMB(cbind(Male.size, Female.size) ~ (firstP + secondP + control + Gen)^3 + ar1(GenID -1|cage) + (1|HostID), 
#                        data =na.omit(TMBs2NA), family = binomial)
#testZeroInflation(GlmS2SexAbun)
## no significance 









print("Stage 2 parasitism rate with Gen, Ar1, by P*H")

BrmS2PrateAR10 <- (brm(Pabundance | trials(totavalar) ~ 0 +  (firstP + secondP + control + mo(Gen) + HostID)^3 + ar(p=1)+
                        (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial, save_all_pars = T))

BrmS2Hres_no_first<-(brm(totdot | trials(Pabundance+totdot) ~ 0 +  (secondP + control + mo(Gen) + HostID)^3 + ar(p=1)+
                               (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial, save_all_pars = T))
bayes_factor(BrmS2Hres0,BrmS2PrateAR10)
BrmS2Hres_no_control<- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (firstP + secondP  + mo(Gen) + HostID)^3 + ar(p=1)+
       (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial, save_all_pars = T))
BrmS2Hres_no_secondP<- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (firstP + control  + mo(Gen) + HostID)^3 + ar(p=1)+
                                  (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
BrmS2Hres_no_Gen<-(brm(totdot | trials(Pabundance+totdot) ~ 0 +  (secondP + control  + HostID)^3 + ar(p=1)+
                               (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
BrmS2Hres_no_host <- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (firstP + secondP + control + mo(Gen) )^3 + ar(p=1)+
                         (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
##
BrmS2Hres02 <- brm(totdot | trials(Pabundance+totdot)~ 0 + (firstP+secondP+control+HostID+mo(Gen))^2 + 
                    ar(p=1)+(1|Cage.), data = TMBs2NA,cores =  4, family = binomial)

BrmS2Hres_no_first2<-(brm(totdot | trials(Pabundance+totdot) ~ 0 +  (secondP + control + mo(Gen) + HostID)^2 + ar(p=1)+
                           (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))

BrmS2Hres_no_contro2l<- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (firstP + secondP  + mo(Gen) + HostID)^2 + ar(p=1)+
                              (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
BrmS2Hres_no_secondP2<- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (firstP + control  + mo(Gen) + HostID)^2 + ar(p=1)+
                              (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
BrmS2Hres_no_Gen2<-(brm(totdot | trials(Pabundance+totdot) ~ 0 +  (secondP + control  + HostID)^2 +
                         (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
BrmS2Hres_no_host2 <- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (firstP + secondP + control + mo(Gen) )^2 + ar(p=1)+
                            (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
BrmS2HresNo122 <- (brm(totdot | trials(Pabundance+totdot) ~ 0 +  (control + mo(Gen) + HostID)^2 + ar(p=1)+
                         (1|Cage.) , data = TMBs2NA, cores = 4,family = binomial))
##
TMBs2NAG2<- subset(TMBs2NA, Gen == 2)

BrmS2HresG2 <- (brm(Pabundance | trials(totavalar) ~ 0 +  (firstP + secondP + control + HostID)^3 + 
                         (1|Cage.) , data = TMBs2NAG2, cores = 4,family = binomial))

BrmS2HresG2No1p <- (brm(Pabundance | trials(totavalar) ~ 0 +  (secondP + control + HostID)^3 +
                      (1|Cage.) , data = TMBs2NAG2, cores = 4,family = binomial))
BrmS2HresG2No2p <- (brm(Pabundance | trials(totavalar) ~ 0 +  (firstP + control + HostID)^3 +
                      (1|Cage.) , data = TMBs2NAG2, cores = 4,family = binomial))
BrmS2HresG2NoCon <- (brm(Pabundance | trials(totavalar) ~ 0 +  (firstP + secondP  + HostID)^3 +
                      (1|Cage.) , data = TMBs2NAG2, cores = 4,family = binomial))

TMBs2NA<- TMBs2NA %>%
  mutate(control = 1- control)
BrmS2Hres0 <- brm(totdot | trials(Pabundance+totdot)~ 0+ (firstP+secondP+control+HostID+mo(Gen))^3 + 
                    ar(p=1)+(1|Cage.), data = TMBs2NA,cores =  4, family = binomial)
hypothesis(BrmS2Hres0, "firstPT + secondPT + firstPT:secondPT = 0")
hypothesis(BrmS2Hres0, "firstPT = 0")

posterior_summary(BrmS2Hres0, pars = "^b_")
fixed <- (fixef(BrmS2Hres0))
rownames(fixed)

conditional_effects(BrmS2Hres0, effects = c("firstP:secondP"))

hypothesis(BrmS2Hres0, "firstPT:secondPT - firstPT:secondPH=0")

pp_check(BrmS2Hres02,type = "hist")

loo(BrmS2Hres02)

r2_full <- bayes_R2(BrmS2Hres02)
r2_nofirst <- bayes_R2(BrmS2Hres_no_first2)
r2_nocon <- bayes_R2(BrmS2Hres_no_contro2l)
r2_nosecond <- bayes_R2(BrmS2Hres_no_secondP2)
r2_noGen <- bayes_R2(BrmS2Hres_no_Gen2)
r2_noHost <- bayes_R2(BrmS2Hres_no_host2)
r2_no12 <- bayes_R2(BrmS2HresNo122)



r2_fullres <- bayes_R2(BrmS2HresG2)
r2_no1res <- bayes_R2(BrmS2HresG2No1p)
r2_no2res <- bayes_R2(BrmS2HresG2No2p)
r2_noconres <- bayes_R2(BrmS2HresG2NoCon)

delta_r2_firstP <- r2_full[,"Estimate"] - r2_nofirst[,"Estimate"]
delta_r2_con <- r2_full[,"Estimate"] - r2_nocon[,"Estimate"]



loo_R2(BrmS2PrateAR10)
cvfit <- cv_varsel(BrmS2PrateAR10)
draws1 <- as_draws_df(BrmS2Hresex)
draws <-
  draws1 %>% 
  select(starts_with("b")) 
draws[,16]<- draws[,16]+2
draws[,21]<- draws[,21]-1
draws[,22]<- draws[,22]+1
f <- function(x) {
  r <- quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}



melted_draws <- melt(draws)

# Calculate lower and upper quartiles for each variable

quartiles<- data.frame()

for (i in unique(melted_draws$variable) ) {
  test <- subset(melted_draws, variable == i)
  quartiles0 <- test %>%
    summarize(
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    mutate(
      asterisk = ifelse((lower > 0 | upper < 0), "*", ""),
      variable = i)
  quartiles <- rbind(quartiles,quartiles0)
}
quartiles$variable <- as.factor(quartiles$variable)

quartiles <- quartiles %>%
  mutate(is_control_interaction = ifelse(grepl("control", variable), "Population-size-mediated", "Population-size-unrelated"))
melted_draws <- melt(draws)
melted_draws<-melted_draws%>%
  left_join(select(quartiles, variable, is_control_interaction), by = "variable")
# Check the output of quartiles
print(quartiles)

melted_draws$is_control_interaction <- as.factor(melted_draws$is_control_interaction)


# Create the box plot and add asterisks
###2nd Stage Host Res fin


# Precompute the y order and shading positions

###reorder parameters 


new_order <- melted_draws %>%
  distinct(variable, is_control_interaction) %>%
  mutate(
    type_order = ifelse(is_control_interaction == "Fixed Population", 1, 2)
  ) %>%
  arrange(type_order) %>%
  pull(variable)

# 1. Set the factor level order (population-size-mediated on bottom)
melted_draws$variable <- factor(melted_draws$variable, levels = new_order)
title <- c(
  # Baseline (population-size-mediated) effects
  expression("1st Parasitoid" ~ italic("L. heterotoma")),
  expression("1st Parasitoid" ~ italic("A. tabida")),
  expression("2nd Parasitoid" ~ italic("A. tabida")),
  expression("Host" ~ italic("D. simulans")),
  expression("1st × 2nd Parasitoid" ~ italic("A. tabida")),
  expression("1st Parasitoid" ~ italic("A. tabida") * " × Host " * italic("D. simulans")),
  expression("2nd Parasitoid" ~ italic("A. tabida") * " × Host " * italic("D. simulans")),
  expression("1st × 2nd Parasitoid" ~ italic("A. tabida") * " × Host " * italic("D. simulans")),
  "Time",
  expression("Time × 1st Parasitoid" ~ italic("A. tabida")),
  expression("Time × 2nd Parasitoid" ~ italic("A. tabida")),
  expression("Time × Host " * italic("D. simulans")),
  expression("Time × 1st × 2nd Parasitoid" ~ italic("A. tabida")),
  expression("Time × 1st Parasitoid" ~ italic("A. tabida") * " × Host " * italic("D. simulans")),
  expression("Time × 2nd Parasitoid" ~ italic("A. tabida") * " × Host " * italic("D. simulans")),
  
  # Fixed population effects
  "Fixed Population",
  expression("Fixed Population × Host " * italic("D. simulans")),
  expression("1st Parasitoid" ~ italic("A. tabida") * " × Fixed Population"),
  expression("1st Parasitoid" ~ italic("A. tabida") * " × Fixed Population × Host " * italic("D. simulans")),
  expression("2nd Parasitoid" ~ italic("A. tabida") * " × Fixed Population"),
  expression("2nd Parasitoid" ~ italic("A. tabida") * " × Fixed Population × Host " * italic("D. simulans")),
  expression("1st × 2nd Parasitoid" ~ italic("A. tabida") * " × Fixed Population"),
  expression("Time × Fixed Population"),
  expression("Time × Fixed Population × Host " * italic("D. simulans")),
  expression("Time × 1st Parasitoid" ~ italic("A. tabida") * " × Fixed Population"),
  expression("Time × 2nd Parasitoid" ~ italic("A. tabida") * " × Fixed Population")
)


param_order <- c(
  # Baseline (population-size-mediated) effects
  "b_firstPH",                     # L. heterotoma
  "b_firstPT",                     # A. tabida 1st
  "b_secondPT",                    # A. tabida 2nd
  "b_HostIDS",
  "b_firstPT:secondPT",
  "b_firstPT:HostIDS",
  "b_secondPT:HostIDS",
  "b_firstPT:secondPT:HostIDS",
  "bsp_moGen",
  "bsp_moGen:firstPT",
  "bsp_moGen:secondPT",
  "bsp_moGen:HostIDS",
  "bsp_moGen:firstPT:secondPT",
  "bsp_moGen:firstPT:HostIDS",
  "bsp_moGen:secondPT:HostIDS",
  
  # Fixed population effects
  "b_control",
  "b_control:HostIDS",
  "b_firstPT:control",
  "b_firstPT:control:HostIDS",
  "b_secondPT:control",
  "b_secondPT:control:HostIDS",
  "b_firstPT:secondPT:control",
  "bsp_moGen:control",
  "bsp_moGen:control:HostIDS",
  "bsp_moGen:firstPT:control",
  "bsp_moGen:secondPT:control"
)


melted_draws$variable <- factor(melted_draws$variable, levels = rev(param_order))
quartiles$variable <- factor(quartiles$variable, levels = rev(param_order))




## final 2nd Stage Host Resistance 
ggplot(data = melted_draws, aes(x = value, y = variable)) +
  stat_summary(aes(color = is_control_interaction), fun = median, geom = "point", size = 3) +  
  stat_summary(aes(color = is_control_interaction), 
               fun.min = ~quantile(.x, 0.025), 
               fun.max = ~quantile(.x, 0.975), 
               geom = "errorbarh", height = 0.2, linewidth = 1) +  
  scale_y_discrete(labels = rev(title)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("2nd Stage Host Resistance") +
  geom_text(data = quartiles, aes(x = max(melted_draws$value) * 0.95, y = variable, label = asterisk), 
            size = 5, color = "red", hjust = 1) + 
  scale_color_manual(name = "Population Size Effect", 
                     values = c("Population-size-mediated" = "red", "Population-size-unrelated" = "black")) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 9))





##
conditions <- expand.grid(firstP = c("T", "H"), HostID = c("S", "M"), secondP = c("T", "H"), control = c(0, 1))
ce <- conditional_effects(BrmS2PrateAR10, effects = "Gen", conditions = conditions)

# Convert the list to a data frame
ce_df <- ce[[1]]

# Combine A and D into a single column for line identity
ce_df$fsP <- interaction(ce_df$firstP, ce_df$secondP, sep = "")

# Create the facet plot
S2corrPrate <- ggplot(ce_df, aes(x = Gen, y = estimate__, color = fsP)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = fsP), alpha = 0.2, color = NA) +
  facet_grid(HostID ~ control, labeller = label_both) +
  labs(x = "Gen", y = "Estimated Prate", title = "Conditional Effects of Gen on Prate") +
  theme_minimal() +
  scale_color_manual(values = c("TT" = "blue", "TH" = "orange", "HT" = "lightgreen", "HH" = "purple")) +
  scale_fill_manual(values = c("TT" = "blue", "TH" = "orange", "HT" = "lightgreen", "HH" = "purple"))

saveRDS(BrmS2PrateAR10, file = "Model/BrmS2PrateAR10.rds")
BrmS2PrateAR1<- readRDS(file = "Model/BrmS2PrateAR1.rds")


sumBrmS2PrateAR1_tbl <-
  posterior_summary(BrmS2PrateAR1) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS2PrateAR1_tbl
plot(BrmS2PrateAR1)

print("Stage 2 parasitism rate with Gen, Ar2, by P*H")

prior1 <- prior(normal(0, 1), class = ar, lb = -1, ub = 1)
prior2 <- prior(normal(0, 0.4), class = ar, lb = -1, ub = 1)


BrmS2PrateAR2 <- (brm(Pabundance | trials(totavalar) ~ (firstP + secondP + control + Gen)^3 + ar(p=2)+(1|Cage.) + 
                         (1|HostID), data = TMBs2NA,prior = prior1, family = binomial, cores = 4))



BrmS2PrateAR2Prior2 <- (brm(Pabundance | trials(totavalar) ~ (firstP + secondP + control + Gen)^3 + ar(p=2, gr = Cage.)+(1|Cage.) + 
                        (1|HostID), data = TMBs2NA,prior = prior2, family = binomial, cores = 4))

plot(BrmS2PrateAR2Prior2)
mcmc_pairs(BrmS2PrateAR2,  c("ar[1]","ar[2]"))





saveRDS(BrmS2PrateAR2, file = "Model/BrmS2PrateAR2.rds")
BrmS2PrateAR2<- readRDS(file = "Model/BrmS2PrateAR2.rds")


sumBrmS2PrateAR2_tbl <-
  posterior_summary(BrmS2PrateAR2) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

sumBrmS2PrateAR2_tbl
plot(BrmS2PrateAR2)

mcmc_pairs(BrmS2PrateAR2,  c("ar[1]","ar[2]"))
prior_summary(BrmS2PrateAR2)

print("Stage 2 tot resistance rate by sex with gen by P*H*control")
TMBs2sex <- subset(TMBs2sex, sex == "M")

TMBs2sex<- TMBs2sex %>%
  mutate(control = 1- control)
BrmS2Hresex <- brm(totdot | trials(Pabundance+totdot)~ 0 + (firstP+secondP+control+HostID+mo(Gen))^3 + 
                    ar(p=1)+(1|Cage.), data = TMBs2sex,cores =  4, family = binomial)
write_rds(BrmS2Hresex, file = "/Users/Apple/Desktop/PhD/Code/BrmS2Hresex.rds")

BrmS2Hresex<- read_rds(file = "/Users/Apple/Desktop/PhD/Code/BrmS2Hresex.rds")

BrmS2HresexAR2 <- brm(totdot | trials(Pabundance+totdot)~ 0 + (firstP+secondP+control+HostID+mo(Gen))^3 + 
                     ar(p=2)+(1|Cage.), data = TMBs2sex,cores =  4, family = binomial)
BrmS2HresAR2Prior2 <- (brm(totdot | trials(Pabundance+totdot) ~ (firstP + secondP + control + Gen)^3 + ar(p=2, gr = Cage.)+(1|Cage.) + 
                              (1|HostID), data = TMBs2sex,prior = prior2, family = binomial, cores = 4))


BrmS2Hres0 <- brm(totdot | trials(Pabundance+totdot)~ 0 + (firstP+secondP+control+HostID+mo(Gen))^3 + 
                    ar(p=1)+(1|Cage.), data = TMBs2NA,cores =  4, family = binomial)
library(reshape2)


draws1 <- as_draws_df(BrmS2Hresex)
draws <-
  draws1 %>% 
  select(starts_with("b")) 
draws[,16:26] <- draws[,16:26] * 7


f <- function(x) {
  r <- quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}



melted_draws <- melt(draws)

# Calculate lower and upper quartiles for each variable

quartiles<- data.frame()

for (i in unique(melted_draws$variable) ) {
  test <- subset(melted_draws, variable == i)
  quartiles0 <- test %>%
    summarize(
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    mutate(
      asterisk = ifelse((lower > 0 | upper < 0), "*", ""),
      variable = i)
  quartiles <- rbind(quartiles,quartiles0)
}


# Check the output of quartiles


# Create the box plot and add asterisks

title <- c("1st Parasitoid L. heterotoma", "1st Parasitoid A. tabida","2nd Parasitoid A. tabida",
           "Fixed Population","Host D.simulans", "1st Parasitoid A.tabida:2nd Parasitoid A.tabida", "1st Parasitoid A.tabida:Fixed Population",
           "1st Parasitoid A.tabida:Host D.simulans","2nd Parasitoid A.tabida:Fixed Population","2nd Parasitoid A.tabida:Host D.simulans",
           "Fixed Population:Host D. simulans", "1st Parasitoid A.tabida:2nd Parasitoid A.tabida:Fixed Population", "1st Parasitoid A.tabida:2nd Parasitoid A.tabida:Host D.simulans",
           "1st Parasitoid A.tabida:Fixed Population:Host D.simulans","2nd Parasitoid A.tabida:Fixed Population:Host D.simulans",
           "Time","Time:1st Parasitoid A.tabida","Time:2nd Parasitoid A.tabida","Time:Fixed Population","Time:Host D.simulans","Time:1st Parasitoid A.tabida:2nd Parasitoid A.tabida",
           "Time:1st Parasitoid A.tabida:Fixed Population","Time:1st Parasitoid A.tabida:Host D.simulans","Time:2nd Parasitoid A.tabida:Fixed Population","Time:2nd Parasitoid A.tabida:Host D.simulans",
           "Time:Fixed Population:Host D.simulans") 




###variance partition


conditions <- expand.grid(
  firstP = c("T", "H"),
  HostID = c("S", "M"),
  secondP = c("T", "H"),
  control = c(0, 1)
)
ce <- conditional_effects(BrmS2Hresex, effects = "Gen", conditions = conditions)
ce_df <- ce[[1]]  # Get the data frame
ce_df$Parasitoid_Identity_across_Stages <- interaction(ce_df$firstP, ce_df$secondP, sep = "")
ce_df$GroupType <- ifelse(ce_df$Parasitoid_Identity_across_Stages %in% c("HH", "TT"), "Same", "Different")
ce_df$Parasitoid_Identity_across_Stages <- interaction(ce_df$firstP, ce_df$secondP, sep = "")
ce_df$Parasitoid_Label <- gsub("(\\w)(\\w)", "\\1\n\\2", ce_df$Parasitoid_Identity_across_Stages)



# Convert the list to a data frame
#ce_df <- ce[[1]]

# Combine A and D into a single column for line identity
ce_df$Parasitoid_ID_2_Stages <- interaction(ce_df$firstP, ce_df$secondP, sep = "")
ce_df$estimate__



ce_df$estimate__ <- scale_values(ce_df$estimate__)
ce_df$upper__ <- scale_values(ce_df$upper__)
ce_df$lower__ <- scale_values(ce_df$lower__)



main_line_plot <- function(data_sub) {
  ggplot(data_sub, aes(x = Gen, y = estimate__, 
                       color = Parasitoid_Identity_across_Stages,
                       linetype = Parasitoid_Identity_across_Stages)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, 
                    fill = Parasitoid_Identity_across_Stages), alpha = 0.1, linewidth = 0,color =) +
    scale_linetype_manual(values = c("TT" = "dashed", "HH" = "dashed", "HT" = "solid", "TH" = "solid")) +
    scale_fill_manual(values = c("TT" = "blue", "TH" = "orange", "HT" = "lightgreen", "HH" = "purple")) +
    scale_color_manual(values = c("TT" = "blue", "TH" = "orange", "HT" = "lightgreen", "HH" = "purple")) +
    labs(x = "Time", y = NULL,
         shape = "Parasitoid\nCombination",
         color = "First and second\nparasitoid species",
         fill = "First and second\nparasitoid species",
         linetype = "First and second\nparasitoid species") +
    theme_minimal() +
    ylim(0, 1)
}
boxplot_side <- function(data, gen_value, show_ylabel = FALSE, show_shape_legend = TRUE) {
  plot_data <- data %>% 
    filter(Gen == gen_value) %>%
    mutate(
      GroupType = ifelse(Parasitoid_Identity_across_Stages %in% c("HH", "TT"), "Conspecific", "Interspecific"),
      Shape = ifelse(GroupType == "Conspecific", 1, 16)  # 1 = hollow, 16 = filled
    )
  
  ggplot(plot_data, aes(x = Parasitoid_Identity_across_Stages, y = estimate__,
                        color = Parasitoid_Identity_across_Stages,
                        shape = GroupType)) +
    geom_pointrange(aes(ymin = lower__, ymax = upper__), 
                    fatten = 3, size = 0.7) +
    scale_color_manual(values = c("TT" = "blue", "TH" = "orange", "HT" = "coral2", "HH" = "purple")) +
    scale_shape_manual(values = c("Conspecific" = 1, "Interspecific" = 16)) +
    labs(
      x = NULL,
      y = if (show_ylabel) "Model Estimate" else NULL,
      title = paste("Gen", gen_value),
      shape = "Parasitoid\nCombination",
      color = "First and second\nparasitoid species",
      fill = "First and second\nparasitoid species",
      linetype = "First and second\nparasitoid species"
    )+
    ylim(0, 1) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = if (show_ylabel) element_text() else element_blank(),
      axis.title.y = if (show_ylabel) element_text() else element_blank(),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = if (show_shape_legend) "right" else "none"
    )
}


panel_combos <- expand.grid(HostID = c("S", "M"), control = c(0, 1))
panel_plots <- list()

for (i in 1:nrow(panel_combos)) {
  host <- panel_combos$HostID[i]
  ctrl <- panel_combos$control[i]
  
  data_sub <- ce_df %>% filter(HostID == host, control == ctrl)
  show_ylabel <- ctrl == 0
  
  p1 <- boxplot_side(data_sub, gen_value = 2, show_ylabel = show_ylabel)
  p2 <- main_line_plot(data_sub)
  p3 <- boxplot_side(data_sub, gen_value = 7, show_ylabel = FALSE)
  
  panel <- p1  + p2+p3 + patchwork::plot_layout(ncol = 3, widths = c(1, 2, 1))
  panel_plots[[i]] <- panel
}
final_plot <- wrap_plots(panel_plots, ncol = 2, guides = "collect") &
  theme(legend.position = "right")

# Add centered y-axis label
ylabel_grob <- grid::textGrob("Model Estimate", rot = 90, gp = gpar(fontsize = 12))
final_plot <- final_plot + inset_element(ylabel_grob, left = -0.04, bottom = 0, right = 0, top = 1)


title_labels <- list(
  wrap_elements(ggplot() + 
                  theme_void() + 
                  annotate("text", x = 0.5, y = 0.5, 
                           label = expression(italic("D. melanogaster")~":Control"), 
                           size = 5, fontface = "bold")),
  
  wrap_elements(ggplot() + 
                  theme_void() + 
                  annotate("text", x = 0.5, y = 0.5, 
                           label = expression(italic("D. melanogaster")~":Fixed Population"), 
                           size = 5, fontface = "bold")),
  
  wrap_elements(ggplot() + 
                  theme_void() + 
                  annotate("text", x = 0.5, y = 0.5, 
                           label = expression(italic("D. simulans")~":Control"), 
                           size = 5, fontface = "bold")),
  
  wrap_elements(ggplot() + 
                  theme_void() + 
                  annotate("text", x = 0.5, y = 0.5, 
                           label = expression(italic("D. simulans")~":Fixed Population"), 
                           size = 5, fontface = "bold"))
)


panels_with_titles <- list()

for (i in 1:4) {
  titled_panel <- title_labels[[i]] / panel_plots[[i]] + 
    plot_layout(heights = c(0.1, 1))  # Adjust spacing if needed
  panels_with_titles[[i]] <- titled_panel
}


final_plot <- wrap_plots(panels_with_titles, ncol = 2, guides = "collect") &
  theme(legend.position = "right")

# Add centered y-axis label on the left
ylabel_grob <- grid::textGrob("Model Estimate", rot = 90, gp = gpar(fontsize = 12))
final_plot <- final_plot + inset_element(ylabel_grob, left = -0.04, bottom = 0, right = 0, top = 1)





# Custom plotting function
boxplot_split_by_treatment <- function(data, gen_value, show_ylabel = FALSE) {
  plot_data <- data %>%
    filter(Gen == gen_value) %>%
    mutate(
      Treatment = ifelse(control == 0, "Control", "Fixed Population"),
      GroupSide = ifelse(Parasitoid_Identity_across_Stages %in% c("HH", "TT"), "Left", "Right"),
      Shape = ifelse(GroupSide == "Left", "Conspecific", "Interspecific"),
      XPos = paste(Treatment, Parasitoid_Identity_across_Stages, sep = "_")
    )
  
  if (nrow(plot_data) == 0) return(ggplot() + theme_void())
  
  plot_data$XPos <- factor(plot_data$XPos,
                           levels = c("Control_TT", "Control_HH", "Control_HT", "Control_TH",
                                      "Fixed Population_TT", "Fixed Population_HH", "Fixed Population_HT", "Fixed Population_TH"))
  
  ggplot(plot_data, aes(x = XPos, y = estimate__,
                        color = Parasitoid_Identity_across_Stages,
                        shape = Shape)) +
    geom_pointrange(aes(ymin = lower__, ymax = upper__), fatten = 3, size = 0.7) +
    scale_color_manual(
      values = c("TT" = "blue", "TH" = "orange", "HT" = "coral2", "HH" = "purple"),
      labels = c(
        "TT" = expression(italic("A. tab") ~ "\u2192" ~ italic("A. tab")),
        "TH" = expression(italic("A. tab") ~ "\u2192" ~ italic("L. het")),
        "HT" = expression(italic("L. het") ~ "\u2192" ~ italic("A. tab")),
        "HH" = expression(italic("L. het") ~ "\u2192" ~ italic("L. het"))
      )
    ) +
    scale_shape_manual(values = c("Conspecific" = 1, "Interspecific" = 16)) +
    ylim(-0.1, 1) +
    labs(
      x = NULL,
      y = if (show_ylabel) "Model Estimate Host Resistance Rate" else NULL,
      shape = "Parasitoid Combination",
      color = "First and second\nparasitoid species"
    )+
    # Labels at the bottom row
    annotate("text", x = 2, y = -0.05, label = "Control", size = 3.5, fontface = "plain") +
    annotate("text", x = 6, y = -0.05, label = "Fixed Population", size = 3.5, fontface = "plain") +
    # Visual divider line
    geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey50") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = if (show_ylabel) element_text() else element_blank(),
      axis.title.y = if (show_ylabel) element_text() else element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
}


# Create all four panels
gens <- c(1, 7)
species <- c("S", "M")
panels_grid <- list()

for (i in seq_along(gens)) {
  g <- gens[i]
  row_panels <- list()
  for (j in seq_along(species)) {
    s <- species[j]
    data_sub <- ce_df %>% filter(HostID == s)
    show_ylabel <- (j == 1)
    
    panel <- boxplot_split_by_treatment(data_sub, gen_value = g, show_ylabel = show_ylabel)
    
    # Add italic species label to top panels only
    if (i == 1) {
      species_label <- if (s == "S") "D. simulans" else "D. melanogaster"
      title_plot <- ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = deparse(bquote(italic(.(species_label)))), parse = TRUE,
                 size = 5, fontface = "bold")
      panel <- wrap_elements(title_plot) / panel + plot_layout(heights = c(0.12, 1))
    }
    
    row_panels[[j]] <- panel
  }
  panels_grid[[i]] <- wrap_plots(row_panels, ncol = 2)
}

# Combine final layout
final_plot <- wrap_plots(panels_grid, nrow = 2, guides = "collect") &
  theme(legend.position = "right")

# Add central vertical y-axis label
ylabel_grob <- grid::textGrob("Model Estimate", rot = 90, gp = gpar(fontsize = 12))
final_plot <- final_plot + inset_element(ylabel_grob, left = -0.04, bottom = 0, right = 0, top = 1)

# Add generation labels on the right side
gen_labels <- list(
  wrap_elements(ggplot() + theme_void() +
                  annotate("text", x = 0.1, y = 0.1, label = "Gen 2", size = 5)),
  wrap_elements(ggplot() + theme_void() +
                  annotate("text", x = 0.1, y = 0.1, label = "Gen 7", size = 5))
)

for (i in 1:2) {
  panels_grid[[i]] <- wrap_plots(panels_grid[[i]], gen_labels[[i]], widths = c(1, 0.08))
}

final_plot <- wrap_plots(panels_grid, nrow = 2, guides = "collect") &
  theme(legend.position = "right")


# Display final result
final_plot








###making contrast figure 

delta_df <- ce_df %>%
  select(Gen, HostID, control, Parasitoid_Identity_across_Stages,
         estimate__, lower__, upper__,se__) %>%
  pivot_wider(names_from = control, values_from = c(estimate__, lower__, upper__,se__),
              names_prefix = "ctrl") %>%
  mutate(
    delta_estimate = estimate___ctrl0 - estimate___ctrl1,
    delta_lower = lower___ctrl0 - lower___ctrl1,
    delta_upper = upper___ctrl0 - upper___ctrl1,
    delta_se = sqrt((se___ctrl0)^2) + sqrt((se___ctrl1)^2)
  )


#delta_df <- delta_df %>%
#  mutate(delta_upper = delta_estimate + delta_se,
#         delta_lower = delta_estimate - delta_se)


plot_delta_lines <- function(host_id, host_label) {
  data_sub <- delta_df %>% filter(HostID == host_id)
  
  ggplot(data_sub, aes(x = Gen, y = delta_estimate, 
                       color = Parasitoid_Identity_across_Stages,
                       linetype = Parasitoid_Identity_across_Stages)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = delta_lower, ymax = delta_upper, 
                    fill = Parasitoid_Identity_across_Stages), alpha = 0.2, linewidth = 0) +
    labs(
      title = host_label,
      x = "Generation",
      y = expression(Delta~"Model Estimate"),
      color = "First and second\nparasitoid species",
      fill = "First and second\nparasitoid species",
      linetype = "First and second\nparasitoid species"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("TT" = "blue", "TH" = "orange", 
                                  "HT" = "lightgreen", "HH" = "purple")) +
    scale_linetype_manual(values = c("TT" = "dotted", "HH" = "dotted", 
                                     "HT" = "solid", "TH" = "solid")) +
        scale_fill_manual(values = c("TT" = "blue", "TH" = "orange", 
                                 "HT" = "lightgreen", "HH" = "purple")) +
    ylim(-1, 1)  # Adjust if needed
}



p_mel <- plot_delta_lines("M", expression(italic("D. melanogaster")))
p_sim <- plot_delta_lines("S", expression(italic("D. simulans")))

delta_plot <- p_mel / p_sim + 
  patchwork::plot_annotation(
    title = "Contrast Between Control and Fixed Populations",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  ) & 
  theme(legend.position = "bottom")

delta_plot



saveRDS(BrmS2Hres0, file = "Model/BrmS2Hres0.rds")
BrmS2Hres<- readRDS(file = "Model/BrmS2Hres0.rds")

BrmS2PrateAR10
BrmS2PrateAR2
BrmS2PrateAR2Prior2

bayes_R2(BrmS2PrateAR10)
loo(BrmS2PrateAR2)


TMBs2NA0  <- TMBs2NA

TMBs2NA0$Pabundance[(TMBs2NA0$Pabundance) == 0] <- NA
TMBs2NA0$totdot[TMBs2NA0$totdot == 0] <- NA

BrmS2Hres0s2 <- brm(totdot | trials(Pabundance+totdot)~ 0+ (firstP+secondP+control+HostID  +mo(Gen))^3 + 
                         ar(p=1)+(1|Cage.), data = TMBs2NA0,cores =  4, family = binomial)


draws1 <- as_draws_df(BrmS2Hres0s2)
draws <-
  draws1 %>% 
  select(starts_with("b")) 
draws[,16:26] <- draws[,16:26] * 7


f <- function(x) {
  r <- quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}



melted_draws <- melt(draws)

# Calculate lower and upper quartiles for each variable

quartiles<- data.frame()

for (i in unique(melted_draws$variable) ) {
  test <- subset(melted_draws, variable == i)
  quartiles0 <- test %>%
    summarize(
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    mutate(
      asterisk = ifelse((lower > 0 | upper < 0), "*", ""),
      variable = i)
  quartiles <- rbind(quartiles,quartiles0)
}
quartiles$variable <- as.factor(quartiles$variable)

quartiles <- quartiles %>%
  mutate(is_control_interaction = ifelse(grepl("control", variable), "Population-size-mediated", "Fixed Population"))
melted_draws <- melt(draws)
melted_draws<-melted_draws%>%
  left_join(select(quartiles, variable, is_control_interaction), by = "variable")
# Check the output of quartiles
print(quartiles)

melted_draws$is_control_interaction <- as.factor(melted_draws$is_control_interaction)


ggplot(data = melted_draws, aes(x = value, y = variable)) +
  # Correct background position
  geom_rect(data = melted_draws, aes(fill = is_control_interaction),
            xmin = -Inf, xmax = Inf, 
            ymin = match(melted_draws$variable, rev(unique(melted_draws$variable))) - 0.5,  
            ymax = match(melted_draws$variable, rev(unique(melted_draws$variable))) + 0.5,  
            alpha = 0.2) +  # Adjust transparency
  
  # Median point (black)
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +  
  
  # Error bar from 25% to 97.5% (black)
  stat_summary(fun.min = function(x) quantile(x, 0.25), 
               fun.max = function(x) quantile(x, 0.975), 
               geom = "errorbarh", height = 0.2, linewidth = 1, color = "black") +  
  
  scale_y_discrete(limits = rev(unique(melted_draws$variable)), labels = rev(title)) + 
  theme_bw() +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Reference line at 0
  ggtitle("2nd Stage Parasitoid Parasitism") +
  
  # Asterisk for significance
  geom_text(data = quartiles, aes(x = 5, y = variable, label = asterisk), 
            size = 5, color = "red") +  
  
  # Define colors for background with a legend
  scale_fill_manual(name = "Ecological Category",  # Legend title
                    values = c("Ecological" = "pink", "Non-Ecological" = "lightgray")) +  
  
  theme(legend.position = "right")  # Show legend


S2Hres <- ggplot(data = melted_draws, aes(x = value, y = variable, fill = is_control_interaction)) +
  geom_violin(alpha = 0.3) +  # Optional: Shows density distribution
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +  # Median as a point
  stat_summary(fun.min = function(x) quantile(x, 0.25), 
               fun.max = function(x) quantile(x, 0.975), 
               geom = "errorbarh", height = 0.2, color = "black", linewidth = 1) +  # Error bar from 25% to 97.5%
  scale_y_discrete(limits = rev, labels = rev(title)) + 
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +  # Reference line at 0
  ggtitle("2nd Stage Host Resistance Rate") +
  geom_text(data = quartiles, aes(x = 5, y = variable, label = asterisk), 
            size = 5, color = "red") +  # Asterisk for significance
  scale_fill_manual(values = c("Control Interaction" = "red", "Other" = "black")) +
  theme(legend.position = "right")  # Show legend
pp_check(BrmS2HresOisNA1)

title <- c("1st Parasitoid H. leptopelina", "1st Parasitoid A. tabida","2nd Parasitoid A. tabida",
           "Fixed Population","Host D. simulans", "1st Parasitoid:2nd Parasitoid A. tabida", "1st Parasitoid A.tabida:Fixed Population",
           "1st Parasitoid A.tabida:Host Identity D. simulans","2nd Parasitoid A.tabida:Fixed Population","2nd Parasitoid A.tabida:Host Identity D. simulans",
           "Fixed Population:Host Identity D. simulans", "1st:2nd Parasitoid A.tabida:Fixed Population", "1st:2nd Parasitoid A.tabida:Host Identity D.simulans",
           "1st Parasitoid A.tabida:Fixed Population:Host Identity D.simulans","2nd Parasitoid A.tabida:Fixed Population:Host Identity D.simulans",
           "Time","Time:1st Parasitoid A.tabida","Time:2nd Parasitoid A.tabida","Time:Fixed Population","Time:Host Identity D.simulans","Time:1st:2nd Parasitoid A.tabida",
           "Time:1st Parasitoid A.tabida:Fixed Population","Time:1st Parasitoid A.tabida:Host Identity D.simulans","Time:2nd Parasitoid A.tabida:Fixed Population","Time:2nd Parasitoid A.tabida:Host Identity D.simulans",
           "Time:Fixed Population:Host Identity D.simulans") 



S2res <- ggplot(data = melt(draws), aes(x=value, y=variable)) + geom_boxplot(width = 0.01, outlier.shape = NA,fill = "black")+     scale_y_discrete(limits=rev)+
  stat_summary(fun.data = f, geom="boxplot")+theme_bw() + geom_vline(xintercept = 0) + ggtitle("2nd Stage Res")

conditions <- expand.grid(firstP = c("T", "H"),HostID = c("M","S"), secondP = c("T", "H"), control = c(0, 1))
ce <- conditional_effects(BrmS2Hres0, effects = "Gen", conditions = conditions)

# Convert the list to a data frame
ce_df <- ce[[1]]
ce_df$control[ce_df$control == 0]<- "Unconstrained"
ce_df$control[ce_df$control == 1]<-"Control"
# Combine A and D into a single column for line identity
ce_df$Parasitoid_Identity_across_Stages <- interaction(ce_df$firstP, ce_df$secondP, sep = "")
ce_df$estimate__


ce_df$estimate__ <- scale_values(ce_df$estimate__)
ce_df$upper__ <- scale_values(ce_df$upper__)
ce_df$lower__ <- scale_values(ce_df$lower__)

host_labels <- c("M" = "D. melanogaster", "S" = "D. simulans")
control_labels <- c("0" = "Unconstrained Populations", "1" = "Control")

# Create the facet plot
#To present on May 9 

#Host Abundance First stage 
sumBrmS1HAbun_tbl
plot(BrmS1HAbun)


#Para abundance AR1 first stage 
sumBrmS1PAbun_tbl
plot(BrmS1PAbunAR1)

#Para abundance AR2 firs stage 
sumBrmS1PAbunAR2_tbl
plot(BrmS1PAbunAR2)



#Parasitism rate AR1 S1
sumBrmS1PrateAR1_tbl
plot(BrmS1PrateAR1)


postS1Prate <- mcmc_areas(
  BrmS1Prate,
  regex_pars = "b_",
  prob = 0.89,
  prob_outer = 1,
  point_est = "median",
  area_method = "equal height"
) +
  geom_vline(xintercept = 0, color = "red", alpha = 0.6, lwd = .8, linetype = "dashed") +
  labs(
    title = "Stage 1 Parasitism rates prediction"
  )


#Parasitisms rate AR2 S1
sumBrmS1PrateAR2_tbl
plot(BrmS1PrateAR2)


#Host res First stage 
sumBrmS1Hres_tbl
plot(BrmS1Hres)

PostS1Hres <- mcmc_areas(
  BrmS1Hres,
  regex_pars = "b_",
  prob = 0.89,
  prob_outer = 1,
  point_est = "median",
  area_method = "equal height"
) +
  geom_vline(xintercept = 0, color = "red", alpha = 0.6, lwd = .8, linetype = "dashed") +
  labs(
    title = "Stage 1 resistance rates prediction"
  )


PostS1Hresmo<- mcmc_areas(
  BrmS2Hres,
  regex_pars = "b_",
  prob = 0.89,
  prob_outer = 1,
  point_est = "median",
  area_method = "equal height"
) +
  geom_vline(xintercept = 0, color = "red", alpha = 0.6, lwd = .8, linetype = "dashed") +
  labs(
    title = "Stage 1 resistance rates prediction"
  )

mcmc_areas()

#Host Abundance stage 2
sumBrmS2HAbun_tbl
plot(BrmS2HAbun)

#Para Abundance stage 2 Ar1
sumBrmS2PAbunAR1_tbl
plot(BrmS2PAbunAR1)

#para Abundance stage 2 Ar2
sumBrmS2PAbunAR2_tbl
plot(BrmS2PAbunAR2)
BrmS2PAbunAR2$fit

#Parasitism rate AR1 stage 2
sumBrmS2PrateAR1_tbl
plot(BrmS2PrateAR1)

#Parasitism rate AR2 stage 2 
sumBrmS2PrateAR2_tbl
plot(BrmS2PrateAR2)


#Res rate stage 2 AR1
sumBrmS2Hres_tbl
plot(BrmS2Hres)

summary(BrmS2Hres)


