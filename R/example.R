library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(sampling)
source('R/models.r')


### Treating survey data as population
pums <- read_csv('Data/ss14pmn.csv')
pums <- pums %>% dplyr::select(PUMA, PWGTP, SEX, RAC1P, POVPIP, HICOV) %>% 
  remove_missing() %>%
  mutate(PWGTP=as.numeric(PWGTP), POVPIP=as.numeric(POVPIP), HICOV=abs(HICOV-2), SEX=factor(SEX), RACE=factor(RAC1P)) %>% 
  mutate(IPR=case_when(
    POVPIP <= 138 ~ 1,
    POVPIP > 138 & POVPIP <= 200 ~ 2, 
    POVPIP > 200 & POVPIP <= 250 ~ 3,
    POVPIP > 250 & POVPIP <= 400 ~ 4, 
    POVPIP > 400 ~ 5, 
  ) %>% factor())
pums$resp <- factor(paste(pums$IPR, pums$HICOV, sep='_'))


pcells <- pums %>% group_by(PUMA, SEX, RACE) %>% summarize(popsize=n())
predX <- model.matrix(~ SEX + RACE - 1, data=pcells)
predPsi <- model.matrix(~as.factor(pcells$PUMA))


### Subsample population 
set.seed(1)
ss <- 10000
prob <- inclusionprobabilities(pums$PWGTP*(1 + (pums$HICOV=='1')), ss)
prob <- prob*ss/sum(prob)
Ind <- UPpoisson(prob)
sample <- pums[as.logical(Ind),]
sample$P <- prob[as.logical(Ind)]
sample$W <- 1/ sample$P
sample$scaledWGT <- (sample$W)*nrow(sample)/sum(sample$W)

### Fit model
modX <- model.matrix(~ SEX + RACE - 1, data=sample)
modPsi <- model.matrix(~ factor(PUMA, levels=levels(as.factor(pcells$PUMA))) - 1, data=sample)
modY <- model.matrix(~ resp -1, data=sample)
modY <- modY[,order(colSums(modY), decreasing = T)]

modMCMC <- PLMMmcmc(X=modX, Psi=modPsi, Y=modY, wgt=sample$scaledWGT, predX=predX, predPsi=predPsi, iter=100, burn=50) 
## note that more iterations should be used in practice...this is just a quick example

modVB <- PLMMvb(X=modX, Psi=modPsi, Y=modY, wgt=sample$scaledWGT, predX=predX, predPsi=predPsi, eps=0.01)
