library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(sampling)
library(ggthemes)
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

truth <- pums %>% group_by(PUMA, IPR) %>% summarize(P=mean(HICOV))

pcells <- pums %>% group_by(PUMA, SEX, RACE) %>% summarize(popsize=n())
predX <- model.matrix(~ SEX + RACE - 1, data=pcells)
predPsi <- model.matrix(~as.factor(pcells$PUMA)) ## Using PUMA incidence vectors


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


### Construct synthetic population and aggregate predictions
predsMCMC <- array(NA, dim=c(length(unique(pums$PUMA)), 5, dim(modMCMC$Preds)[3]))
for(j in 1:dim(modMCMC$Preds)[3]){
  temp <- data.frame(PUMA=pcells$PUMA, N=pcells$popsize, P=I(modMCMC$Preds[,,j]))
  temp2 <- t(apply(cbind(temp$N, temp$P), 1, function(x) rmultinom(1,x[1],x[-1])))
  temp2 <- data.frame(PUMA=pcells$PUMA, temp2) %>% group_by(PUMA) %>% summarize_all(sum)
  temp3 <- as.matrix(temp2[,c(4,5,6,3,2)]/(temp2[,c(7,8,11,9,10)]+temp2[,c(4,5,6,3,2)]))
  predsMCMC[,,j] <- temp3
}
plotMCMC <- data.frame(PUMA=rep(unique(pcells$PUMA), 5), IPRCAT=sort(rep(1:5, 43)), matrix(predsMCMC, ncol=dim(predsMCMC)[3])) %>% 
  mutate(Type="MCMC") %>%
  relocate(Type) %>%
  pivot_longer(-c(1:3), values_to="Est.")

predsVB <- array(NA, dim=c(length(unique(pums$PUMA)), 5, dim(modVB$Preds)[3]))
for(j in 1:dim(modVB$Preds)[3]){
  temp <- data.frame(PUMA=pcells$PUMA, N=pcells$popsize, P=I(modVB$Preds[,,j]))
  temp2 <- t(apply(cbind(temp$N, temp$P), 1, function(x) rmultinom(1,x[1],x[-1])))
  temp2 <- data.frame(PUMA=pcells$PUMA, temp2) %>% group_by(PUMA) %>% summarize_all(sum)
  temp3 <- as.matrix(temp2[,c(4,5,6,3,2)]/(temp2[,c(7,8,11,9,10)]+temp2[,c(4,5,6,3,2)]))
  predsVB[,,j] <- temp3
}
plotVB <- data.frame(PUMA=rep(unique(pcells$PUMA), 5), IPRCAT=sort(rep(1:5, 43)), matrix(predsVB, ncol=dim(predsVB)[3])) %>%
  mutate(Type="VB") %>% 
  relocate(Type) %>%
  pivot_longer(-c(1:3), values_to="Est.")

plotFinal <- rbind(plotMCMC, plotVB)


### Compare results (for subgroup of PUMAS)
set.seed(1)
plotRed <- plotFinal %>% filter(PUMA %in% sample(unique(pcells$PUMA), 5))
truthRed <- truth %>% filter(PUMA %in% unique(plotRed$PUMA)) %>% mutate(IPRCAT=IPR)
ggplot(plotRed)+
  geom_violin(aes(x=Type, y=Est.)) +
  facet_grid(PUMA ~ IPRCAT, scales='free', labeller = 'label_both') +
  geom_hline(data=truthRed, color='red', linetype='dashed', aes(yintercept=P)) +
  theme_classic() 
