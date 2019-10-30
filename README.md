# Ensemble package in R
Functions for working with and combining ensemble models

# Example R code

```r
#install.packages("devtools")
#devtools::install_github("allanhicks/Ensemble")

library(Ensemble)
library(r4ss)
```

```r
model1 <- SS_output("Models/model1",verbose=F)
model2 <- SS_output("Models/model2",verbose=F)
model3 <- SS_output("Models/model3",verbose=F)
model4 <- SS_output("Models/model4",verbose=F)
models <- list(model1, model2, model3, model4)
```

```r
xx <- combineEnsemble.fn(models,param="SSB",element="derived_quants",yrs=1994:2019,totN=4000,useCov=TRUE)
apply(xx,2,median)
sum(xx[,"SPB_2017"]<xx[,"SPB_2016"])/nrow(xx)
sum(xx[,"SPB_2017"]<0.95*xx[,"SPB_2016"])/nrow(xx)
sum(xx[,"SPB_2019"]<xx[,"SPB_2016"])/nrow(xx)
sum(xx[,"SPB_2019"]<0.95*xx[,"SPB_2016"])/nrow(xx)

brat <- combineEnsemble.fn(models,param="Bratio",element="derived_quants",yrs=2015:2019,totN=4000,useCov=TRUE)
apply(brat,2,median)
sum(brat[,"Bratio_2017"]<0.3)/nrow(brat)

spr2016 <- combineEnsemble.fn(models,param="SPRratio",element="derived_quants",yrs=2016,totN=4000,useCov=TRUE)
quantile(1-spr2016,probs=c(0.025,0.5,0.975))

oflCatch <- combineEnsemble.fn(models,param="OFLCatch",element="derived_quants",yrs=2016:2018,totN=4000,useCov=TRUE)
quantile(oflCatch,probs=c(0.025,0.5,0.975))

#do this for multiple quantities at the same time
multQuants <- combineEnsemble.fn(models,param=c("SSB","Bratio"),element="derived_quants",yrs=2016:2018,totN=4000,useCov=TRUE)
apply(multQuants,2,quantile,probs=c(0.025,0.5,0.975))

#create a list separated by models for viewing individual results
indModels <- combineEnsemble.fn(models,param="SSB",element="derived_quants",yrs=2016:2018,totN=4000,useCov=TRUE,asList=T)
lapply(indModels,function(x){apply(x,2,quantile,probs=c(0.025,0.5,0.975))})

xx <- combineEnsemble.fn(models,param="NatM_p_1_Mal_GP_1",element="parameters",yrs=NULL,totN=4000,useCov=TRUE)
round(quantile(xx, probs=c(0.025, 0.5, 0.975)),4)

#note, when there are fixed parameters, the code will not work unless paired with estimated paramaters.
#you can also get creative with 'param' and 'yrs' arguments since they are combined
xx <- combineEnsemble.fn(models,param="NatM_p_1",element="parameters",yrs=c("Fem_GP_1","Mal_GP_1"),totN=4000,useCov=TRUE)
round(apply(xx,2,quantile,probs=c(0.025,0.5,0.975)),4)

```
