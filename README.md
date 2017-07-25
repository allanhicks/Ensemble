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
xx <- combineEnsemble.fn(models,param="SPB",element="derived_quants",yrs=1994:2019,totN=4e6,useCov=TRUE)
apply(xx,2,median)
sum(xx[,"SPB_2017"]<xx[,"SPB_2016"])/nrow(xx)
sum(xx[,"SPB_2017"]<0.95*xx[,"SPB_2016"])/nrow(xx)
sum(xx[,"SPB_2019"]<xx[,"SPB_2016"])/nrow(xx)
sum(xx[,"SPB_2019"]<0.95*xx[,"SPB_2016"])/nrow(xx)

brat <- combineEnsemble.fn(models,param="Bratio",element="derived_quants",yrs=2015:2019,totN=4e6,useCov=TRUE)
apply(brat,2,median)
sum(brat[,"Bratio_2017"]<0.3)/nrow(brat)

spr2016 <- combineEnsemble.fn(models,param="SPRratio",element="derived_quants",yrs=2016,totN=4e6,useCov=TRUE)
quantile(1-spr2016,probs=c(0.025,0.5,0.975))

foreCatch <- combineEnsemble.fn(models,param="ForeCatch",element="derived_quants",yrs=2016:2018,totN=4e6,useCov=TRUE)
quantile(foreCatch,probs=c(0.025,0.5,0.975))
```
