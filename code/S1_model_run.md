---
title: 'Yaguarundí IDM'
author: 'Florencia Grattarola, Diana Bowler & Petr Keil'
date: '2022-09-08'
output:
  html_document:
    keep_md: yes
    toc: yes
    toc_float: yes
    highlight: pygments
    theme: flatly
---

Analysis for the *Herpailurus yagouaroundi*, with new covariates: `npp`, `elev`, `bio7`, and `bio15`.

  - Libraries


```r
library(sf)
library(R2jags) # interface to JAGS
library(tidyverse)
```

# Data


```r
# Presence-absence data
PA_time1 <- readRDS('data/data_hyagouaroundi_PA_time1.rds') %>%
  filter(!is.na(env.npp) &!is.na(env.bio_15) & !is.na(env.bio_7) & !is.na(env.elev)) # remove NA's
PA_time2 <- readRDS('data/data_hyagouaroundi_PA_time2.rds') %>%
  filter(!is.na(env.npp) &!is.na(env.bio_15) & !is.na(env.bio_7) & !is.na(env.elev)) # remove NA's


# Presence-only data
PO_time1 <- readRDS('data/data_hyagouaroundi_PO_time1.rds') %>%
  filter(!is.na(env.npp) &!is.na(env.bio_15) & !is.na(env.bio_7) & !is.na(env.elev) & !is.na(acce) & !is.na(count)) # remove NA's
PO_time2  <- readRDS('data/data_hyagouaroundi_PO_time2.rds')  %>%
  filter(!is.na(env.npp) &!is.na(env.bio_15) & !is.na(env.bio_7) & !is.na(env.elev) & !is.na(acce) & !is.na(count)) # remove NA's


PA_time1_time2 <- rbind(PA_time1 %>% mutate(time=1), PA_time2 %>% mutate(time=2)) 
PO_time1_time2 <- rbind(PO_time1 %>% mutate(time=1), PO_time2 %>% mutate(time=2))
```


## Splines

Before going fully Bayesian in JAGS, we will do some GAMs in `mgcv` to (1) extract spline basis functions (using `mgcv::jagam`) that will then be used in JAGS, and (2) get some reasonable starting values for coefficients for the MCMC algorithm in JAGS.

  -   Set the **`k` parameter**: the number of basis functions.  
  -   Set the **model formula**.  
  -   Fit the **GAM model**.  


```r
k = 10

gam.formula <- formula(count ~ env.npp + 
                               env.bio_15 + 
                               env.bio_7 + 
                               env.elev +
                               offset(log(area)) +
                               s(X, Y, by=as.factor(time), k=k))

model.gam <- mgcv::gam(gam.formula, 
                       data = PO_time1_time2,
                       family = poisson)

summary(model.gam)
```

```
## 
## Family: poisson 
## Link function: log 
## 
## Formula:
## count ~ env.npp + env.bio_15 + env.bio_7 + env.elev + offset(log(area)) + 
##     s(X, Y, by = as.factor(time), k = k)
## 
## Parametric coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -29.4024     0.3161 -93.030  < 2e-16 ***
## env.npp      -0.2125     0.1028  -2.067   0.0387 *  
## env.bio_15    1.1392     0.2715   4.196 2.72e-05 ***
## env.bio_7    -1.1839     0.1288  -9.188  < 2e-16 ***
## env.elev      1.0939     0.2588   4.227 2.37e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                           edf Ref.df Chi.sq p-value    
## s(X,Y):as.factor(time)1 8.974      9  212.6  <2e-16 ***
## s(X,Y):as.factor(time)2 8.987      9  555.3  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0424   Deviance explained = 38.7%
## UBRE = -0.47236  Scale est. = 1         n = 4366
```

  -   Get the **jagam splines** values: using `mgcv::jagam`.  
  
To get spline values for both types of data we first run the code for one type of data (presence-only in this case, given the data are more widely distributed), and then use the model to predict on the other type of data.


```r
jagam.out <- mgcv::jagam(gam.formula, 
                         data = PO_time1_time2,
                         family = poisson,
                         file='') # we will not need this file
```

```
## model {
##   eta <- X %*% b + offset ## linear predictor
##   for (i in 1:n) { mu[i] <-  exp(eta[i]) } ## expected response
##   for (i in 1:n) { y[i] ~ dpois(mu[i]) } ## response 
##   ## Parametric effect priors CHECK tau=1/8.8^2 is appropriate!
##   for (i in 1:5) { b[i] ~ dnorm(0,0.013) }
##   ## prior for s(X,Y):as.factor(time)1... 
##   K1 <- S1[1:9,1:9] * lambda[1]  + S1[1:9,10:18] * lambda[2]
##   b[6:14] ~ dmnorm(zero[6:14],K1) 
##   ## prior for s(X,Y):as.factor(time)2... 
##   K2 <- S2[1:9,1:9] * lambda[3]  + S2[1:9,10:18] * lambda[4]
##   b[15:23] ~ dmnorm(zero[15:23],K2) 
##   ## smoothing parameter priors CHECK...
##   for (i in 1:4) {
##     lambda[i] ~ dgamma(.05,.005)
##     rho[i] <- log(lambda[i])
##   }
## }
```

```r
smooth1 <- jagam.out$pregam$smooth[[1]]
smooth2 <- jagam.out$pregam$smooth[[2]]

# presence-only data
jagam.PO.time1 <- mgcv::PredictMat(object = smooth1, data = PO_time1_time2)
jagam.PO.time2 <- mgcv::PredictMat(object = smooth2, data = PO_time1_time2)

# presence-absence data
jagam.PA.time1 <- mgcv::PredictMat(object = smooth1, data = PA_time1_time2)
jagam.PA.time2 <- mgcv::PredictMat(object = smooth2, data = PA_time1_time2)
```

# Prepare data for JAGS


```r
PA.X <- cbind('(Intercept)'= 1, 
              'env.npp' = PA_time1_time2$env.npp,
              'env.bio_15' = PA_time1_time2$env.bio_15,
              'env.bio_7' = PA_time1_time2$env.bio_7,
              'env.elev' = PA_time1_time2$env.elev,
              jagam.PA.time1, 
              jagam.PA.time2)

PO.X <- cbind('(Intercept)'= 1, 
              'env.npp' = PO_time1_time2$env.npp,
              'env.bio_15' = PO_time1_time2$env.bio_15,
              'env.bio_7' = PO_time1_time2$env.bio_7,
              'env.elev' = PO_time1_time2$env.elev,
              jagam.PO.time1, 
              jagam.PO.time2)

# number of all columns in X
n.X = ncol(PA.X)

# number of columns in X of spline basis functions
n.spl = ncol(jagam.PA.time1)

# number of columns in X of env. predictors + intercept
n.par = n.X - ncol(jagam.PA.time1)*2

# number of factors of time in X 
n.fac = length(unique(as.factor(PO_time1_time2$time)))*2

# country as an indexed variable
PO.country  <- as.numeric(as.factor(PO_time1_time2$country))

#number of countries
n.country <- length(unique(PO.country))

jags.data <- list(n.PA = nrow(PA_time1_time2),
                  y.PA = PA_time1_time2$presabs, 
                  X.PA = PA.X,
                  area.PA = PA_time1_time2$area,
                  effort = PA_time1_time2$effort,
                  n.PO = nrow(PO_time1_time2),
                  n.PO.half = nrow(PO_time1_time2)/2,
                  y.PO = PO_time1_time2$count, 
                  X.PO = PO.X,
                  area.PO = PO_time1_time2$area, 
                  acce = PO_time1_time2$acce,
                  country = PO.country,
                  n.X = n.X,
                  n.cntr = n.country,
                  n.par = n.par,
                  n.fac = n.fac,
                  n.spl = n.spl,
                  Z = rep(0, length(jagam.out$jags.data$zero)),
                  S.time1 = jagam.out$jags.data$S1,
                  S.time2 = jagam.out$jags.data$S2) 
```

# Specify the model in BUGS language


```r
cat('model
  { 
    # PRIORS --------------------------------------------------
    
    ## Thinning at locations with complete accessibility in PO data
      
      # intercept of the decay function for each country of origin. 
      # It needs a flat prior between 0 and 1
      for (c in 1:n.cntr) 
      {
        alpha0[c] ~ dbeta(1, 1)  
      }

      # steepness of the decaying distance-P.ret relationship in PO data
      alpha1 ~ dgamma(0.5, 0.05)   
      
    ## Effect of sampling effort in PA data
      beta ~ dnorm(0, 0.01)
  
    ## Parametric effects of environment driving the point process intensity
     # (it also includes an intercept)
     
      for (r in 1:n.par)
      { 
        b[r] ~ dnorm(0,0.01) 
      }
      
    ## Splines (imported and adjusted form output of mgcv::jagam)
    
      ## prior for s(X,Y):as.factor(time)1 
      sigma.time1 <- S.time1[1:n.spl, 1:n.spl] * gamma[1]  + 
                     S.time1[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[2]
      b[(n.par+1):(n.spl + n.par)] ~ dmnorm(Z[(n.par+1):(n.spl + n.par)], sigma.time1) 
     
      ## prior for s(X,Y):as.factor(time)2
      sigma.time2 <- S.time2[1:n.spl, 1:n.spl] * gamma[3]  + 
                     S.time2[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[4]
      b[(n.X - n.spl + 1):(n.X)] ~ dmnorm(Z[(n.X - n.spl + 1):(n.X)], sigma.time2) 
     
      ## Priors for smoothing parameter 
      for (f in 1:n.fac) 
      {
        gamma[f] ~ dgamma(.5,.5)
        rho[f] <- log(gamma[f])
      }

    # LIKELIHOOD --------------------------------------------------
    
      ## --- Presence-Absence (PA) data ---
      
       eta.PA <- X.PA %*% b ## linear predictor
        
       for (i in 1:n.PA) 
       { 
         # the probability of presence
         cloglog(psi[i]) <- eta.PA[i] + log(area.PA[i]) + beta*log(effort[i]) 
        
         # presences and absences come from a Bernoulli distribution
         y.PA[i] ~ dbern(psi[i]*0.9999) 
         
       } 
  
      ## --- Presence-Only (PO) data --- 
  
      eta.PO <- X.PO %*% b  ## linear predictor

      for (j in 1:n.PO)
      {
        # cell-specific probability of retainin (observing) a point is a function of accessibility
        P.ret[j] <- alpha0[country[j]] * exp( (-alpha1) * acce[j]) 
        
        # true mean number (nu) of points per cell i is the true intensity multiplied by cell area
        nu[j] <- area.PO[j] * exp(eta.PO[j]) 

        # thinning: the true lambda
        lambda[j] <- nu[j] * P.ret[j]

        # counts of observed points come from a Poisson distribution
        y.PO[j] ~ dpois(lambda[j]) 
      }
  
    # PREDICTIONS -------------------------------------------------

    eta.pred <- X.PO %*% b

    for (j in 1:n.PO)
    {
      # predicted probability of occurrence in grid cell j
      cloglog(P.pred[j]) <- eta.pred[j] + log(area.PO[j])
    }

    # POSTERIOR PREDICTIVE CHECK  --------------------------------
    
    # for PA
    for (i in 1:n.PA)
    {
      # Fit assessments: Tjur R-Squared (fit statistic for logistic regression)
      pres[i] <- ifelse(y.PA[i] > 0, psi[i], 0)
      absc[i] <- ifelse(y.PA[i] == 0, psi[i], 0)
    }
    
    # Discrepancy measures for entire PA data set
    pres.n <- sum(y.PA[] > 0)
    absc.n <- sum(y.PA[] == 0)
    r2_tjur <- abs(sum(pres[])/pres.n - sum(absc[])/absc.n)

    # for PO
    for (j in 1:n.PO)
    {
      # Fit assessments: Posterior predictive check and data for DHARMA
      y.PO.new[j] ~ dpois(lambda[j]) 

    }

    # DERIVED QUANTITIES ------------------------------------------
    
    # area in each time period, and temporal change of area
    A.time1 <- sum(P.pred[1:n.PO.half])
    A.time2 <- sum(P.pred[(n.PO.half+1):n.PO])
    delta.A <- A.time2 - A.time1
    
    # uncertainty for the temporal change
    for (j in 1:n.PO.half)
    {
      delta.Grid[j] <- P.pred[n.PO.half+j] - P.pred[j]
    }
  }
', file = 'model/yaguarundi_model.txt')
```


## Inits function


```r
jags.inits <- function(model = model.gam)
{
  return(list(b=rnorm(n.X, mean=model$coefficients, sd=1)))
}
```

# Fit the model


```r
start.time <- Sys.time()

yaguarundi_model <- R2jags::jags(data=jags.data,
                                 model.file='model/yaguarundi_model.txt',
                                 parameters.to.save=c('b', 'P.pred', 
                                                      'A.time1', 'A.time2', 'delta.A',
                                                      'alpha0', 'alpha1', 'beta',
                                                      'lambda', 'P.ret', 'psi',
                                                      'y.PO.new', 'r2_tjur', 
                                                      'delta.Grid'),
                                 inits = jags.inits,
                                 n.chains=3,
                                 n.iter=100000,
                                 n.thin=10,
                                 n.burnin=10000,
                                 DIC = FALSE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# saveRDS(yaguarundi_model, 'big_data/yaguarundi_model_fit.rds')
```
