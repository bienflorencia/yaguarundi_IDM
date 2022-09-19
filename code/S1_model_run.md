---
title: 'Yaguarundí IDM'
author: 'Florencia Grattarola, Diana Bowler & Petr Keil'
date: '30/04/2022'
output: 
  html_document: 
    keep_md: yes
    toc: yes
    highlight: pygments
    theme: flatly
    number_sections: yes
---



# Libraries


```r
library(R2jags) # interface to JAGS
library(tidyverse)
```

# Data

Loading the pre- and post-2014 data for presence-absence (PA) and presence-only (PO).

**Covariates**

 - **env.elev**: Elevation (*)  
 - **env.npp**: Net Primary Productivity (*)  
 - **env.bio_7**: Temperature Annual Range (*)  
 - **env.bio_15**: Precipitation Seasonality (Coefficient of Variation) (*)  


```r
# Presence-absence data
PA_pre <- readRDS('data/PA_pre.rds') %>% 
  filter(!is.na(env.elev) &!is.na(env.bio_4) & !is.na(env.npp)) # remove NA's
PA_post <- readRDS('data/PA_post.rds') %>% 
  filter(!is.na(env.elev) &!is.na(env.bio_4) & !is.na(env.npp)) # remove NA's

# Presence-only data
PO_pre <- readRDS('data/PO_pre.rds') %>% 
  filter(!is.na(env.elev) &!is.na(env.bio_4) & !is.na(env.npp) & !is.na(acce) & !is.na(count)) # remove NA's
PO_post  <- readRDS('data/PO_post.rds') %>% 
  filter(!is.na(env.elev) &!is.na(env.bio_4) & !is.na(env.npp) & !is.na(acce) & !is.na(count)) # remove NA's

PA_pre_post <- rbind(PA_pre %>% mutate(time=0), PA_post %>% mutate(time=1)) 
PO_pre_post <- rbind(PO_pre %>% mutate(time=0), PO_post %>% mutate(time=1)) 
```

# Generalized Additive models from `mgcv`

Before going fully Bayesian in JAGS, we will do some GAMs in `mgcv` to (1) extract spline basis functions (using `mgcv::jagam`) that will then be used in JAGS, and (2) get some reasonable starting values for coefficients for the MCMC algorithm in JAGS.

## Setting the scene for GAM

The **`k` parameter** setting the number of basis functions:


```r
k = 9
```

**Model formula**. we set it here so that we can then use it at different places
in the code, if needed.


```r
gam.formula <- formula(count ~ env.elev + 
                               env.bio_7 + 
                               env.bio_15 + 
                               env.npp +
                               offset(log(area)) +
                               s(X, Y, by=as.factor(time), k=k))
```

## Fitting the GAM model
                               

```r
model.gam <- mgcv::gam(gam.formula, 
                       data = PO_pre_post,
                       family = poisson)

summary(model.gam)
```

```
## 
## Family: poisson 
## Link function: log 
## 
## Formula:
## count ~ env.elev + env.bio_7 + env.bio_15 + env.npp + offset(log(area)) + 
##     s(X, Y, by = as.factor(time), k = k)
## 
## Parametric coefficients:
##              Estimate Std. Error  z value Pr(>|z|)    
## (Intercept) -28.02812    0.20154 -139.068  < 2e-16 ***
## env.elev      1.71599    0.23118    7.423 1.15e-13 ***
## env.bio_7    -0.93753    0.11493   -8.157 3.42e-16 ***
## env.bio_15    2.32194    0.25441    9.127  < 2e-16 ***
## env.npp       0.80207    0.08286    9.680  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                           edf Ref.df Chi.sq p-value    
## s(X,Y):as.factor(time)0 7.892  7.996  290.4  <2e-16 ***
## s(X,Y):as.factor(time)1 7.932  7.999  604.6  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0417   Deviance explained = 35.1%
## UBRE = -0.44329  Scale est. = 1         n = 4366
```


# Jagam splines

To get spline values for both types of data we use the `mgcv::jagam`. We first run the code for one type of data (presence-only in this case, given the data are more widely distributed), and then use the model to predict on the other type of data.


```r
jagam.out <- mgcv::jagam(gam.formula, 
                         data = PO_pre_post,
                         family = poisson,
                         file='') # we will not need this file
```

```
## model {
##   eta <- X %*% b + offset ## linear predictor
##   for (i in 1:n) { mu[i] <-  exp(eta[i]) } ## expected response
##   for (i in 1:n) { y[i] ~ dpois(mu[i]) } ## response 
##   ## Parametric effect priors CHECK tau=1/8.5^2 is appropriate!
##   for (i in 1:5) { b[i] ~ dnorm(0,0.014) }
##   ## prior for s(X,Y):as.factor(time)0... 
##   K1 <- S1[1:8,1:8] * lambda[1]  + S1[1:8,9:16] * lambda[2]
##   b[6:13] ~ dmnorm(zero[6:13],K1) 
##   ## prior for s(X,Y):as.factor(time)1... 
##   K2 <- S2[1:8,1:8] * lambda[3]  + S2[1:8,9:16] * lambda[4]
##   b[14:21] ~ dmnorm(zero[14:21],K2) 
##   ## smoothing parameter priors CHECK...
##   for (i in 1:4) {
##     lambda[i] ~ dgamma(.05,.005)
##     rho[i] <- log(lambda[i])
##   }
## }
```

```r
smooth0 <- jagam.out$pregam$smooth[[1]]
smooth1 <- jagam.out$pregam$smooth[[2]]

# presence-only data
jagam.PO.time0 <- mgcv::PredictMat(object = smooth0, data = PO_pre_post)
jagam.PO.time1 <- mgcv::PredictMat(object = smooth1, data = PO_pre_post)

# presence-absence data
jagam.PA.time0 <- mgcv::PredictMat(object = smooth0, data = PA_pre_post)
jagam.PA.time1 <- mgcv::PredictMat(object = smooth1, data = PA_pre_post)
```

# The model in JAGS

## Prepare data for JAGS


```r
PA.X <- cbind('(Intercept)'= 1, 
              'env.elev' = PA_pre_post$env.elev,
              'env.bio_7' = PA_pre_post$env.bio_7,
              'env.bio_15' = PA_pre_post$env.bio_15,
              'env.npp' = PA_pre_post$env.npp,
              jagam.PA.time0, 
              jagam.PA.time1)

PO.X <- cbind('(Intercept)'= 1, 
              'env.elev' = PO_pre_post$env.elev,
              'env.bio_7' = PO_pre_post$env.bio_7,
              'env.bio_15' = PO_pre_post$env.bio_15,
              'env.npp' = PO_pre_post$env.npp,
              jagam.PO.time0, 
              jagam.PO.time1)

# number of all columns in X
n.X = ncol(PA.X)

# number of columns in X of spline basis functions
n.spl = ncol(jagam.PA.time0)

# number of columns in X of env. predictors + intercept
n.par = n.X - ncol(jagam.PA.time0)*2

# number of factors of time in X 
n.fac = length(unique(as.factor(PO_pre_post$time)))*2

# country as an indexed variable
PO.country  <- as.numeric(as.factor(PO_pre_post$country))

#number of countries
n.country <- length(unique(PO.country))

jags.data <- list(n.PA = nrow(PA_pre_post),
                  y.PA = PA_pre_post$presabs, 
                  X.PA = PA.X,
                  area.PA = PA_pre_post$area,
                  effort = PA_pre_post$effort,
                  n.PO = nrow(PO_pre_post),
                  n.PO.half = nrow(PO_pre_post)/2,
                  y.PO = PO_pre_post$count, 
                  X.PO = PO.X,
                  area.PO = PO_pre_post$area, 
                  acce = PO_pre_post$acce,
                  country = PO.country,
                  n.X = n.X,
                  n.cntr = n.country,
                  n.par = n.par,
                  n.fac = n.fac,
                  n.spl = n.spl,
                  Z = rep(0, length(jagam.out$jags.data$zero)),
                  S.pre = jagam.out$jags.data$S1,
                  S.post = jagam.out$jags.data$S2) 
```


## Specify the model in BUGS language


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
    
      ## prior for s(X,Y):as.factor(time)0 
      sigma.pre <- S.pre[1:n.spl, 1:n.spl] * gamma[1]  + 
            S.pre[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[2]
      b[(n.par+1):(n.spl + n.par)] ~ dmnorm(Z[(n.par+1):(n.spl + n.par)], sigma.pre) 
     
      ## prior for s(X,Y):as.factor(time)1
      sigma.post <- S.post[1:n.spl, 1:n.spl] * gamma[3]  + 
            S.post[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[4]
      b[(n.X - n.spl + 1):(n.X)] ~ dmnorm(Z[(n.X - n.spl + 1):(n.X)], sigma.post) 
     
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
      y.PA.new[i] ~ dbern(psi[i]*0.9999)   # replicate (new) data set
      
      pres[i] <- ifelse(y.PA[i] > 0, y.PA.new[i], 0)
      absc[i] <- ifelse(y.PA[i] == 0, y.PA.new[i], 0)
    }
    
    # Discrepancy measures for entire PA data set
    pres.n <- sum(y.PA.new[] > 0)
    absc.n <- sum(y.PA.new[] == 0)
    r2_tjur <- abs(sum(pres[])/pres.n - sum(absc[])/absc.n)

    # for PO
    for (j in 1:n.PO)
    {
      # Fit assessments: Freeman-Tukey test and posterior predictive check
      ppft[j] <- (sqrt(y.PO[j]) - sqrt(lambda[j]))^2          # observed
      y.PO.new[j] ~ dpois(lambda[j])                          # replicate (new) data set
      ppft.new[j] <- (sqrt(y.PO.new[j]) - sqrt(lambda[j]))^2  # expected
    }
    
    # Discrepancy measures for entire data set
    mean.abs.lambda.diff <- mean(abs(lambda[] - y.PO[]))
    
    # Add up discrepancy measures for entire data set
    fit.PO <- sum(ppft[])                     
    fit.PO.new <- sum(ppft.new[])             
    
    # DERIVED QUANTITIES ------------------------------------------
    
    # area in each time period, and temporal change of area
    A.pre <- sum(P.pred[1:n.PO.half])
    A.post <- sum(P.pred[(n.PO.half+1):n.PO])
    delta.A <- A.post - A.pre
    
    # uncertainty for the temporal change
    for (j in 1:n.PO.half)
    {
      delta.Grid[j] <- P.pred[n.PO.half+j] - P.pred[j]
    }
  }
', file = 'model/yaguarundi_model.txt')
```

## Inits function

We are taking the parameters estimated in the `mgcv` GAM object, plus we add a bit
of noise - just to give the algorithm some reasonable values to start with.


```r
jags.inits <- function(model = model.gam)
{
  return(list(b=rnorm(n.X, mean=model$coefficients, sd=1)))
}
```


## Fit the model


```r
start.time <- Sys.time()

yaguarundi_model <- R2jags::jags(data=jags.data,
                                          model.file='model/yaguarundi_model.txt',
                                          parameters.to.save=c('b', 'P.pred', 
                                                      'A.pre', 'A.post', 'delta.A',
                                                      'alpha0', 'alpha1', 'beta',
                                                      'lambda', 'P.ret', 'psi',
                                                      'y.PO.new', 'y.PA.new',
                                                      'r2_tjur', 
                                                      'fit.PO', 'fit.PO.new',
                                                      'mean.abs.lambda.diff',
                                                      'delta.Grid'),
                                          inits = jags.inits,
                                          n.chains=3,
                                          n.iter=100000,
                                          n.thin=10,
                                          n.burnin=10000,
                                          DIC = FALSE)
```

```
## module glm loaded
```

```
## module dic loaded
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 5318
##    Unobserved stochastic nodes: 5358
##    Total graph size: 210845
## 
## Initializing model
## 
##   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
##   |**************************************************| 100%
```

```r
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
```

```
## Time difference of 3.425978 hours
```

```r
# save the fitted model to a file
# saveRDS(yaguarundi_model, 'data/yaguarundi_model_fit.rds')
```

