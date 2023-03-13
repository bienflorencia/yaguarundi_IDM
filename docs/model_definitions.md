---
title: "Model Definitions and Equations"
author: 'Florencia Grattarola, Diana Bowler & Petr Keil'
date: "2023-03-13"
output: 
  html_document: 
    keep_md: yes
    toc: no
    highlight: pygments
    theme: cerulean
    number_sections: no
---



# DAG

![DAG](docs/figs/DAG.png)

# Definitions (variables and parameters)

**Presence-absence data (camera-trap surveys)**  

  - `n.PA`: number of blobs for both time periods (time1 and time2) -- **notation:** $n_{PA}$   
  - `i`: index identifying blobs, $i$, where $i \in 1:n_{PA}$ -- **notation:** $i$ 
  - `y.PA[i]`: presence (1) or absence (0) value in each $i$-th blob (overlapping surveys' area), can be for time1 or time2 -- **notation:** $y_{PA_i}$   
  - `X.PA`: design matrix including vector of 1s (for intercept) and all the covariates and spline bases for each blob, for both time periods -- **notation:** $\mathbf{X_{PA}}$  
  - `area.PA[i]`: area of $i$-th blob in meters for both time periods -- **notation:** $area_{PA_i}$  
  - `effort[i]`: sampling effort for $i$-th blob in the given period for both time periods -- **notation:** $effort_i$   
  
**Presence-only data (occurrence records)**  

  - `n.PO`: number of grid-cells for both time periods (time1 and time2) -- **notation:** $n_{PO}$    
  - `j`: index identifying grid cells -- **notation** $j$, where $j \in 1:{n_{PO}}$
  - `n.PO.half`: number of grid-cells for one time period -- **notation** $n_{PO/2}$  
  - `y.PO[j]`: count of observed points in $j$-th grid-cell, can be for time1 or time2 -- **notation:** $y_{PO_j}$
  - `X.PO`: design matrix including vector of 1s (for intercept) and all the covariates and spline bases for each grid-cell for both time periods -- **notation:** $X_{PO}$  
  - `area.PO[j]`: area of each grid-cell in meters for both time periods -- **notation:** $area_{PO_j}$  
  - `acce[j]`: accessibility from urban areas based on travel time for $j$-th grid-cell for both time periods -- **notation:** $acce_j$  
  - `country[j]`: country name for $j$-th grid-cell for both time periods -- **notation:** $country_j$  

**Auxiliary variables** 

  - `n.X`: total number of columns in X  (`X.PA` or `X.PO`) -- **notation:** $n_X$  
  - `n.cntr`: total number of countries -- **notation:** $n_{cntr}$  
  - `c`: index identifying countries, where $c \in 1:n_{cntr}$
  - `n.par`: number of parameters considered (intercept and covariates) -- **notation:** $n_{par}$
  - `r`: index identifying parameters, where $r \in 1:n_{par}$
  - `n.fac`: number of factors of time in X  (`X.PA` or `X.PO`) -- **notation:** $n_{fac}$ 
  - `f`: index identifying factors, where $f \in 1:n_{fac}$  

**Splines**  

  - `k`: number of spline basis functions (**not** in the model, only used in `mgcv::jagam` analysis) -- **notation:** `k`  
  - `n.spl`: number of spline basis functions in in X  (`X.PA` or `X.PO`) -- **notation:** $n_{spl}$   
  - `S.time1`: spline values for the first time period (time1) -- **notation:** $S_{time1}$  
  - `S.time2`: spline values for the second time period (time2) -- **notation:** $S_{time2}$  
  - `Z`: a vector of zeros (0) of the length of the splines -- **notation:** $Z$  
  - `sigma.time1`: variance of splines for the first time period (time1) -- **notation:** $\sigma_{time1}$  
  - `sigma.time2`: variance of splines for the second time period (time2) -- **notation:** $\sigma_{time2}$ 

**Model hyperparameters**  

  - `b`: vector of parametric effects of covariates driving the point process intensity (it also includes an intercept) -- **notation:** $b_r \in \mathbf{b}$ and $r \in 1:n_{par}$ 
  - `alpha0`: intercept of the thinning process -- **notation:** $\alpha_0$    
  - `alhpa1`: slope -steepness- of the thinning process (decaying distance~P.ret relationship) in presence-only data -- **notation:** $\alpha_1$  
  - `beta`: coefficient of the effect of sampling effort in the presence-absence data -- **notation:** $\beta$  
  - `gamma`: prior for smoothing parameter -- **notation:** $\gamma$ 

**Calculated variables**

  - `eta.PA`: linear predictor for presence-absence data -- **notation:** $\mathbf{\eta_{PA}}$      
  - `eta.PA[i]`: expected presence-absence for the $i$-th blob -- **notation:** $\eta_{PA_i}$  
  - `eta.PO`: linear predictor for presence-only data -- **notation:** $\mathbf{\eta_{PO}}$    
  - `eta.PO[j]`: expected count points for the $j$-th grid-cell -- **notation:** $\eta_{PO_j}$
  - `psi[i]`: blob-specific probability of presence -- **notation:** $\psi_i$ 
  - `P.ret[j]`: cell-specific probability of retaining (observing) a point as a function of accessibility and country of origin -- **notation:** $P_{ret_j}$   
  - `nu[j]`: true mean number of points per grid-cell (the true intensity) -- **notation:** $\nu_j$    
  - `lambda[j]`: thinning of the true intensity -- **notation:** $\lambda_j$

**Predicted parameters**  

  - `eta.pred`: linear predictor for the predicted probability of occurrence -- **notation:** $\mathbf{\eta_{pred}}$   
  - `eta.pred[j]`: predicted count points for the j-th grid-cell -- **notation:** $\mathbf{\eta_{PO}}$ 
  - `P.pred[j]`: predicted probability of occurrence for the j-th grid-cell -- **notation:** $P_{pred_j}$  

**Derived quantities**  

  - `A.time1`, range area in the first time period (time1) -- **notation:** $A_{time1}$
  - `A.time2`: range area in the second time period (time2) -- **notation:** $A_{time2}$  
  - `delta.A`: temporal change of range area (time2-time1) -- **notation:** $\Delta A$ 
  - `delta.Grid`: uncertainty (SD) of the temporal change (time2-time1) -- **notation:** $\Delta SD$  
  
  
# Equations

## Priors

 - Thinning at locations with complete accessibility in PO data

  Intercept of the decay function varying for country of origin (it needs a flat prior between 0 and 1)  
  
  **JAGS code**:  `for (i in 1:n.cntr) { alpha0[c] ~ dbeta(1, 1) }`  
  **Equation**:  $\alpha_{0_c} \sim \textsf{Beta}(1,1)$ where $c \in 1:n_{cntr}$  

  Steepness of the decaying distance-P.ret relationship in PO data   
  
  **JAGS code**:  `alpha1 ~ dgamma(0.5, 0.05)`  
  **Equation** :  $\alpha_1 \sim \textsf{Gamma}(0.5, 0.05)$  

 - Effect of sampling effort in presence-absence (PA) data  

  **JAGS code**:  `beta ~ dnorm(0, 0.01)`  
  **Equation**:  $\beta \sim \textsf{Normal}(0, 0.01)$  
      
 - Parametric effects of environment driving the point process intensity (it also includes an intercept)

  **JAGS code**:  `for (r in 1:n.par) { b[r] ~ dnorm(0,0.01) }`  
  **Equation**:  $b_r \sim \textsf{Normal}(0, 0.01)$ where $r \in 1:n_{par}$   

 - Splines (imported and adjusted form output of `mgcv::jagam`)  
  
  prior for s(X,Y):as.factor(time)0  
  **JAGS code:** `sigma.time1 <- S.time1[1:n.spl, 1:n.spl] * gamma[1] + S.time1[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[2]`  
  **Equation:** $\sigma_{time1} = S_{time1_{1:n.spl, 1:n.spl}} \times \gamma_1  + S_{time1_{1:n.spl, n.spl+1:n.spl \times 2}} \times \gamma_2$  
  
  **JAGS code:** `b[(n.par+1):(n.spl + n.par)] ~ dmnorm(Z[(n.par+1):(n.spl + n.par)], sigma.time1)`  
  **Equation:** $b_{n.par+1:n.spl+n.par} \sim \textsf{Normal}(Z_{n.par+1:n.spl+n.par}, \sigma_{time1})$  

  prior for s(X,Y):as.factor(time)1  
  **JAGS code:** `sigma.time2 <- S.time2[1:n.spl, 1:n.spl] * gamma[3] + S.time2[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[4]`  
  **Equation:** $\sigma_{time2} = S_{time2_{1:n.spl, 1:n.spl}} \times \gamma_3  + S_{time2_{1:n.spl, n.spl+1:n.spl \times 2}} \times \gamma_4$  
  
  **JAGS code:** `b[(n.X - n.spl + 1):(n.X)] ~ dmnorm(Z[(n.X - n.spl + 1):(n.X)], sigma.time2)`  
  **Equation:** $b_{n.X-n.spl+1:n.X} \sim \textsf{Normal}(Z_{n.X-n.spl+1:n.X}, \sigma_{time2})$  
      
 - Priors for smoothing parameter  
 
  **JAGS code:** `for (f in 1:n.fac) { gamma[f] ~ dgamma(.5,.5)}`   
  **Equation:** $\gamma_{f} \sim \textsf{Gamma}(0.5, 0.5)$   

  **JAGS code:** `for (f in 1:n.fac) { rho[f] <- log(gamma[f])}`   
  **Equation:** $\rho_{f} = \log(gamma_f)$   


## Likelihood

 - Presence-Absence (PA) data 
  
  linear predictor       

  **JAGS code**: `eta.PA <- X.PA %*% b`  
  **Equation**: $\mathbf{\eta_{PA}} = \mathbf{X_{PA}} \times \mathbf{b}$  

  **JAGS code**: `for (i in 1:n.PA) { cloglog(psi[i]) <- eta.PA[i] + log(area.PA[i]) + beta*log(effort[i]) }`  
  **Equation**: $cloglog(\psi_i) = \eta_{PA_i} + \log(area_{PA_i}) + \beta \times log(effort_i)$  

  **JAGS code**: `for (i in 1:n.PA) { y.PA[i] ~ dbern(psi[i]*0.9999) }`  
  **Equation**: $y_{PA_i} \sim Bernoulli(\psi_i \times 0.9999)$  

 - Presence-Only (PO) data
  
  linear predictor  
  **JAGS code:** `eta.PO <- X.PO %*% b`  
  **Equation:** $\mathbf{\eta_{PO}} = \mathbf{X_{PO}} \times \mathbf{b}$   

  Cell-specific probability of retaining (observing) a point is a function of accessibility  
  **JAGS code:** `for (j in 1:n.PO) { P.ret[j] <- alpha0[country[j]] * exp( (-alpha1) * acce[j]) }`  
  **Equation:** $P_{ret_j} = \alpha_{0c} \times \exp^{-\alpha_1 \times acce_j}$  
        
  true mean number (lambda) of points per cell i is the true intensity multiplied by cell area  
  **JAGS code:** `for (j in 1:n.PO) { nu[j] <- area.PO[j] * exp(eta.PO[j]) }`  
  **Equation:** $\nu_j = area_{PO_j}  \times \exp^{\eta_{PO_{j}}}$  

  thinning the true lambda  
  **JAGS code:** `for (j in 1:n.PO) { lambda[j] <- nu[j] * P.ret[j] }`  
  **Equation:** $\lambda_j = \nu_j \times P_{ret_j}$  

  counts of observed points come from Poisson distribution  
  **JAGS code:** `for (j in 1:n.PO) { y.PO[j] ~ dpois(lambda[j]) }`  
  **Equation:** $y_{PO_j} \sim Poisson(\lambda_j)$  
  

## Predictions
  
  linear predictor  
  **JAGS code:** `eta.pred <- X.PO %*% b`  
  **Equation:** $\mathbf{\eta_{pred}} =  \mathbf{X_{PO}} \times \mathbf{b}$  

  Predicted probability of occurrence in grid cell i  
  **JAGS code:** `for (j in 1:n.PO) { cloglog(P.pred[j]) <- eta.pred[j] + log(area.PO[j]) }`  
  **Equation:** $cloglog(P_{pred_j}) =  \eta_{pred_j} + \log(area_{PO_j})$  

 - Derived Quantities 
    
  range size in each time period, and temporal change of range size  
  **JAGS code:** `A.time1 <- sum(P.pred[1:n.PO.half])`  
  **Equation:** $A_{time1} = \sum \eta_{pred_j}$ where $j \in 1:n_{PO}/2$  
  
  **JAGS code:** `A.time2 <- sum(P.pred[(n.PO.half+1):n.PO])`  
  **Equation:** $A_{time2} = \sum \eta_{pred_j}$ where $j \in n_{PO}/2:n_{PO}$  
  
  **JAGS code:** `delta.A <- A.time2 - A.time1`  
  **Equation:** $\Delta A  = A_{time2} - A_{time1}$  

  **JAGS code:** `delta.Grid[j] <- P.pred[n.PO.half+j] - P.pred[j]`  
  **Equation:** $\Delta SD_j  = P_{pred_{n_{PO}/2+j}} - P_{pred_j}$  

## Posterior predictive checks

 - Presence-Absence (PA) data 
 
  fit assessments: Tjur R-Squared (fit statistic for logistic regression)
  **JAGS code:** `y.PA.new[i] ~ dbern(psi[i]*0.9999)`   # replicate (new) data set  
  **Equation:**
  **JAGS code:** `pres[i] <- ifelse(y.PA[i] > 0, y.PA.new[i], 0)`  
  **Equation:**
  **JAGS code:** `absc[i] <- ifelse(y.PA[i] == 0, y.PA.new[i], 0)`  
  **Equation:**
  
  discrepancy measures for entire PA data set
  **JAGS code:** `pres.n <- sum(y.PA.new[] > 0)`  
  **Equation:**
  **JAGS code:** `absc.n <- sum(y.PA.new[] == 0)`   
  **Equation:**
  **JAGS code:** `r2_tjur <- abs(sum(pres[])/pres.n - sum(absc[])/absc.n)`  
  **Equation:**

 - Presence-Only (PO) data
 
  discrepancy measures for entire data set
  **JAGS code:** `mean.abs.lambda.diff <- mean(abs(lambda[] - y.PO[]))`
  **JAGS code:** `fit.PO <- sum(ppft[])`                    
  **JAGS code:** `fit.PO.new <- sum(ppft.new[])`
  
# Updated model in BUGS language


```r
 'model
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
      sigma.time1 <- S.time1[1:n.spl, 1:n.spl] * gamma[1]  + 
                   S.time1[1:n.spl, (n.spl + 1):(n.spl * 2)] * gamma[2]
            b[(n.par+1):(n.spl + n.par)] ~ dmnorm(Z[(n.par+1):(n.spl + n.par)], sigma.time1) 
     
      ## prior for s(X,Y):as.factor(time)1
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
'
```
