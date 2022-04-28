---
title: "Quick Reference of Model Definitions"
author: 'Florencia Grattarola, Diana Bowler & Petr Keil'
date: '30/04/2022'
output: md_document
---

## Quick Reference of Model Definitions

See more details [here](/model_definitions.md)


<table class=" lightable-paper lightable-striped lightable-hover" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> Model term </th>
   <th style="text-align:left;"> Definition </th>
   <th style="text-align:left;"> Equation notation </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> n.PA </td>
   <td style="text-align:left;"> number of blobs for both time periods (pre and pos) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{PA}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{PA}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> i </td>
   <td style="text-align:left;"> index identifying blobs </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large i \text{where} i \in 1:n_{PA}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large i \text{where} i \in 1:n_{PA}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> y.PA[i] </td>
   <td style="text-align:left;"> presence (1) or absence (0) value in each i-th blob (overlapping surveys' area), can be for pre- or post- period </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large y_{PA_i}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large y_{PA_i}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> X.PA </td>
   <td style="text-align:left;"> design matrix including vector of 1s (for intercept) and all the covariates and spline bases for each blob, for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{X_{PA}}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{X_{PA}}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> area.PA[i] </td>
   <td style="text-align:left;"> area of i-th blob in meters for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large area_{PA_i}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large area_{PA_i}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> effort[i] </td>
   <td style="text-align:left;"> sampling effort for i-th blob in the given period for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large effort_{PA_i}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large effort_{PA_i}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.PO </td>
   <td style="text-align:left;"> number of grid-cells for both time periods concatenated (pre and pos) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{PO}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{PO}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> j </td>
   <td style="text-align:left;"> index identifying grid cells </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large j \text{where} j \in j:n_{PO}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large j \text{where} j \in j:n_{PO}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.PO.half </td>
   <td style="text-align:left;"> number of grid-cells for one time period </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{PO/2}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{PO/2}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> y.PO[j] </td>
   <td style="text-align:left;"> count of observed points in j-th grid-cell, can be for pre- or post- period </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large y_{PO_j}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large y_{PO_j}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> X.PO </td>
   <td style="text-align:left;"> design matrix including vector of 1s (for intercept) and all the covariates and spline bases for each grid-cell for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{X_{PO}}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{X_{PO}}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> area.PO[j] </td>
   <td style="text-align:left;"> area of each grid-cell in meters for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large area_{PO_j}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large area_{PO_j}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> acce[j] </td>
   <td style="text-align:left;"> accessibility from urban areas based on travel time for j-th grid-cell for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large acce_j}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large acce_j}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> country[j] </td>
   <td style="text-align:left;"> country name for j-th grid-cell for both time periods </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large country_j}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large country_j}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.X </td>
   <td style="text-align:left;"> total number of columns in X (`X.PA` or `X.PO`) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_X}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_X}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.cntr </td>
   <td style="text-align:left;"> total number of countries </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{cntr}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{cntr}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> c </td>
   <td style="text-align:left;"> index identifying countries </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large c \text{where} c \in 1:n_{cntr}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large c \text{where} c \in 1:n_{cntr}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.par </td>
   <td style="text-align:left;"> number of parameters considered (intercept and covariates) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{par}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{par}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> r </td>
   <td style="text-align:left;"> index identifying parameters </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large r \text{where} r \in 1:n_{par}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large r \text{where} r \in 1:n_{par}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.fac </td>
   <td style="text-align:left;"> number of factors of time in X (`X.PA` or `X.PO`) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{fac}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{fac}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> f </td>
   <td style="text-align:left;"> index identifying factors </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large f \text{where} f \in 1:n_{fac}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large f \text{where} f \in 1:n_{fac}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> n.spl </td>
   <td style="text-align:left;"> number of spline bases functions in in X (`X.PA` or `X.PO`) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large n_{spl}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large n_{spl}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> S.pre </td>
   <td style="text-align:left;"> spline values for the first time period (pre) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large S_{pre}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large S_{pre}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> S.post </td>
   <td style="text-align:left;"> spline values for the second time period (post) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large S_{post}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large S_{post}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> Z </td>
   <td style="text-align:left;"> a vector of zeros (0) of the length of the splines </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large Z}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large Z}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> sigma.pre </td>
   <td style="text-align:left;"> variance of splines for the first time period (pre) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \sigma_{pre}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \sigma_{pre}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> sigma.post </td>
   <td style="text-align:left;"> variance of splines for the second time period (post) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \sigma_{post}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \sigma_{post}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> b </td>
   <td style="text-align:left;"> vector of parametric effects of covariates driving the point process intensity (it also includes an intercept) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large b_r \in \mathbf{b}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large b_r \in \mathbf{b}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> alpha0 </td>
   <td style="text-align:left;"> intercept of the thinning process in the presence-only data </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \alpha_0}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \alpha_0}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> alhpa1 </td>
   <td style="text-align:left;"> slope -steepness- of the thinning process in the presence-only data (decaying distance~P.ret relationship) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \alpha_1}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \alpha_1}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> beta </td>
   <td style="text-align:left;"> coefficient of the effect of sampling effort in the presence-absence data </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \beta}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \beta}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma </td>
   <td style="text-align:left;"> prior for splines smoothing parameter </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \gamma}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \gamma}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> eta.PA </td>
   <td style="text-align:left;"> linear predictor for presence-absence data </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{\eta_{PA}}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{\eta_{PA}}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> eta.PA[i] </td>
   <td style="text-align:left;"> expected presence-absence for the i-th blob </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \eta_{PA_i}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \eta_{PA_i}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> eta.PO </td>
   <td style="text-align:left;"> linear predictor for presence-only data </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{\eta_{PO}}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{\eta_{PO}}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> eta.PO[j] </td>
   <td style="text-align:left;"> expected count points for the j-th grid-cell </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \eta_{PO_j}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \eta_{PO_j}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> psi[i] </td>
   <td style="text-align:left;"> blob-specific probability of presence </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large psi_i}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large psi_i}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> P.ret[j] </td>
   <td style="text-align:left;"> cell-specific probability of retaining (observing) a point as a function of accessibility and country of origin </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large P_{ret_j}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large P_{ret_j}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> nu[j] </td>
   <td style="text-align:left;"> true mean number of points per grid-cell (the true intensity) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \nu_j}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \nu_j}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> lambda[j] </td>
   <td style="text-align:left;"> thinning of the true intensity </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \lambda_j}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \lambda_j}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> eta.pred </td>
   <td style="text-align:left;"> linear predictor for the predicted probability of occurrence </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{\eta_{pred}}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \mathbf{\eta_{pred}}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> eta.pred[j] </td>
   <td style="text-align:left;"> predicted count points for the j-th grid-cell </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \eta_{pred_j}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \eta_{pred_j}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> P.pred[j] </td>
   <td style="text-align:left;"> predicted probability of occurrence for the j-th grid-cell </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large P_{pred_j}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large P_{pred_j}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> A.pre </td>
   <td style="text-align:left;"> range area in the first time period (pre) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large A_{pre}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large A_{pre}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> A.post </td>
   <td style="text-align:left;"> range area in the second time period (post) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large A_{post}}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large A_{post}}#gh-dark-mode-only"></td>
  </tr>
  <tr>
   <td style="text-align:left;"> delta.A </td>
   <td style="text-align:left;"> temporal change of range area (post-pre) </td>
   <td style="text-align:left;"> <img src="https://render.githubusercontent.com/render/math?math={\Large \Delta A}#gh-light-mode-only">
                                 <img src="https://render.githubusercontent.com/render/math?math={\Large \Delta A}#gh-dark-mode-only"></td>
  </tr>
</tbody>
</table>

