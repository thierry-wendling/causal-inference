## Estimation of the Average Treament Effect (tau) of a binary exposure variable (Z) on a Gaussian outcome (Y) using matching.

## In this script, we simulate observational data and evaluate the performance of different matching techniques.
## Observational data means that some subjects are more likely to be exposed to the treatment than others (i.e. exposure mechanism is not random like in randomised experiments)
## To estimate the ATE, we want to compare the outcome between the two treatment groups (control Z=0 and treated Z=1). However, we need to make sure that we are comparing similar subjects to avoid confounding. 

## Matching enables to extract a sample of subjects that have similar covariate distribution. This enables to elminate confounding and get an unbiased estimate of the ATE. 
## Matching can be done directly on observed covariates or on the propensity score depending on the knowledge of the true confounding factors and the dimension of the dataset.

rm(list=ls())
library(Matching)
setwd('/Users/mq20160354/Google Drive/R')
source('df_generate.R')

## Loop over data samples (assess impact of sampling noise)
Nsamp = 50
tau_relBias = matrix(0, Nsamp, 5)
for (k in 1:Nsamp) {

  ## ------------------------------------------------------------------------------------------------------
  ## Simualte data
  ## X: n x p matrix of pre-treatment covariates
  ## Z: Binary exposure indicator
  ## Y: Continuous outcome variable
  ## ------------------------------------------------------------------------------------------------------
  n = 500
  p = 9
  seed = k #= sample(1:50,1)
  ls = df_generate(n, p, seed)
  df = ls[['df']]
  ps_true = ls[['ps']]
  tau_true = ls[['tau']]
  rm(ls)
  
  ## ------------------------------------------------------------------------------------------------------
  ## Explore covariate balance and outcome difference between treatment groups 
  ## ------------------------------------------------------------------------------------------------------
  ## Outcome difference-in-means (Naive estimation of tau)
  tau_hat = mean(df$Y[df$Z==1]) - mean(df$Y[df$Z==0])
  tau_relBias_naive = abs(tau_hat - tau_true)/tau_true * 100
  
  ## Covariate difference-in-means
  lapply(1:p, function(j) {
    t.test(df[,j] ~ df$Z)
  })
  par(mfrow=c(3,3))
  for (j in 1:p) {
    hist(df[df$Z==0,j], main = paste(colnames(df)[j]), xlab = '')
    hist(df[df$Z==1,j], add = T, col = 'grey')
  }
  
  
  ## ------------------------------------------------------------------------------------------------------
  ## Estimate propensity score (probability of being exposed conditional on the baseline covariates)
  ## ------------------------------------------------------------------------------------------------------
  ## Logistic regression of Z on X
  fml = as.formula(paste('Z~', paste(colnames(df[,1:p]), sep = '', collapse = '+')))
  glm_fit = glm(fml, df, family = 'binomial')
  summary(glm_fit)
  ps_hat = glm_fit$fitted.values
  # plot(ps_hat, ps_true)
  # abline(0, 1, col='red')
  
  ## Examine common support
  table(df$Z)
  par(mfrow=c(1,1))
  hist(ps_hat[df$Z==0], main = 'Probability of being exposed', xlab = '')
  hist(ps_hat[df$Z==1], add = T, col = 'grey')
  legend('topright', c('Control','Exposed'), fill=c('white','grey'))
  
  ## ------------------------------------------------------------------------------------------------------
  ## Matching individuals based on observed confounders, or propensity score
  ## Propensity score: probability of being exposed to treatment conditional on baseline characteristics
  ## ------------------------------------------------------------------------------------------------------
  ## CASE 1: confounders are known (X1 and X2)
  match_out = Match(Y=df$Y, Tr=df$Z, X=df[,1:2], estimand = 'ATE')
  tau_hat = as.numeric(match_out$est)
  tau_CI_m1 = round(c(tau_hat-1.96*match_out$se, tau_hat+1.96*match_out$se), 2)
  tau_relBias_m1 = abs(tau_hat - tau_true)/tau_true * 100
  
  ## CASE 2: confounders are not known --> matching on all covariates
  match_out = Match(Y=df$Y, Tr=df$Z, X=df[,1:10], estimand = 'ATE')
  tau_hat = as.numeric(match_out$est)
  tau_CI_m2 = round(c(tau_hat-1.96*match_out$se, tau_hat+1.96*match_out$se), 2)
  tau_relBias_m2 = abs(tau_hat - tau_true)/tau_true * 100
  
  ## CASE 3: confounders are not known --> matching on propensity score (perfectly specified)
  match_out = Match(Y=df$Y, Tr=df$Z, X=ps_true, estimand = 'ATE') #, caliper = 0.2)
  tau_hat = as.numeric(match_out$est)
  tau_CI_m3 = round(c(tau_hat-1.96*match_out$se, tau_hat+1.96*match_out$se), 2)
  tau_relBias_m3 = abs(tau_hat - tau_true)/tau_true * 100
  
  ## CASE 4: confounders are not known --> matching on propensity score (perfectly specified)
  match_out = Match(Y=df$Y, Tr=df$Z, X=ps_hat, estimand = 'ATE') #, caliper = 0.2)
  tau_hat = as.numeric(match_out$est)
  tau_CI_m4 = round(c(tau_hat-1.96*match_out$se, tau_hat+1.96*match_out$se), 2)
  tau_relBias_m4 = abs(tau_hat - tau_true)/tau_true * 100
  
  ## Print results
  tau_relBias_k = round(c(tau_relBias_naive, tau_relBias_m1, tau_relBias_m2,
                          tau_relBias_m3, tau_relBias_m4), 3)
  out = data.frame(Method = c('Naive', 'Matching on known confounders', 'Matching on all variables',
                              'Matching on PS', 'Matching on estimated PS'),
                   Relative_Bias = tau_relBias_k)
  out = out[order(out$Relative_Bias),]
  knitr::kable(out)
  
  tau_relBias[k,] = tau_relBias_k


}

## Visualise results
par(mfrow=c(1,1))
boxplot(tau_relBias[,1], tau_relBias[,2], tau_relBias[,3], tau_relBias[,4], tau_relBias[,5],
        names = c('Naive', 'X1-2', 'Xall', 'PStrue', 'PShat'), ylab = 'Absolute Relative Bias in tau (%)',
        border = c('black','black','black','grey60','black'), main = 'Performance in estimating the ATE')
apply(tau_relBias, 2, mean)

