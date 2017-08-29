df_generate = function(n, p, seed) { 
  
  library(arm) # for the invlogit() function
  
  ## Generate covariate matrix (X)
  ## n = number of subjects
  ## p = number of baseline covariates
  X = matrix(rnorm(n*p, mean=0, sd=1), n, p) 
  
  ## Simulate a binary treatment indicator (Z)
  ## ps (propensity score): probability of being exposed conditonal on X
  ps = invlogit(-2 + 1.5*X[,1] + 1.75*X[,2] + 1.25*X[,3]) # logistic regression model
  # hist(ps_true)
  # seed = sample(1:500, 1)
  set.seed(seed)
  Z = rbinom(n, 1, ps)  # sample Bernoulli variable based on conditional probabilities
  # table(Z)
  
  ## Simulate a Gaussian outcome variable (Y)
  tau = 2
  Y = 10 + 1.75*X[,1] + 1.5*X[,2] + 1.25*X[,4] + tau*Z  # linear regression model with treatment effect
  # hist(Y)
  
  ## Bind data
  df = data.frame(X, Z, Y)
  return(list(df=df, ps=ps, tau=tau))

}