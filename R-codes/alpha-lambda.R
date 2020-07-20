# # setting proper working dir and sourcing FastProx alg (works only in RStudio)
# library("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("FastProxSL1.R")

F <- function(tau, alpha, delta, prior, sigma2=0, iter = 40){
  p <- length(alpha)
  res <- 0
  for(i in 1:iter){
    res <- (res*(i-1) + mean( (FastProxSL1(prior + tau*rnorm(p), alpha*tau)-prior)**2)/delta)/i
  }
  
  return(sigma2+res)
}

ZeroAstNorm <- function(foo) length(unique(abs(foo[which(foo != 0)])))

alpha_to_tau_ast <- function(alpha, prior, delta, sigma2 = 0, iter = 25){
  tau <- sqrt(sigma2 + mean(prior**2)/delta)
  res <- tau
  
  for(i in 1:iter){
    tau <- sqrt(F(tau, alpha, delta, prior, sigma2))
    res <- (res*(i-1) + tau)/i
  }
  
  return(res)
}

alpha_to_lambda <- function(alpha, prior, delta, sigma2 = 0, max_iter = 20){
  tau_ast <-  alpha_to_tau_ast(alpha, prior, delta, sigma2)
  p <-  length(alpha)
  
  res <-  0
  for(i in 1:max_iter){
    res <-  (res*(i-1) + ZeroAstNorm(FastProxSL1(prior + tau_ast*rnorm(p), tau_ast*alpha))/(delta*p) )/i
  }
  
  return((1-res)*alpha*tau_ast)
}

# by courtesy of https://github.com/woodyx218/SLOPE_AMP/blob/master/lambda_to_alpha.R
lambda_to_alpha <- function(lambda, prior, delta, sigma2 = 0, tol = 10**(-2)){
  l <-  lambda/lambda[1]
  
  alpha1 <- l/2
  alpha2 <- l*2
  
  lambda1 <- alpha_to_lambda(alpha1,prior=prior,delta=delta,sigma2=sigma2)
  lambda2 <- alpha_to_lambda(alpha2,prior=prior,delta=delta,sigma2=sigma2)
  
  while (((lambda1[1]<lambda[1]) && (lambda2[1]>lambda[1])) == FALSE){
    if (lambda1[1]<lambda[1]){
      alpha1 <- alpha1*2
      alpha2 <- alpha2*2
    } else{
      alpha1 <- alpha1/2
      alpha2 <- alpha2/2
    }
    lambda1 <- alpha_to_lambda(alpha1,prior=prior,delta=delta,sigma2=sigma2)
    lambda2 <- alpha_to_lambda(alpha2,prior=prior,delta=delta,sigma2=sigma2)
  }
  
  ### bisection to find the alpha_seq which is parallel to lambda_seq
  while ((alpha2[1]-alpha1[1])>tol){
    middle_alpha <- (alpha1+alpha2)/2
    middle_lambda <- alpha_to_lambda(middle_alpha, prior=prior, delta=delta, sigma2=sigma2)
    if (middle_lambda[1]>lambda[1]){
      alpha2 <- middle_alpha
    }else if (middle_lambda[1]<lambda[1]){
      alpha1 <- middle_alpha
    }else{
      break
    }
  }
  
  return(middle_alpha)
}
