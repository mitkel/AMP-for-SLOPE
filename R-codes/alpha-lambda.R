# # setting proper working dir and sourcing FastProx alg (works only in RStudio)
# library("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("FastProxSL1.R")

F <- function(tau, alpha, delta, prior, sigma2=0, iter = 100){
  p <- length(alpha)
  res <- 0
  for(i in 1:iter){
    foo <- prior()
    res <- (res*(i-1) + mean( (FastProxSL1(foo + tau*rnorm(p), alpha*tau) - foo)^2) )/i
  }
  
  return(sigma2+res/delta)
}

# psi <- function(tau, lambda, delta, prior, iter=100){
#   p <- length(lambda)
#   res <- 0
#   for(i in 1:iter){
#     foo <- prior
#     res <- (res*(i-1) + mean( (FastProxSL1(foo + tau/sqrt(delta)*rnorm(p), alpha*tau)-foo)**2))/i
#   }
#   return(res)
# }

ZeroAstNorm <- function(foo) length(unique(abs(foo[which(foo != 0)])))

alpha_to_tau_ast <- function(alpha, prior, delta, sigma2 = 0, iter1 = 50, iter2 = 20, 
                             tau = NA_real_, verbose=FALSE, ergodic=TRUE){
  
  if(is.na(tau)){
    foo <- 0
    for(i in 1:10) foo <- (foo*(i-1)+mean(prior()^2))/i
    tau <- sqrt(sigma2 + foo/delta)
  }
  res <- tau
  
  for(i in 1:iter1){
    tau <- sqrt(F(tau, alpha, delta, prior, sigma2, iter=iter2))
    res <- c(res,tau)
    
    if(verbose) print(paste0(i,": ",tau))
  }
  if(ergodic) res <- mean(res)
  else res <- tail(res,1)
  
  return(res)
}

alpha_to_lambda <- function(alpha, x, delta, sigma2 = 0, max_iter = 100, 
                            tau_ast = NA_real_, tau=NA_real_, prior = NULL, verbose = FALSE){
  if(is.na(tau_ast) && is.null(prior)) {
    warning("neither prior nor tau_ast specified")
    return(0)
  }
  
  if(is.na(tau_ast)){
    tau_ast <-  alpha_to_tau_ast(alpha, prior, delta, sigma2, tau=tau)
  }
  p <-  length(alpha)
  
  res <-  0
  for(i in 1:max_iter){
    res <-  (res*(i-1) + ZeroAstNorm(FastProxSL1(x + tau_ast*rnorm(p), tau_ast*alpha))/(delta*p) )/i
    if(verbose) print(paste0(i,": ",res))
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
