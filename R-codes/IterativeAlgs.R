# # setting proper working dir and sourcing FastProx alg (works only in RStudio)
# library("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("FastProxSL1.R")

# auxilary function returning sorted l1 loss function ||Ax-y||_2^2/2 + SL1(x,lambda)
# SL1Loss <- function(foo) {
#   sum( (foo@A %*% foo@x - foo@y)**2 )/2 + 
#     as.numeric(MapToDelta(foo@x)$w %*% MapToDelta(foo@lambda)$w) 
# }
SL1Loss <- function(A,x,y,lambda){
  sum( (A %*% x - y)^2)/2 + as.numeric(MapToDelta(x)$w %*% MapToDelta(lambda)$w)
}

# for working with S4 objects
library(methods)

# Generalized First Order Method Algorithm class
# iteration = number of the current iteration,
# loss = loss at the current state
# lossFun = loss function we want to minimize
# tol and max_iter - stopping condition parameters
setClass("GFOMA", 
         representation(iteration = "numeric",
                        loss = "numeric",
                        backtrack = "logical",
                        lossFun = "function",
                        tol = "numeric",
                        max_iter = "numeric",
                        init_params = "logical",
                        verbose = "logical"),
         prototype(iteration = 0,
                   loss = Inf,
                   backtrack = TRUE,
                   lossFun = SL1Loss,
                   tol = 10**(-4),
                   max_iter = 10**3,
                   init_params = TRUE,
                   verbose = FALSE)
         )

# (different for different methods)
setGeneric("algUpdate", function(x) standardGeneric("algUpdate"))

# backTracking procedure -- by default does not do anything
setGeneric("backTrack", function(x) standardGeneric("backTrack"))
setMethod("backTrack", "GFOMA", function(x) x)

# initParams -- needed for initializing intrinsic parameters of the AMP, by default does not change x
setGeneric("initParams", function(x) standardGeneric("initParams"))
setMethod("initParams", "GFOMA", function(x) x)

# algorithm stopping condition
setGeneric("stopCond", function(x) standardGeneric("stopCond"))
# compares loss at previous state with the current loss and stops if the improvements is no bigger than some threshold 
# or # of iterations does not exceeds some other threshold
setMethod("stopCond", "GFOMA", function(x){
  # isTRUE(abs(x@lossFun(x) - tail(as.vector(x@loss), n =1)) < x@tol | x@iteration > x@max_iter)
  i <- x@iteration
  if(i < 2){
    return(FALSE)
  } else{
    return(isTRUE(abs(x@loss[i]-x@loss[i-1])/x@loss[i-1] < x@tol || x@iteration > x@max_iter))
  }
})

setGeneric("runAlg", function(x) standardGeneric("runAlg"))
setMethod("runAlg", "GFOMA", function(x){
  if(x@verbose) print("begin algorithm")
  
  # update some defaulr params (needed for AMP tau^ast, alpha etc.)
  if(x@init_params){
    if(x@verbose) print("begin initParams")
    x <- initParams(x)
    if(x@verbose) print("initialized params")
  }
  
  
  x <- algUpdate(x)

  return(x)
})

# -------------------------------------------------------------------
# # ISTA
# y: observed values
# A: data matrix
# x: estimated vector (its starting value)
# lambda: parameter of the penalty function
# stepSize: length of the step (chosen via backtracking)
setClass("ISTA-SLOPE", 
         representation(y = "numeric",
                        A = "matrix",
                        x = "numeric",
                        lambda = "numeric",
                        stepSize = "numeric"),
         prototype(y = NA_real_,
                   A = as.matrix(NA_real_),
                   x = NA_real_,
                   lambda = NA_real_,
                   stepSize = 0.1),
         contains = "GFOMA")

setMethod("backTrack", "ISTA-SLOPE", function(x){
  x_old <- x@x
  stepSize <- x@stepSize
  repeat{
    foo <- x_old - stepSize * t(x@A) %*% (x@A %*% x_old - x@y)
    x@x <- as.numeric(FastProxSL1(foo, stepSize * x@lambda))
    Qxy <- sum( (x@A%*%x_old-x@y)**2 )/2 +
      (x@x-x_old) %*% t(x@A) %*% (x@A %*% x_old - x@y) +
      0.5*sum((x@x-x_old)**2)/stepSize +
      MapToDelta(x@lambda)$w %*% MapToDelta(x@x)$w
    
    if(x@lossFun(x@A, x@x, x@y, x@lambda) <= Qxy || x@stepSize < 10^(-6)){
      break
    } else {
      stepSize <- stepSize * 0.9
      if(x@verbose) print(stepSize)
    }
    x@x <- x_old
    x@stepSize <- stepSize
  }
  return(x)
})

setMethod("algUpdate", "ISTA-SLOPE", function(x){
  while(!stopCond(x)){
    x@loss <- c(x@loss, x@lossFun(x@A, x@x, x@y, x@lambda))
    
    if(x@backtrack) {
      x <- backTrack(x)}
    
    foo <- x@x - x@stepSize * t(x@A) %*% (x@A %*% x@x - x@y)
    x@x <- as.numeric(FastProxSL1(foo, x@stepSize * x@lambda))
    
    x@iteration <- x@iteration + 1
    if(x@iteration %% 5 == 1 && x@verbose) print(paste0("ISTA-iter: ",x@iteration))
  }  
  return(x)
})

# -------------------------------------------------------
# # FISTA
# x_old: previous state (x_old(0)=x(0))
# t: acceleration parameter
# t_old: as above (t_old(0)=t(0))
setClass("FISTA-SLOPE", 
         representation(x_old = "numeric",
                        u = "numeric",
                        t = "numeric",
                        t_old = "numeric"),
         prototype(x_old = NA_real_,
                   u = NA_real_,
                   t = 1,
                   t_old = 1),
         contains = "ISTA-SLOPE")

setMethod("initParams", "FISTA-SLOPE", function(x){
  x@x_old <- x@x
  x@u <- x@x
  return(x)
})

setMethod("backTrack", "FISTA-SLOPE", function(x){
  u_old <- x@u
  stepSize <- x@stepSize
  repeat{
    foo <- u_old - stepSize * t(x@A) %*% (x@A %*% u_old - x@y)
    u_new <- as.numeric(FastProxSL1(foo, stepSize * x@lambda))
    Qxy <- sum( (x@A%*%u_old-x@y)**2 )/2 +
      (u_new-u_old) %*% t(x@A) %*% (x@A %*% u_old - x@y) +
      0.5*sum((u_new-u_old)**2)/stepSize +
      MapToDelta(x@lambda)$w %*% MapToDelta(u_new)$w
    
    Fxy <- sum( (x@A %*% u_new - x@y)**2 )/2 +
      as.numeric(MapToDelta(u_new)$w %*% MapToDelta(x@lambda)$w) 
    
    if(Fxy <= Qxy || x@stepSize < 10^(-6)){
      break
    } else{
      stepSize <- stepSize * 0.9
      if(x@verbose) print(x@stepSize)
    }
    
    x@stepSize <- stepSize
  }
  return(x)
})

setMethod("algUpdate", "FISTA-SLOPE", function(x){
  while(!stopCond(x)){
    x@loss <- c(x@loss, x@lossFun(x@A, x@x, x@y, x@lambda))
    
    if(x@backtrack) {
      x <- backTrack(x)}
    
    x@x_old <- x@x
    
    foo <- x@u - x@stepSize * t(x@A) %*% (x@A %*% x@u - x@y)
    x@x <- as.numeric(FastProxSL1(foo, x@stepSize * x@lambda))
    
    x@t_old <- x@t
    x@t <- ((1+(1+4*(x@t_old)**(-2))**(1/2))/2)**(-1)
    
    x@u <- x@x + x@t*((x@t_old)**(-1)-1) * (x@x - x@x_old)
    
    x@iteration <- x@iteration + 1
    if(x@iteration %% 5 == 1 && x@verbose) print(paste0("FISTA-iter: ",x@iteration))
  }
  return(x)
})

#----------------------------------------------------------------
# AMP-SLOPE
# source("alpha-lambda.R")

setClass("AMP-SLOPE",
         representation(y = "numeric",
                        A = "matrix",
                        x = "numeric",
                        x_old = "numeric",
                        v = "numeric",
                        alpha = "numeric",
                        lambda = "numeric",
                        prior = "function",
                        n = "numeric",
                        p = "numeric",
                        delta = "numeric",
                        tau = "numeric",
                        sigma2 = "numeric",
                        F_iter = "numeric"),
         prototype(y = NA_real_,
                   A = as.matrix(NA_real_),
                   x = NA_real_,
                   x_old = NA_real_,
                   v = NA_real_,
                   lambda = NA_real_,
                   n = NA_real_,
                   p = NA_real_,
                   delta = NA_real_,
                   tau = NA_real_,
                   alpha = NA_real_,
                   sigma2 = 0,
                   F_iter = 10^3),
         contain = "GFOMA")

setMethod("initParams", "AMP-SLOPE", function(x) {
  if(x@verbose) print("initParams::begin")
  
  if(anyNA(x@n)) x@n <- dim(x@A)[1]
  if(anyNA(x@p)) x@p <- dim(x@A)[2]
  if(anyNA(x@delta)) x@delta <- x@n/x@p
  if(anyNA(x@x)) x@x = rep(0,x@p)
  if(anyNA(x@x_old)) x@x_old = rep(0,x@p)
  if(anyNA(x@v)) x@v <- as.numeric( t(x@A) %*% (x@A %*% x@x - x@y))
  
  if(is.na(x@tau)){
    foo <- 0
    for(i in 1:10) foo <- (foo*(i-1)+mean(x@prior()**2))/i
    x@tau <- sqrt(x@sigma2 + foo/delta)
  }
  
  if(anyNA(x@alpha)){
    warning("alpha not specified!")
    return(NULL)
  } 
  
  if(anyNA(x@lambda)){
    x@lambda <- alpha_to_lambda(x@alpha, x@prior, x@delta, x@tau)
  }
  
  if(x@verbose) print("initParams::end")
  return(x)
})

setMethod("algUpdate", "AMP-SLOPE", function(x){
  while(!stopCond(x)){
    if(x@iteration %% 5 == 0 && x@verbose) print(paste0("Begin AMP-iter: ",x@iteration))
    
    x@loss <- c(x@loss, x@lossFun(x@A, x@x, x@y, x@lambda))
    
    if(x@backtrack) {
      x <- backTrack(x)}
    
    x@x_old <- x@x
    x@x <- as.numeric(FastProxSL1(x@x_old - x@v, x@alpha * x@tau))
    
    x@v <- as.numeric( t(x@A) %*% (x@A %*% x@x - x@y) + x@v*ZeroAstNorm(x@x)/x@n )
    
    tau_old <- x@tau
    x@tau <- sqrt( F(x@tau, x@alpha, x@delta, x@prior, x@sigma2, iter= x@F_iter) )
    # x@t <- x@t*x@tau/tau_old
    
    if(x@iteration %% 5 == 0 && x@verbose) print(paste0("End AMP-iter: ", x@iteration))
    x@iteration <- x@iteration + 1
  }
  
  return(x)
})
