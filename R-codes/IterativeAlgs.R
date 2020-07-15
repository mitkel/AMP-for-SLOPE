# setting proper working dir and sourcing FastProx alg (works only in RStudio)
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("FastProxSL1.R")

# for working with S4 objects
library(methods)

# Generalized First Order Method Algorithm class
# iteration = number of the current iteration,
# loss = loss at the current state
# lossFun = loss function we want to minimize
setClass("GFOMA", 
         representation(iteration = "numeric",
                        loss = "numeric",
                        backtrack = "logical",
                        lossFun = "function"),
         prototype(iteration = 0,
                   loss = Inf,
                   backtrack = FALSE,
                   lossFun = function(lambda) 0)
         )
# (different for different methods)
setGeneric("algUpdate", function(x) standardGeneric("algUpdate"))
setGeneric("backTrack", function(x) standardGeneric("backTrack"))

# algorithm stopping condition
setGeneric("stopCond", function(x) standardGeneric("stopCond"))
# compares loss at previous state with the current loss and stops if the improvements is no bigger than some threshold 
# or # of iterations does not exceeds some other threshold
setMethod("stopCond", "GFOMA", function(x){
  isTRUE(abs(x@lossFun(x) - tail(as.vector(x@loss), n =1)) < 10**(-10) | x@iteration > 10**4)
})

setGeneric("runAlg", function(x) standardGeneric("runAlg"))
setMethod("runAlg", "GFOMA", function(x){
  while(!stopCond(x)){
    # save loss at the previous state (before update)
    x@loss <- c(x@loss, x@lossFun(x))
    
    if(x@backtrack) {
      x <- backTrack(x)}
    
    x <- algUpdate(x)
    
    x@iteration <- x@iteration + 1
  }
  return(x)
})

# auxilary function returning sorted l1 loss function ||Ax-y||_2^2/2 + SL1(x,theta)
SL1Loss <- function(foo) { sum( (foo@A %*% foo@x - foo@y)**2 )/2 + as.numeric(MapToDelta(foo@x)$w %*% MapToDelta(foo@theta)$w) }

# # ISTA
# y: observed values
# A: data matrix
# x: estimated vector (its starting value)
# theta: parameter of the penalty function
# stepSize: length of the step (chosen via backtracking)
setClass("ISTA-SLOPE", 
         representation(y = "numeric",
                        A = "matrix",
                        x = "numeric",
                        theta = "numeric",
                        stepSize = "numeric"),
         prototype(y = NA_real_,
                   A = as.matrix(NA_real_),
                   x = NA_real_,
                   theta = NA_real_,
                   stepSize = 1),
         contains = "GFOMA")

setMethod("algUpdate", "ISTA-SLOPE", function(x){
  foo <- x@x - x@stepSize * t(x@A) %*% (x@A %*% x@x - x@y)
  x@x <- as.numeric(FastProxSL1(foo, x@stepSize * x@theta))
  return(x)
})

# to do...
setMethod("backTrack", "ISTA-SLOPE", function(x){
  return(x)
})


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

setMethod("algUpdate", "FISTA-SLOPE", function(x){
  x@x_old <- x@x
  
  foo <- x@u - x@stepSize * t(x@A) %*% (x@A %*% x@u - x@y)
  x@x <- as.numeric(FastProxSL1(foo, x@stepSize * x@theta))
  
  x@t_old <- x@t
  x@t <- ((1+(1+4*(x@t_old)**(-2))**(1/2))/2)**(-1)
  
  x@u <- x@x + x@t*((x@t_old)**(-1)-1) * (x@x - x@x_old)
  
  return(x)
})

setClass("AMP-SLOPE",
         representation(x_old = "numeric",
                        z = "numeric",
                        z_old = "numeric",
                        theta_t = ),
         prototype(),
         contain = "ISTA")

setMethod("algUpdate", "AMP-SLOPE", function(x){
  0
})