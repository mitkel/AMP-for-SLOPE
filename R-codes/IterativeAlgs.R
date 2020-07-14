# setting proper working dir and sourcing FastProx alg (works only in RStudio)
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("FastProxSL1.R")

# for working with S4 objects
library(methods)

# Generalized First Order Method Algorithm class
# params = vector of parameters,
# loss = loss function we want to minimize
setClass("GFOMA", 
         representation(params = "list",
                        loss = "function"),
         prototype(params = list(),
                   loss = function(lambda) 0)
         )
# algorithm update rule (different for different methods)
setGeneric("algUpdate", function(x) standardGeneric("algUpdate"))

# algorithm stopping condition
setGeneric("stopCond", function(x) standardGeneric("stopCond"))
# compares loss at previous state with the current loss and stops if the improvements is no bigger than some threshold 
# or # of iterations does not exceeds some other threshold
setMethod("stopCond", "GFOMA", function(x){
  isTRUE(abs(x@loss(x@params) - x@params$loss) < 10**(-9) | x@params$iteration > 10**4)
})

setGeneric("runAlg", function(x) standardGeneric("runAlg"))
setMethod("runAlg", "GFOMA", function(x){
  while(!stopCond(x)){
    # save loss at the previous state (before update)
    x@params$loss <- x@loss(x@params)
    
    # state update
    x@params <- algUpdate(x)
    
    # increment iteration
    x@params$iteration <- x@params$iteration + 1
  }
  return(x@params)
})

# # ISTA
# params - required to have fields:
#   - y: observed values
#   - A: data matrix 
#   - x: estimated vector (its starting value)
#   - theta: parameter of the penalty function
#   - stepSize: length of the step (chosen via backtracking)
#   - iteration: number of the current iteration of the algorithm
#   - loss: loss function at the previous state
setClass("ISTA-SLOPE", contains = "GFOMA")
setMethod("algUpdate", "ISTA-SLOPE", function(x){
  
  # here should be the backtracking procedure
  # while( x@loss( list("x" = bar, "theta" = x@params$theta*x@params$stepSize) ) > 0 ){
  #   x@params$stepSize <- x@params$stepSize * 0.9
  # }
  
  # update rule for ISTA
  foo <- x@params$x - x@params$stepSize * t(x@params$A) %*% (x@params$A %*% x@params$x - x@params$y)
  x@params$x <- FastProxSL1(foo, x@params$stepSize * x@params$theta)
  return(x@params)
})


# # FISTA
# params - required to have fields:
#   - y: observed values
#   - A: data matrix
#   - x: estimated vector (its starting value)
#   - x_old: previous state (x_old(0)=x(0))
#   - theta: parameter of the penalty function
#   - t: acceleration parameter
#   - t_old: as above (t_old(0)=t(0))
#   - stepSize: length of the step (chosen via backtracking)
#   - iteration: number of the current iteration of the algorithm
#   - loss: loss function at the previous state
setClass("FISTA-SLOPE", representation(theta = "numeric"), contains = "GFOMA")
setMethod("algUpdate", "FISTA-SLOPE", function(x){
  # update rule for FISTA
  x@params$x_old <- x@params$x
  
  foo <- x@params$u - x@params$stepSize * t(x@params$A) %*% (x@params$A %*% x@params$u - x@params$y)
  x@params$x <- FastProxSL1(foo, x@params$stepSize * x@params$theta)
  
  x@params$t_old <- x@params$t
  x@params$t <- ((1+(1+4*(x@params$t_old)**(-2))**(1/2))/2)**(-1)
  
  x@params$u <- x@params$x + x@params$t*((x@params$t_old)**(-1)-1) * (x@params$x - x@params$x_old)
  
  return(x@params)
})

# sanity check
ISTA <- new("ISTA-SLOPE",
            params = list("y" = c(6,6),
                          "A" = diag(c(3,2)),
                          "x" = c(23,-10),
                          "theta" = c(1,1),
                          "stepSize" = 0.01,
                          "iteration" = 0,
                          "loss" = Inf),
            loss = SL1Loss)
FISTA <- new("FISTA-SLOPE",
            params = list("y" = c(6,6),
                          "A" = diag(c(3,2)),
                          "x" = c(23,-10),
                          "x_old" = c(23,-10),
                          "u" = c(23,-10),
                          "theta" = c(1,1),
                          "stepSize" = 0.01,
                          "t" = 1,
                          "t_old" = 1,
                          "iteration" = 0,
                          "loss" = Inf),
            loss = SL1Loss)
print(ISTA_res <- runAlg(ISTA))
print(FISTA_res <- runAlg(FISTA))

