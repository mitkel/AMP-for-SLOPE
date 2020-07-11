# setting proper working dir and sourcing FastProx alg (works only in RStudio)
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("FastProxSL1.R")

# Generalized First Order Method Algorithm
setClass("GFOMA", representation(params = "list"))
setClass("ISTA", representation(), contains = "GFOMA")
setClass("FISTA", representation(), contains = "GFOMA")