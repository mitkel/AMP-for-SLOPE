library(gsubfn) # for multivalued functions handling

# MapToDelta:
#   Input: vector v
#   Output: vector w, vector sgn, permutation vector P
#     such that w in Delta_p, w = (sgn * v)[P]
MapToDelta <- function( v ){
  sgn <- ((v>0) - 1/2)*2
  foo <- sgn*v
  bar <- sort(foo, decreasing = TRUE, index.return = TRUE)
  return( list("w" = bar$x, "sgn" = sgn, "P" = bar$ix) )
}

# RetrV:
#   - is an inverse of MapToDelta()
#   Input: vector w, vector sgn, permutation vector P
#   Output: vector v 
#     such that w = (sgn * v)[P], i.e. v = w[InvP]*sgn, where InvP is the inverse permutation
RetrV <- function( w, sgn, P ){
  # retrieves the inverse permutation to P, i.e. 
  # finds [ j: P[j]=i ]_i
  foo <- as.numeric( 
    lapply(c(1:length(w)), function(x) (which(P == x)))
  )
  return( w[foo]*sgn )
}

# sanity check
foo <- rnorm(20)
bar <- MapToDelta(foo)
print( data.frame("v" = foo, "w" = bar$w) )
print( all.equal((foo*bar$sgn)[bar$P], bar$w) )
print( all.equal(RetrV(bar$w, bar$sgn, bar$P), foo) )

FastProxSL1 <- function( y, theta ){
  
}