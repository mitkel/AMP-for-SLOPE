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

# DeltaCheck:
#   - since MapToDelta() requires a sorting step, we want a faster method of checking if v in Delta_p
#   Input: vector v
#   Output: TRUE iff v in Delta_p
DeltaCheck <- function( v ){
  # return( isTRUE(all.equal(v, MapToDelta(v)$w )) )
  return( max(diff(v)) <= 0 && min(v) >= 0  )
}

# FastProxSL1:
#   Input: vect. v, vect. theta
#   Output: argmin_x { |x-y|^2/2 + phi_theta(x) } = (y-theta)_+ if (y-theta) in Delta_p
FastProxSL1 <- function( y, theta){
  if( !DeltaCheck(theta) ) warning("Arg. theta is not in Delta_p!") 
  foo <- MapToDelta(y)
  y <- foo$w
  bar <- as.numeric( lapply( y-theta, function(x) max(x,0) ) )
  
  while( !DeltaCheck(bar) ){
   i <- which(diff(y-theta) > 0)[1] # begining of the first nondecreasing subsequence in y-theta
   j <- i + which(diff( (y-theta)[-(1:i)] ) < 0)[1] # end of the first nondecreasing subsequence in y-theta
   if(is.na(j)) {j <- length(y)} # to handle the case when the tail is increasing
   y[i:j] <- mean(y[i:j])
   theta[i:j] <- mean(theta[i:j])
   bar <- as.numeric( lapply( y-theta, function(x) max(x,0) ) )
  }
  return( RetrV(bar, foo$sgn, foo$P))
}

# ===============================================================
# # sanity check
# foo <- rnorm(5)
# bar <- MapToDelta(foo)
# print( data.frame("v" = foo, "w" = bar$w) )
# print( isTRUE( all.equal((foo*bar$sgn)[bar$P], bar$w) ) )
# print( isTRUE( all.equal(RetrV(bar$w, bar$sgn, bar$P), foo) ) )
# 
# # for theta = const, ProxSL1 is a shrinkage operator
# print(foo <- sample(seq(from=-4, to=4, by=1), 9))
# print(FastProxSL1(foo, theta=rep(2, times=9)))

