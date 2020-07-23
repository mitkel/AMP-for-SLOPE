source("R-codes/IterativeAlgs.R")
source("R-codes/alpha-lambda.R")
source("R-codes/FastProxSL1.R")

library(ggplot2, reshape2)

set.seed(351759)

p=500;
delta=0.5;
eps=0.1;
sigma2=0.2;
n=p*delta
A=matrix(rnorm(n*p,mean=0,sd=1/sqrt(p*delta)), n,p) 
alpha=c(rep(1,p*delta),rep(0,p*(1-delta)))
# alpha <- rep(1,p)
stepSize=100;
maxIter=10**3

prior= function() rnorm(p)*rbinom(p,1,eps)

x0=prior()
x=prior() 
y=as.numeric(A%*%x+sqrt(sigma2)*rnorm(n))
tau <- sqrt(sigma2+eps/delta)
tau_ast <- alpha_to_tau_ast(alpha, x, delta, tau=tau, sigma2 = sigma2, verbose= TRUE, ergodic = TRUE)
lambdahat <- alpha_to_lambda(alpha, x, delta, tau_ast = tau_ast, sigma2 = sigma2)


print(paste0("True loss:",loss <- sum(A %*% x - y)^2/2 + MapToDelta(x)$w %*% MapToDelta(lambdahat)$w))

AMP <-   new("AMP-SLOPE", 
             y=y, A=A, prior=prior, tau=tau, lambda=lambdahat,
             alpha=alpha, sigma2=sigma2, backtrack=FALSE, verbose=FALSE)
AMP_res <- runAlg(AMP)
print(paste0("AMP loss=", paste0(tail(AMP_res@loss,1), collapse = "; "), " after ",AMP_res@iteration," iterations"))

ISTA <- new("ISTA-SLOPE", y=y, A=A, x=x0, lambda=lambdahat, stepSize=stepSize, 
            max_iter=maxIter, init_params=FALSE, verbose=FALSE)
ISTA_res <- runAlg(ISTA)
print(paste0("ISTA loss=", paste0(tail(ISTA_res@loss,1), collapse = "; "), " after ",ISTA_res@iteration," iterations"))

FISTA <- new("FISTA-SLOPE", y=y, A=A, x=x0, lambda=lambdahat, stepSize=stepSize, 
             max_iter=maxIter, verbose=FALSE)
FISTA_res <- runAlg(FISTA)
print(paste0("FISTA loss=", paste0(tail(FISTA_res@loss,1), collapse = "; "), " after ",FISTA_res@iteration," iterations"))

# plotting the results
t <- c(1:100)
df <- data.frame(t, ISTA_res@loss[t], FISTA_res@loss[t], AMP_res@loss[t])
colnames(df) <- c("iteration","ISTA","FISTA","AMP")
df <- melt(df, id="iteration")
ggplot(data = df,
       aes(x=iteration, y=value, colour=variable))+
  geom_line(size=.5) +
  ggtitle("Loss evolution")