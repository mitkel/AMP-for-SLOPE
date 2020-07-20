# # setting proper working dir and sourcing algorithms (works only in RStudio)
# library("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("IterativeAlgs.R")

library(ggplot2, respape2)

set.seed(351759)

p=500;
delta=0.5;
eps=0.1;
sigma=0;
n=p*delta
A=matrix(rnorm(n*p,mean=0,sd=1/sqrt(p*delta)), n,p) 
alpha=c(rep(1,p*delta),rep(0,p*(1-delta)))
stepSize=100;
maxIter=10**3

x=rnorm(p)*rbinom(p,1,eps)
y=as.numeric(A%*%x+sigma*rnorm(n))
lambda <- alpha_to_lambda(alpha, x, delta)

x0=rnorm(p)*rbinom(p,1,eps)

ISTA <- new("ISTA-SLOPE", y=y, A=A, x=x0, lambda=lambda, stepSize=stepSize, max_iter=maxIter)
FISTA <- new("FISTA-SLOPE", y=y, A=A, x=x0, lambda=lambda, stepSize=stepSize, max_iter=maxIter)
AMP <-   new("AMP-SLOPE", y=y, A=A, prior=x, lambda=lambda, alpha=alpha, max_iter=10**2, backtrack=FALSE)

ISTA_res <- runAlg(ISTA)
FISTA_res <- runAlg(FISTA)
AMP_res <- runAlg(AMP)

print(paste0("ISTA: x=", paste0(head(ISTA_res@x,3), collapse = "; "), " after ",ISTA_res@iteration," iterations"))
print(paste0("FISTA: x=", paste0(head(FISTA_res@x,3), collapse = "; "), " after ",FISTA_res@iteration," iterations"))
print(paste0("AMP: x=", paste0(head(AMP_res@x,3), collapse = "; "), " after ",AMP_res@iteration," iterations"))

# plotting the results
t <- c(1:100)
df <- data.frame(t, ISTA_res@loss[t], FISTA_res@loss[t], AMP_res@loss[t])
colnames(df) <- c("iteration","ISTA","FISTA","AMP")
df <- melt(df, id="iteration")
ggplot(data = df,
       aes(x=iteration, y=value, colour=variable))+
  geom_line(size=.5) +
  ggtitle("Loss evolution")
  