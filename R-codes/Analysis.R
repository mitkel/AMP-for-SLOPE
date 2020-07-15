# setting proper working dir and sourcing algorithms (works only in RStudio)
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("IterativeAlgs.R")

library(ggplot2, respape2)


ISTA <- new("ISTA-SLOPE",
            y = c(6,6),
            A = diag(c(3,2)),
            x = c(23,-10),
            theta = c(2,1),
            stepSize = 0.001,
            iteration = 0,
            loss = Inf,
            lossFun = SL1Loss)

FISTA <- new("FISTA-SLOPE",
             y = c(6,6),
             A = diag(c(3,2)),
             x = c(23,-10),
             x_old = c(23,-10),
             u = c(23,-10),
             theta = c(2,1),
             stepSize = 0.001,
             t = 1,
             t_old = 1,
             iteration = 0,
             loss = Inf,
             lossFun = SL1Loss)

ISTA_res <- runAlg(ISTA)
FISTA_res <- runAlg(FISTA)

print(paste0("ISTA: x=", paste0(ISTA_res@x, collapse = "; "), " after ",ISTA_res@iteration," iterations"))
print(paste0("FISTA: x=", paste0(FISTA_res@x, collapse = "; "), " after ",FISTA_res@iteration," iterations"))

# plotting the results
library(ggplot2)
library(reshape2)
t <- c(1:500)
df <- data.frame(t, ISTA_res@loss[t], FISTA_res@loss[t])
colnames(df) <- c("iteration","ISTA","FISTA")
df <- melt(df, id="iteration")
ggplot(data = df,
       aes(x=iteration, y=value, colour=variable))+
  geom_line(size=.5) +
  ggtitle("Loss evolution")
  