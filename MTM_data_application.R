##  Copyright 2023 Annamaria Guolo (University of Padova) 
##  Permission to use, copy, modify and distribute this software and
##  its documentation, for any purpose and without fee, is hereby granted,
##  provided that:
##  1) this copyright notice appears in all copies and
##  2) the source is acknowledged through a citation to the paper
##     Guolo A. (2023). Hierarchical multinomial processing tree models for meta-analysis of diagnostic accuracy studies. arXiv.
##  The Authors make no representation about the suitability of this software
##  for any purpose.  It is provided "as is", without express or implied warranty
###########################################

## APPLICATION

## Confusion assessment data from
## Shi et al. (Neuropsych Dis Treat, 2013)
###########################################

source('MTM_software.R')
dyn.load('MTM_auxiliary_functions.so')

## nodes and weights for 21 points Gauss-Hermite quadrature
objGH <- gauss.quad(21, 'hermite') 

## data
TP <- c(21,35,77,16,27,22,33,26,19,22,137,225,14,9,15,15,39,15,56,34)
FP <- c(4,2,0,3,0,0,3,8,0,2,11,12,3,2,0,0,13,0,6,0)
FN <- c(3,40,3,1,3,3,13,6,6,2,42,60,62,12,2,2,16,6,5,3)
TN <- c(43,104,16,80,93,19,70,41,75,76,36,706,344,131,71,37,64,58,59,29)

P <- TP + FN
N <- TN + FP
n <- length(TP)
study <- 1:length(TP)
data.bin <- data.frame('study'=study, 'TP'=TP, 'FP'=FP, 'FN'=FN, 'TN'=TN, 'P'=P, 'N'=N)

## compute the observed SP and SE, logit link
sp.obs <- se.obs <- eta.obs <- xi.obs <- var.eta <- var.xi <- rep(NA, n)
for(i in 1:n){
    sp.obs[i] <- 1-data.bin$FP[i]/data.bin$N[i]
    se.obs[i] <- data.bin$TP[i]/data.bin$P[i]
    var.eta[i] <- 1/(data.bin$TP[i]) + 1/(data.bin$P[i]-data.bin$TP[i])
    var.xi[i] <- 1/(data.bin$FP[i]) + 1/(data.bin$N[i]-data.bin$FP[i])
    if(sp.obs[i]==1 | sp.obs[i]==0){ ## continuity  correction, if necessary
        sp.obs[i] <- 1 - (data.bin$FP[i]+0.5)/(data.bin$N[i]+1)
        var.xi[i] <- 1/(data.bin$FP[i]+0.5) + 1/(data.bin$N[i]-data.bin$FP[i]+0.5)
    }
    if(se.obs[i]==1 | se.obs[i]==0){ ## continuity  correction, if necessary
        se.obs[i] <- (data.bin$TP[i]+0.5)/(data.bin$P[i]+1)
        var.eta[i] <- 1/(data.bin$TP[i]+0.5) + 1/(data.bin$P[i]-data.bin$TP[i]+0.5)
    }
    eta.obs[i] <- log(se.obs[i]/(1-se.obs[i]))  ## logit of SE
    xi.obs[i] <-log( sp.obs[i]/(1-sp.obs[i]) ) ## logit of SP  
}
data.all <- data.frame(data.bin, eta.obs, xi.obs, var.eta, var.xi)

## Approximate likelihood  

theta.start.lik <- c('eta.mean'=mean(data.all[,'eta.obs']),
                 'xi.mean'=mean(data.all[,'xi.obs']),
                 'eta.var'=var(data.all[,'eta.obs']),
                 'xi.var'=var(data.all[,'xi.obs']),
                 'rho'=cor(data.all[,'eta.obs'], data.all[,'xi.obs']))
ris.lik <- lik(theta.start.lik, data.approx = data.all, hessian=TRUE)
   

## MTM, logit link
## estimate of the prevalences
prev <- rlik.MTM(data.all)
ris.mtm.logit <- lik.MTM(theta.start.lik, data.bin=data.all, objGH=objGH, method='Nelder-Mead', link='logit')

## compute SE and SP (mle + s.e.), using delta method
f.SE <- function(par, var.matrix){
    fn.to.deriv <- function(x)
        plogis(x[1])
    est <- plogis(par[1])
    J <- matrix(grad(fn.to.deriv, par), nrow=1)
    variance <- J%*%var.matrix%*%t(J)
    return(c(est, sqrt(variance)))
}

f.SP <- function(par, var.matrix){
    fn.to.deriv <- function(x)
        1-plogis(x[2])
    est <- 1-plogis(par[2])
    J <- matrix(grad(fn.to.deriv, par), nrow=1)
    variance <- J%*%var.matrix%*%t(J)
    return(c(est, sqrt(variance)))
}

round(f.SE(ris.lik$theta.hat, ris.lik$theta.sand), 2)
round(f.SP(ris.lik$theta.hat, ris.lik$theta.sand) *c(-1,1) + c(1,0), 2)

round(f.SE(ris.mtm.logit$theta.hat, ris.mtm.logit$theta.sand), 2)
round(f.SP(ris.mtm.logit$theta.hat, ris.mtm.logit$theta.sand) *c(-1,1) + c(1,0), 2)

