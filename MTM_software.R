##  Copyright 2023 Annamaria Guolo (University of Padova) 
##  Permission to use, copy, modify and distribute this software and
##  its documentation, for any purpose and without fee, is hereby granted,
##  provided that:
##  1) this copyright notice appears in all copies and
##  2) the source is acknowledged through a citation to the paper
##     Guolo A. (2023). Hierarchical multinomial processing tree models for meta-analysis of diagnostic accuracy studies. arXiv.
##  The Authors make no representation about the suitability of this software
##  for any purpose.  It is provided "as is", without express or implied warranty



library(nlme)
library(numDeriv)
library(statmod)
library(mvtnorm)
library(VGAM)

## Parameter vector
## theta= c(1 = mu.eta, 2 = mu.xi,
##          3 = var.eta, 4 = var.xi, 
##          5 = rho.etaxi)
## prevalences do not enter theta given the parameter separability

## data.approx: dataset with information about
## eta.obs, xi.obs, var.eta, var.xi, rho.etaxi

## data.bin: dataset with information about
## TP, P, FP, FN, TN, N

## ========================================
## APPROXIMATE MODEL: LIKELIHOOD FUNCTION
## maximum likelihood estimation
## for logit, probit, cloglog link function
## data.approx
## ========================================

fn.lik <- function(theta, data.approx){
    mu.eta <- theta[1]
    mu <- theta[2]
    cond <- (theta[5] > 1 | theta[5] < -1)
    if( theta[3] <= 0 | theta[4] <= 0 | cond==TRUE)
        return(NA)
    sigma2.eta <- theta[3]
    sigma2.xi <- theta[4]
    rho <- theta[5]
    n <- NROW(data.approx)
    f <- 0
    for(i in 1:n){
        eta.i <- data.approx[i,'eta.obs']
        xi.i <- data.approx[i,'xi.obs']
        var.eta <- sigma2.eta + data.approx[i,'var.eta'] 
        var.xi <- sigma2.xi + data.approx[i,'var.xi']
        cov.etaxi <- rho*sqrt(var.eta*var.xi)
        var.matrix <- matrix(c(var.eta, cov.etaxi, cov.etaxi, var.xi), ncol=2, nrow=2)
        f <- f + dmvnorm(c(data.approx[i,'eta.obs'], data.approx[i,'xi.obs']), mean=c(mu.eta, mu), sigma=var.matrix, log=TRUE)
    }
    return(f)
}

## maximize the likelihood function
lik <- function(theta.start.lik, data.approx, method='Nelder-Mead', hessian=TRUE){
    ans <- list()
    ris <- optim(theta.start.lik, fn.lik, data.approx=data.approx, control=list(maxit=50000, fnscale=-1), hessian=hessian, method=method)
    ## mle
    theta.hat <- ris$par
    ##standard error
    theta.sand <- sand.lik(theta.hat, data.approx)
    ans$theta.hat <- theta.hat
    ans$theta.var <- solve(-ris$hessian)
    ans$theta.se <- sqrt(diag(solve(-ris$hessian)))
    ans$convergence <- ris$convergence
    ans$theta.sand <- theta.sand
    ans$theta.sand.se  <- sqrt(diag(theta.sand))
    return(ans)
}

## compute the sandwhich standard error
sand.lik <- function(theta, this.data){
    a <- fdHess(theta, fn.lik, data.approx=this.data[1,])
    values.gradient <- a$gradient  
    G <- values.gradient%*%t(values.gradient)
        H <- a$Hessian    
    for(i in 2:nrow(this.data)){
        a <- fdHess(theta, fn.lik, data.approx=this.data[i,])            
        values.gradient <- a$gradient 
        G <- G + values.gradient%*%t(values.gradient)
        H <- H + a$Hessian      
    }
    return(solve(H)%*%G%*%solve(H))
}

#########################################################
## MTM with random-effects with separable parameters
## LOGIT, PROBIT, or CLOGLOG link
#########################################################

## Compute the mle for all the parameters but the prevalences
fn.MTM.logit <- function(theta, data.bin, objGH){
    p <- length(theta)
    n <- nrow(data.bin)
    cond <- any(theta[(3):(4)] < 0) | theta[(5)] < -1 | theta[(5)] > 1
    if(cond)
        return(NA)
    nodes1 <- nodes2 <- objGH$nodes
    weights1 <- weights2 <- objGH$weights
    num <- length(nodes1)
    f <- .C("int_MTM_logit",
            as.double(theta),
            as.integer(n),
            as.double(data.bin[,'P']),
            as.double(data.bin[,'TP']), 
            as.double(data.bin[,'N']),
            as.double(data.bin[,'FP']), 
            as.double(data.bin[,'FN']),  
            as.double(data.bin[,'TN']), 
            as.double(nodes1),
            as.double(nodes2),
            as.double(weights1),
            as.double(weights2),
            as.integer(num),
            lik=as.double(0.0))$lik
    return(f)
}


fn.MTM.probit <- function(theta, data.bin, objGH){
    p <- length(theta)
    n <- nrow(data.bin)
    cond <- any(theta[(3):(4)] < 0) | theta[(5)] < -1 | theta[(5)] > 1
    if(cond)
        return(NA)
    nodes1 <- nodes2 <- objGH$nodes
    weights1 <- weights2 <- objGH$weights
    num <- length(nodes1)
    f <- .C("int_MTM_probit",
            as.double(theta),
            as.integer(n),
            as.double(data.bin[,'P']),
            as.double(data.bin[,'TP']),  
            as.double(data.bin[,'N']),
            as.double(data.bin[,'FP']),  
            as.double(data.bin[,'FN']),  
            as.double(data.bin[,'TN']),  
            as.double(nodes1),
            as.double(nodes2),
            as.double(weights1),
            as.double(weights2),
            as.integer(num),
            lik=as.double(0.0))$lik
    return(f)
}


fn.MTM.cloglog <- function(theta, data.bin, objGH){
    p <- length(theta)
    n <- nrow(data.bin)
    cond <- any(theta[(3):(4)] < 0) | theta[(5)] < -1 | theta[(5)] > 1
    if(cond)
        return(NA)
    nodes1 <- nodes2 <- objGH$nodes
    weights1 <- weights2 <- objGH$weights
    num <- length(nodes1)
    f <- .C("int_MTM_cloglog",
            as.double(theta),
            as.integer(n),
            as.double(data.bin[,'P']),
            as.double(data.bin[,'TP']),  
            as.double(data.bin[,'N']),
            as.double(data.bin[,'FP']),  
            as.double(data.bin[,'FN']), 
            as.double(data.bin[,'TN']),  
            as.double(nodes1),
            as.double(nodes2),
            as.double(weights1),
            as.double(weights2),
            as.integer(num),
            lik=as.double(0.0))$lik
    return(f)
}

lik.MTM <- function(theta, data.bin, objGH, prev.GS, method='Nelder-Mead', link='logit'){
    if(link=='logit'){
        ris <- optim(theta, fn.MTM.logit, data.bin=data.bin, objGH=objGH, control=list(maxit=50000, fnscale=-1), method=method)
        theta.hat <- ris$par
        theta.var <- solve(-fdHess(theta.hat, fn.MTM.logit, data.bin=data.bin, objGH=objGH)$Hessian)
        theta.sand <- sand.MTM(theta.hat, data.bin=data.bin, objGH=objGH)
    }
    if(link=='probit'){
        ris <- optim(theta, fn.MTM.probit, data.bin=data.bin, objGH=objGH, control=list(maxit=50000, fnscale=-1), method=method)
        theta.hat <- ris$par
        theta.var <- solve(-fdHess(theta.hat, fn.MTM.probit, data.bin=data.bin, objGH=objGH)$Hessian)
        theta.sand <- sand.MTM(theta.hat, data.bin=data.bin, objGH=objGH, link='probit')
    }
   if(link=='cloglog'){
        ris <- optim(theta, fn.MTM.cloglog, data.bin=data.bin, objGH=objGH, control=list(maxit=50000, fnscale=-1), method=method)
        theta.hat <- ris$par
        theta.var <- solve(-fdHess(theta.hat, fn.MTM.cloglog, data.bin=data.bin, objGH=objGH)$Hessian)
        theta.sand <- sand.MTM(theta.hat, data.bin=data.bin, objGH=objGH, link='cloglog')
    }
    ## exit with mle, s.e., sandwich, vcov matrix, convergence
    return(list(theta.hat=theta.hat,
                theta.se=sqrt(diag(theta.var)),
                theta.var=theta.var,
                theta.sand=theta.sand,
                theta.sand.se =sqrt(diag(theta.sand)),
                convergence=ris$convergence))
}

## Restricted mle for prevalences
rlik.MTM <- function(data.bin){
    n <- NROW(data.bin)
    ni <- data.bin$P+data.bin$N
    prev.hat <- data.bin$P/ni
    prev.se <- 1/sqrt(ni^3/(data.bin$P*data.bin$N))
    ## exit with mle, s.e., and vcov martix
    return(list(prev.hat=prev.hat,
                prev.se=prev.se))
}


## standard error based on sandwich for MTM procedure
sand.MTM <- function(theta, data.bin, objGH, link='logit'){
    if(link=='logit')
        a <- fdHess(theta, fn.MTM.logit, data.bin=data.bin[1,], objGH=objGH)
    if(link=='probit')
        a <- fdHess(theta, fn.MTM.probit, data.bin=data.bin[1,], objGH=objGH)
    if(link=='cloglog')
        a <- fdHess(theta, fn.MTM.cloglog, data.bin=data.bin[1,], objGH=objGH)    
    values.gradient <- a$gradient   
    G <- values.gradient%*%t(values.gradient)
    H <- a$Hessian     
    for(i in 2:nrow(data.bin)){
        if(link=='logit')
            a <- fdHess(theta, fn.MTM.logit, data.bin=data.bin[i,], objGH=objGH)
        if(link=='probit')
            a <- fdHess(theta, fn.MTM.probit, data.bin=data.bin[i,], objGH=objGH)
        if(link=='cloglog')
            a <- fdHess(theta, fn.MTM.cloglog, data.bin=data.bin[i,], objGH=objGH)
        values.gradient <- a$gradient 
        G <- G + values.gradient%*%t(values.gradient)
        H <- H + a$Hessian       
    }
    return(solve(H)%*%G%*%solve(H))
}
