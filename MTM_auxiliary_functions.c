/*
####################################################################
  Copyright 2023 Annamaria Guolo (University of Padova) 
  Permission to use, copy, modify and distribute this software and
  its documentation, for any purpose and without fee, is hereby granted,
  provided that:
  1) this copyright notice appears in all copies and
  2) the source is acknowledged through a citation to the paper
     Guolo A. (2023). Hierarchical multinomial processing tree models for meta-analysis of diagnostic accuracy studies. arXiv.
  The Authors make no representation about the suitability of this software
  for any purpose.  It is provided "as is", without express or implied warranty
####################################################################
 * Parameter vector
 theta= c(0 = mu.eta, 1 = mu.xi, 2 = var.eta, 3 = var.xi, 4 = rho.etaxi)
 prevalences do not enter theta given the parameter separability

 * data.bin: dataset with information about
 * TP=x11, P=n1, FP=x01, FN=x10, TN=x00, N=n0
 * nodes and weights as in the chosen Gauss-Hermite quadrature
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>

void int_MTM_logit(const double *theta,
	       const int *n,
	       const double *n1,
	       const double *x11,
	       const double *n0,
	       const double *x01,
	       const double *x10,
	       const double *x00,
	       const double *nodes1,
	       const double *nodes2,
	       const double *weights1,
	       const double *weights2,
	       const int *num,
	       double *lik)
{
  int i,l,j;
  double ftmp, se, sp, a;
  
  for(i=0; i<(*n); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	sp = exp(nodes2[l])/(1+exp(nodes2[l]));
	se = exp(nodes1[j])/(1+exp(nodes1[j]));
	a = pow(se, x11[i])*
	  pow(1-se,x10[i])*
	  pow(1-sp, x01[i])*
	  pow(sp, x00[i]);	  
	ftmp += a *
	  dnorm(nodes2[l], theta[(1)], sqrt(theta[(3)]), 0) *
	  dnorm(nodes1[j], theta[0]+theta[(4)]*sqrt(theta[(2)]/theta[(3)])*(nodes2[l]-theta[(1)]), sqrt(theta[(2)]*(1-theta[(4)]*theta[(4)])), 0) *
  	  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * weights1[j] * weights2[l];
      }
    }
    (*lik) += log(ftmp);
  }
}

void int_MTM_probit(const double *theta,
	       const int *n,
	       const double *n1,
	       const double *x11,
	       const double *n0,
	       const double *x01,
	       const double *x10,
	       const double *x00,
	       const double *nodes1,
	       const double *nodes2,
	       const double *weights1,
	       const double *weights2,
	       const int *num,
	       double *lik)
{
  int i,l,j;
  double ftmp, se, sp, a;
  
  for(i=0; i<(*n); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	sp = pnorm(nodes2[l], 0, 1, 1, 0);   
	se = pnorm(nodes1[j], 0, 1, 1, 0);
	a = pow(se, x11[i])*
	  pow(1-se,x10[i])*
	  pow(1-sp, x01[i])*
	  pow(sp, x00[i]);	  
	ftmp += a *
	  dnorm(nodes2[l], theta[(1)], sqrt(theta[(3)]), 0) *
	  dnorm(nodes1[j], theta[0]+theta[(4)]*sqrt(theta[(2)]/theta[(3)])*(nodes2[l]-theta[(1)]), sqrt(theta[(2)]*(1-theta[(4)]*theta[(4)])), 0) *
  	  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * weights1[j] * weights2[l];
      }
    }
    (*lik) += log(ftmp);
  }
}

void int_MTM_cloglog(const double *theta,
	       const int *n,
	       const double *n1,
	       const double *x11,
	       const double *n0,
	       const double *x01,
	       const double *x10,
	       const double *x00,
	       const double *nodes1,
	       const double *nodes2,
	       const double *weights1,
	       const double *weights2,
	       const int *num,
	       double *lik)
{
  int i,l,j;
  double ftmp, se, sp, a;
  
  for(i=0; i<(*n); i++){
    ftmp = 0.0;
    for(l=0; l<(*num); l++){
      for(j=0; j<(*num); j++){
	sp = 1- exp(-exp(nodes2[l]));   
	se = 1- exp(-exp(nodes1[j]));
	a = pow(se, x11[i])*
	  pow(1-se,x10[i])*
	  pow(1-sp, x01[i])*
	  pow(sp, x00[i]);	  
	ftmp += a *
	  dnorm(nodes2[l], theta[(1)], sqrt(theta[(3)]), 0) *
	  dnorm(nodes1[j], theta[0]+theta[(4)]*sqrt(theta[(2)]/theta[(3)])*(nodes2[l]-theta[(1)]), sqrt(theta[(2)]*(1-theta[(4)]*theta[(4)])), 0) *
  	  exp(nodes1[j]*nodes1[j]) * exp(nodes2[l]*nodes2[l]) * weights1[j] * weights2[l];
      }
    }
    (*lik) += log(ftmp);
  }
}

