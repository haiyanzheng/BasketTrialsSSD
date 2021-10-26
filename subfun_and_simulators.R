# load r packages
library(ggplot2)
library(reshape2)
library(gtools)
library(RColorBrewer)

library(patchwork)

#----------------------------- sub-functions -----------------------------#
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#------------------------------- No borrowing ---------------------------------#
NoBrwfun = function(SS, sig02, Ri, targEff, eta = 0.95, zeta = 0.80, s02 = 100){
  sum(abs(SS*Ri*(1-Ri)/sig02 - (qnorm(eta)+qnorm(zeta))^2/targEff^2 - 1/s02))
}


NoBrwNi = function(vari, Ri, eta, zeta, targEff, s02 = 100){
  vari/(Ri*(1-Ri))*((qnorm(eta)+qnorm(zeta))^2/targEff^2-1/s02)
}

#---------------------- Proposed approach of borrowing ------------------------#
pq = function(wk, s0){
  exp(-wk^2/s0)/sum(exp(-wk^2/s0))
}

# remove the diagnal 
RmDiag = function(Matr){
  diag(Matr) = NA
  return(matrix(Matr[is.na(Matr)==0], 
                nrow = nrow(Matr)-1, ncol = ncol(Matr))
  )
}


MyBrwfun = function(MyNi, Ri, sig02, s02, wiq, cctrpar = 0.1, 
                    dw, br, targEff, eta = 0.95, zeta = 0.80){
  
  cplNi = RmDiag(replicate(length(MyNi), MyNi))
  cplsig02 = RmDiag(replicate(length(sig02), sig02))
  tsWiq = RmDiag(wiq)
  tsPiq = apply(tsWiq, 2, pq, s0 = cctrpar)
  obj = MyNi*Ri*(1-Ri)/sig02 + 1/colSums(tsPiq^2*(1/(1/s02 + cplNi*Ri*(1-Ri)/cplsig02) + 
                                                    tsWiq*dw[2]/(dw[1]-1) + (1-tsWiq)*br[2]/(br[1]-1))) - 
    (qnorm(eta) + qnorm(zeta))^2/targEff^2
  
  return(sum(abs(obj)))
}

#----------------------------- other subfunctions -----------------------------#
# Squared Hellinger distance between any two normal distributions
NormHd2 = function(Norm1, Norm2){
  1 - sqrt(2*sqrt(Norm1[2]*Norm2[2])/(Norm1[2]+Norm2[2]))*exp(-1/4*(Norm1[1]-Norm2[1])^2/(Norm1[2]+Norm2[2]))
}

# computation of the pairwise Hellinger distance
HdMat = function(mu.hist, var.hist){
  histdat = cbind(1:length(mu.hist), mu.hist, var.hist)
  colnames(histdat)[1] = "id"
  
  ## pairwise indicators
  prs = t(combn(1:length(mu.hist), 2))
  colnames(prs) = c("id1", "id2")
  
  myNorm1 = subset(histdat[prs[,1],], select = c("mu.hist", "var.hist"))
  myNorm2 = subset(histdat[prs[,2],], select = c("mu.hist", "var.hist"))
  
  a = numeric(length(mu.hist)*(length(mu.hist)-1)/2)
  for(i in 1:(length(mu.hist)*(length(mu.hist)-1)/2)) a[i] = NormHd2(Norm1 = myNorm1[i,], Norm2 = myNorm2[i,])
  
  myHdMat = cbind(prs, round(sqrt(a), 3))
  
  return(myHdMat)
}

#-------------------------- basket trials simulator ---------------------------#
#- which also produces the analysis result based on the proposed methodology -#
SimBrw = function(Sc, wMatr, Ni, effsize, Ri = 0.5, opprior = c(0, 100), 
                  dw = c(2, 2), br = c(54, 3), cctrpar = 0.1){
  
  respoA = sapply(1:length(Ni), function(x){
    rnorm(Ni[x]/2, mean = Sc[x, 1], sd = sqrt(Sc[x, 2]))
  })
  
  respoB = sapply(1:length(Ni), function(x){
    rnorm(Ni[x]/2, mean = 0, sd = sqrt(Sc[x, 2]))
  })
  
  if(is.array(respoA)){
    smpdiff = colMeans(respoA) - colMeans(respoB)
  } else{
    smpdiff = sapply(respoA, mean)-sapply(respoB, mean) 
  }
  
  cplNi = RmDiag(replicate(length(Ni), Ni))
  cplsig02 = RmDiag(replicate(nrow(Sc), Sc[,2]))
  cplsmpdiff = RmDiag(replicate(length(smpdiff), smpdiff))
  
  etaq = opprior[1]/(1 + opprior[2]/cplsig02*Ni*Ri*(1-Ri)) + 
    cplsmpdiff/(1 + cplsig02/(opprior[2]*cplNi*Ri*(1-Ri)))
  var.etaq = 1/(1/opprior[2] + cplNi*Ri*(1-Ri)/cplsig02)
  
  tsWiq = RmDiag(wMatr)
  tsPiq = apply(tsWiq, 2, pq, s0 = cctrpar)
  
  xi2iq = var.etaq + tsWiq*dw[2]/(dw[1]-1) + (1-tsWiq)*br[2]/(br[1]-1)
  
  dmui = 1/(colSums(tsPiq^2*xi2iq)+Sc[,2]/(Ni*Ri*(1-Ri)))*(colSums(tsPiq^2*xi2iq)*smpdiff + 
                                                             Sc[,2]/(Ni*Ri*(1-Ri))*colSums(tsPiq*etaq))
  sig2mui = 1/(1/colSums(tsPiq^2*xi2iq) + Ni*Ri*(1-Ri)/Sc[,2])
  
  if(effsize > 0){
    ProbRej = pnorm((dmui-0)/sqrt(sig2mui))
    ProbNoRej = pnorm((effsize-dmui)/sqrt(sig2mui))
  } else{
    ProbRej = pnorm((0-dmui)/sqrt(sig2mui))
    ProbNoRej = pnorm((dmui-effsize)/sqrt(sig2mui))
  }
  
  return(list(ProbRej, ProbNoRej))
  
}

#-------------------------- basket trials simulator ---------------------------#
# which also produces the analysis result based on the approach of no borrowing #
SimNoBrw = function(Sc, Ni, effsize, Ri = 0.5, opprior = c(0, 100)){
  
  respoA = sapply(1:length(Ni), function(x){
    rnorm(Ni[x]/2, mean = Sc[x, 1], sd = sqrt(Sc[x, 2]))
  })
  
  respoB = sapply(1:length(Ni), function(x){
    rnorm(Ni[x]/2, mean = 0, sd = sqrt(Sc[x, 2]))
  })
  
  if(is.array(respoA)){
    smpdiff = colMeans(respoA) - colMeans(respoB)
  } else{
    smpdiff = sapply(respoA, mean)-sapply(respoB, mean) 
  }
  
  etai =  opprior[1]/(1 + opprior[2]/Sc[,2]*Ni*Ri*(1-Ri)) + 
    smpdiff/(1 + Sc[,2]/(opprior[2]*Ni*Ri*(1-Ri)))
  
  sig2etai = 1/(1/opprior[2] + Ni*Ri*(1-Ri)/Sc[,2])
  
  if(effsize > 0){
    ProbRej = pnorm((etai-0)/sqrt(sig2etai))
    ProbNoRej = pnorm((effsize-etai)/sqrt(sig2etai))
  } else{
    ProbRej = pnorm((0-etai)/sqrt(sig2etai))
    ProbNoRej = pnorm((etai-effsize)/sqrt(sig2etai))
  }
  
  return(list(ProbRej, ProbNoRej))
  
}

#-------------------------- Simulation scenarios --------------------------#
# Configuration 1: Divergent data across subtrials; taken from the SUMMIT
MySc1 = cbind(c(-0.489, 0.226, -0.181, 0.293, 0.329, -0.275, -0.136),
              c( 0.587, 0.345,  0.380, 0.347, 0.344,  0.392,  0.392)^2)  ## for the experimental arm


# Configuration 2: Consistent data across subtrials
MySc2 = cbind(c(-0.489, -0.226, -0.281, -0.293, -0.329, -0.275, -0.236),
              c( 0.587,  0.345,  0.380,  0.347,  0.344,  0.392,  0.392)^2)


# Configuration 3: Consistent data across subtrials with homoscedasticity
MySc3 = cbind(c(-0.289, -0.226, -0.281, -0.293, -0.329, -0.275, -0.236),
              rep(0.30, 7))


# Configuration 4: Consistent data across subtrials
MySc4 = cbind(rep(-0.4, 7),
              rep(0.30, 7))


# Configuration 5: Mixed null
MySc5 = cbind(c(-0.289, 0,    -0.181, 0,     0,    -0.275,  0),
              c( 0.587, 0.345, 0.380, 0.347, 0.344, 0.392,  0.392)^2)


# Configuration 6: Global null
MySc6 = cbind(rep(0, 7), rep(0.30, 7))
