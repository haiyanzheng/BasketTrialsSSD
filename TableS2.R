source("subfun_and_simulators.R")


#-------------------------- FigS3: sensitivity analysis -------------------------#
# Simulation results for Scenario 1
# For reproducing other plots of FigS3
# please set w

nsim = 100000
MytargEff = -0.4


SimSc = MySc1
# SimSc = MySc6

HDistmat = array(0, dim = c(nrow(SimSc), nrow(SimSc)))
MyHD = HdMat(mu.hist = SimSc[, 1], var.hist = SimSc[, 2])

for(i in 1:nrow(MyHD)){
  HDistmat[MyHD[i,1], MyHD[i,2]] = MyHD[i,3]
  HDistmat[MyHD[i,2], MyHD[i,1]] = MyHD[i,3]
}

MyBrwFit = nlm(MyBrwfun, Ri = 0.5, sig02 = SimSc[,2], s02 = 100, wiq = HDistmat, 
               cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = -0.4, 
               eta = 0.95, zeta = 0.80,
               rep(0.0001, nrow(SimSc)))

MyBrwNi = MyBrwFit$estimate; MyBrwNi

#--------------------- misspecified w_qk ---------------------#
# AnaWqk1 = array(0.05, dim = c(nrow(SimSc), nrow(SimSc)))
# diag(AnaWqk2) <- 0

AnaWqk2 = array(0.1, dim = c(nrow(SimSc), nrow(SimSc)))
diag(AnaWqk2) <- 0

AnaWqk3 = array(0.3, dim = c(nrow(SimSc), nrow(SimSc)))
diag(AnaWqk3) <- 0

AnaWqk4 = array(0.5, dim = c(nrow(SimSc), nrow(SimSc)))
diag(AnaWqk4) <- 0
#-------------------------------------------------------------#
# substitute AnaWqk1, ..., AnaWqk4 in the argument of wMatr 
set.seed(123)
a = lapply(1:nsim, function(x){
  if(x%%1000 == 0) print(x)
  SimBrw(Sc = SimSc, wMatr = AnaWqk4, Ni = MyBrwNi, effsize = MytargEff, 
         Ri = 0.5, opprior = c(0, 100), dw = c(1.1, 1.1), br = c(54, 3), cctrpar = 0.05)
})


MyBrwprob1 = array(0, dim = c(nsim, nrow(SimSc)))
MyBrwprob2 = array(0, dim = c(nsim, nrow(SimSc)))


for(i in 1:nsim){
  MyBrwprob1[i,] = unlist(a[[i]][1])
  MyBrwprob2[i,] = unlist(a[[i]][2])
}

# Proportion of trials declaring EFFECTIVENESS based on the proposed methodology
round(colSums(MyBrwprob1 >= 0.95)/nsim*100, 1)
# Proportion of trials declaring FUTILITY based on the proposed methodology
round(colSums(MyBrwprob2 >= 0.80)/nsim*100, 1)

