source("subfun_and_simulators.R")


#-------------------------- Fig2: trial operating characteristics -------------------------#
# Simulation results for Scenario 1
# For reproducing results for other plots of Fig2
# please set SimSc = MySc2, MySc3, MySc4, MySc5, MySc6, respectively, on line 13

nsim = 100000
MytargEff = -0.4


SimSc = MySc1

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

set.seed(123)
a = lapply(1:nsim, function(x){
  if(x%%1000 == 0) print(x)
  SimBrw(Sc = SimSc, wMatr = HDistmat, Ni = MyBrwNi, effsize = MytargEff, 
         Ri = 0.5, opprior = c(0, 100), dw = c(1.1, 1.1), br = c(54, 3), cctrpar = 0.05)
})


set.seed(123)
b = lapply(1:nsim, function(x){
  if(x%%1000 == 0) print(x)
  SimNoBrw(Sc = SimSc, Ni = MyBrwNi, effsize = MytargEff, Ri = 0.5, opprior = c(0, 100))
})


MyBrwprob1 = array(0, dim = c(nsim, nrow(SimSc)))
MyBrwprob2 = array(0, dim = c(nsim, nrow(SimSc)))

NoBrwprob1 = array(0, dim = c(nsim, nrow(SimSc)))
NoBrwprob2 = array(0, dim = c(nsim, nrow(SimSc)))

for(i in 1:nsim){
  MyBrwprob1[i,] = unlist(a[[i]][1])
  MyBrwprob2[i,] = unlist(a[[i]][2])
  NoBrwprob1[i,] = unlist(b[[i]][1])
  NoBrwprob2[i,] = unlist(b[[i]][2])
}

# Proportion of trials declaring EFFECTIVENESS based on the *proposed methodology*
colSums(MyBrwprob1 >= 0.95)/nsim
# Proportion of trials declaring FUTILITY based on the *proposed methodology*
colSums(MyBrwprob2 >= 0.80)/nsim

# Proportion of trials declaring EFFECTIVENESS based on the approach of *no borrowing*
colSums(NoBrwprob1 >= 0.95)/nsim
# Proportion of trials declaring FUTILITY based on the approach of *no borrowing*
colSums(NoBrwprob2 >= 0.80)/nsim
