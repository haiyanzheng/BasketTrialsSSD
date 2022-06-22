library("reshape2")

source("subfun_and_simulators.R")

nsim = 100000
MytargEff = -0.4
#------------------------------------------------------------------------------#
#                        True positive & False positive                        #
#------------------------------------------------------------------------------#
MyVar = seq(0.2, 0.5, by = 0.05)
MyNk = rep(0, length(MyVar))


HDistmat = array(0, dim = c(7, 7))  # 7 subgroups

for(i in 1:length(MyVar)){
  MyNk[i] = nlm(MyBrwfun, Ri = 0.5, sig02 = rep(MyVar[i], 7), s02 = 100, wiq = HDistmat, 
                cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = MytargEff, 
                eta = 0.95, zeta = 0.80,
                rep(0.0001, 7))$estimate[1]
}


set.seed(123)
for(i in 1:length(MyVar)){
assign(paste0("Sc4a", i), lapply(1:nsim, function(x){
  if(x%%10000 == 0) print(c(i, x))
  SimBrw(Sc = cbind(rep(-0.4, 7), rep(MyVar[i], 7)), wMatr = HDistmat, Ni = rep(MyNk[i], 7), 
         effsize = MytargEff, Ri = 0.5, opprior = c(0, 100), dw = c(1.1, 1.1), br = c(54, 3), 
         cctrpar = 0.05)
}))
}


set.seed(123)
for(i in 1:length(MyVar)){
assign(paste0("Sc6a", i), lapply(1:nsim, function(x){
  if(x%%10000 == 0) print(c(i, x))
  SimBrw(Sc = cbind(rep(0, 7), rep(MyVar[i], 7)), wMatr = HDistmat, Ni = rep(MyNk[i], 7), 
         effsize = MytargEff, Ri = 0.5, opprior = c(0, 100), dw = c(1.1, 1.1), br = c(54, 3), 
         cctrpar = 0.05)
}))
}

#------------------------------------------------------------------------------#

for(i in 1:length(MyVar)){
  assign(paste0("Sc4BrwTP", i), array(0, dim = c(nsim, 7)))
  assign(paste0("Sc6BrwFP", i), array(0, dim = c(nsim, 7)))
}


for(i in 1:nsim){
  Sc4BrwTP1[i,] = unlist(Sc4a1[[i]][1])
  Sc4BrwTP2[i,] = unlist(Sc4a2[[i]][1])
  Sc4BrwTP3[i,] = unlist(Sc4a3[[i]][1])
  Sc4BrwTP4[i,] = unlist(Sc4a4[[i]][1])
  Sc4BrwTP5[i,] = unlist(Sc4a5[[i]][1])
  Sc4BrwTP6[i,] = unlist(Sc4a6[[i]][1])
  Sc4BrwTP7[i,] = unlist(Sc4a7[[i]][1])
  
  Sc6BrwFP1[i,] = unlist(Sc6a1[[i]][1])
  Sc6BrwFP2[i,] = unlist(Sc6a2[[i]][1])
  Sc6BrwFP3[i,] = unlist(Sc6a3[[i]][1])
  Sc6BrwFP4[i,] = unlist(Sc6a4[[i]][1])
  Sc6BrwFP5[i,] = unlist(Sc6a5[[i]][1])
  Sc6BrwFP6[i,] = unlist(Sc6a6[[i]][1])
  Sc6BrwFP7[i,] = unlist(Sc6a7[[i]][1])
}

  
subTP = array(0, dim = c(7, length(MyVar)))
subFP = array(0, dim = c(7, length(MyVar)))

for(i in 1:length(MyVar)){
  subTP[,i] = colSums(get(paste0("Sc4BrwTP", i)) >= 0.95)/nsim
  subFP[,i] = colSums(get(paste0("Sc6BrwFP", i)) >= 0.95)/nsim
}

colnames(subTP) <- paste0("TPvar", 1:length(MyVar))
colnames(subFP) <- paste0("TPvar", 1:length(MyVar))

subTPtr <- melt(subTP)
subFPtr <- melt(subFP)

subTPtr$sig2k <- rep(MyVar, each = 7)
subFPtr$sig2k <- rep(MyVar, each = 7)


#------------------------------------------------------------------------------#

TPfig = ggplot(data = subTPtr, aes(x = sig2k, y = value*100)) + 
  geom_hline(yintercept = 80, linetype = 2, colour = "gray") + 
  geom_line() + geom_point() + theme_bw() +
  scale_y_continuous(name = "True positive (%)", limits = c(60, 100)) +
  scale_x_continuous(name = bquote(sigma[k]^2),  breaks = MyVar,
                     sec.axis = sec_axis(~.*29.51824, name = bquote(n[k]),
                                         breaks = c(5.9, 7.4, 8.9, 10.3, 11.8, 13.3, 14.8)))

FPfig = ggplot(data = subFPtr, aes(x = sig2k, y = value*100)) + 
  geom_hline(yintercept = 5, linetype = 2, colour = "gray") + 
  geom_line() + geom_point() + theme_bw() +
  scale_y_continuous(name = "False positive (%)", limits = c(0, 10)) +
  scale_x_continuous(name = bquote(sigma[k]^2),  breaks = MyVar,
                     sec.axis = sec_axis(~.*29.51824, name = bquote(n[k]),
                                         breaks = c(5.9, 7.4, 8.9, 10.3, 11.8, 13.3, 14.8)))


# eps: 600*350
cowplot::plot_grid(TPfig, FPfig, align = 'v',
                   labels = c("(i)", "(ii)"),
                   hjust = -0.5,
                   ncol = 2
)
