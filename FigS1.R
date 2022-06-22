source("subfun_and_simulators.R")


#-------------------------- Fig1: subtrial sample sizes -------------------------#
MySc01 = cbind(c(0.489, -0.226),
               c(0.587,  0.345)^2)       # Divergent

MySc02 = cbind(c(-0.489, -0.181),
               c( 0.587,  0.380)^2)       # Consistent


SimSc = MySc01

# Sample size computed based on the approach of no borrowing
round(nlm(NoBrwfun, Ri = 0.5, sig02 = SimSc[,2], targEff = 0.4, 
          eta = 0.95, zeta = 0.80, s02 = 100, 
          rep(0.0001, nrow(SimSc)))$estimate, 1)

NoBrwNk = sapply(SimSc[,2], function(vari){
  NoBrwNi(vari, Ri = 0.5, eta = 0.95, zeta = 0.80, targEff = 0.4, s02 = 100)
})

round(NoBrwNk, 1)


w12 = c(0, 0.1, 0.3, 0.5, 1)
w21 = c(0, 0.1, 0.3, 0.5, 1)
data = expand.grid(w12 = w12, w21 = w21)

MyWiq = array(0, dim = c(nrow(SimSc), nrow(SimSc)))


data$SS1 = numeric(nrow(data))
data$SS2 = numeric(nrow(data))

for(i in 1:nrow(data)){
  
  MyWiq[1, 2] = data[i, 1]
  MyWiq[2, 1] = data[i, 2]
  
  data[i, 3:4] = nlm(MyBrwfun, Ri = 0.5, sig02 = SimSc[,2], s02 = 100, 
                     wiq = MyWiq, cctrpar = 0.05, dw = c(2, 2), br = c(54, 3), 
                     targEff = 0.4, eta = 0.95, zeta = 0.80,
                     rep(0.0001, nrow(SimSc)))$estimate
  
}

data$total = data$SS1 + data$SS2

sbtSS = melt(data, id = c("w12", "w21"))
sbtSS$variable <- factor(sbtSS$variable, 
                         labels = c("Subtrial 1", "Subtrial 2", "Total")
)


ggplot(sbtSS, aes(factor(w12), factor(w21), fill = value)) + geom_tile(color = "black") +
  scale_x_discrete(bquote(w[12]), breaks = w12) + scale_y_discrete(bquote(w[21]), breaks = w21) + 
  scale_fill_gradient2(low = "white", high = "white", mid = "white",
                       midpoint = 0, limit = c(5, 75), space = "Lab",
                       name = "Sample size") +
  geom_text(aes(x = factor(w12), y = factor(w21), label = round(value, 1)), color = "black", size = 3) +
  facet_wrap(~variable, ncol = 3) + theme_bw() + theme(legend.position = "none")    
