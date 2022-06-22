source("subfun_and_simulators.R")


#-------------------------- Fig2S: sample sizes -------------------------#
# delta*/sigma2 = c(0.4/0.25, 0.4/0.2, 0.4/0.15, 0.4/0.1)
SSeff = array(0, dim = c(6, 7))  # dim = c(#wiq, #MyEff)
MyEff = c(0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30)


for(p in 1:7){
  
  SSeff[1, p] = sum(sapply(rep(0.25, 7), function(vari){
    NoBrwNi(vari, Ri = 1/2, eta = 0.95, zeta = 0.80, targEff = MyEff[p])
  }))
  
  SSeff[2, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(0.25, 7), s02 = 100,
                        wiq = diag(x = -1, nrow = 7, ncol = 7) + 1, 
                        cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = MyEff[p], 
                        eta = 0.95, zeta = 0.80,
                        rep(0.0001, 7))$estimate)
  
  SSeff[3, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(0.25, 7), s02 = 100,
                        wiq = diag(x = -0.3, nrow = 7, ncol = 7) + 0.3, 
                        cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = MyEff[p], 
                        eta = 0.95, zeta = 0.80,
                        rep(0.0001, 7))$estimate)
  
  SSeff[4, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(0.25, 7), s02 = 100,
                        wiq = diag(x = -0.1, nrow = 7, ncol = 7) + 0.1, 
                        cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = MyEff[p], 
                        eta = 0.95, zeta = 0.80,
                        rep(0.0001, 7))$estimate)
  
  SSeff[5, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(0.25, 7), s02 = 100,
                        wiq = diag(x = -0.05, nrow = 7, ncol = 7) + 0.05, 
                        cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = MyEff[p], 
                        eta = 0.95, zeta = 0.80,
                        rep(0.0001, 7))$estimate)
  
  SSeff[6, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(0.25, 7), s02 = 100,
                        wiq = diag(x = 0, nrow = 7, ncol = 7), 
                        cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = MyEff[p], 
                        eta = 0.95, zeta = 0.80,
                        rep(0.0001, 7))$estimate)
}


for(p in 1:7){
  for(k in 1:6){
    assign(paste0("cMyEff", p, "Brw", k), 
           data.frame(SS = SSeff[k,p],
                      Appr = paste0("Approach ", k), 
                      Alp = paste0("delta = ", MyEff[p])),
    )
  }
}

mats <- grep(x= ls(pos=1), pattern="cMyEff", value=TRUE)
Fig2SSeff <- do.call(rbind, mget(mixedsort(mats)))


m2c <- ggplot(data = Fig2SSeff, aes(x = Alp, y = SS, group = Appr)) + theme_bw() +
  geom_line(aes(color = Appr), size = 0.8) + 
  geom_point(aes(shape = Appr, color = Appr), size = 2.1) + 
  scale_x_discrete(bquote(delta/sigma["k"]^2), labels = MyEff/0.25) + 
  scale_y_continuous(name = "Total sample size", 
                     breaks = seq(0, 500, by = 100), limits = c(0, 500)) +
  # scale_colour_discrete(labels = c("Commensurability", "Normalised weight")) +
  theme(legend.position='right', legend.title = element_blank())

# (eps) 450*350
fig2c <- m2c + scale_color_manual(name = element_blank(), 
                                  labels = c("No borrowing",
                                             bquote(w["qk"] ~"="~ 1),
                                             bquote(w["qk"] ~"="~ 0.3),
                                             bquote(w["qk"] ~"="~ 0.1),
                                             bquote(w["qk"] ~"="~ 0.05),
                                             bquote(w["qk"] ~"="~ 0)),
                                  # values = c("#0072B2", "#31a354", "#56B4E9", "#a1d99b", 
                                  # "#bcbddc", "lightgray", "pink"))  + 
                                  values = c("lightgray", "#0072B2", "#f4a582", 
                                             "#ca0020", "#a6dba0", "#008837")) + 
  scale_shape_manual(name = element_blank(), 
                     labels = c("No borrowing",
                                bquote(w["qk"] ~"="~ 1),
                                bquote(w["qk"] ~"="~ 0.3),
                                bquote(w["qk"] ~"="~ 0.1),
                                bquote(w["qk"] ~"="~ 0.05),
                                bquote(w["qk"] ~"="~ 0)), 
                     values = c(13, 11, 16, 17, 0, 1)) 


fig2c
