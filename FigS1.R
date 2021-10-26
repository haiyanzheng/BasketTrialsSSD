require(gtools)

source("subfun_and_simulators.R")

#------------------------ figure 2a ------------------------#
# delta*/sigma2 = c(0.4/0.25, 0.4/0.2, 0.4/0.15, 0.4/0.1)
SSeta = array(0, dim = c(4, 6, 5))
MySig = c(0.25, 0.2, 0.15, 0.1)
MyEta = c(0.875, 0.90, 0.925, 0.95, 0.975)


for(p in 1:5){
  for(i in 1:4){
    SSeta[i, 1, p] = sum(sapply(rep(MySig[i], 7), function(vari){
      NoBrwNi(vari, Ri = 1/2, eta = MyEta[p], zeta = 0.80, targEff = 0.4)
    }))
    
    SSeta[i, 2, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                             wiq = diag(x = -1, nrow = 7, ncol = 7) + 1, 
                             cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                             eta = MyEta[p], zeta = 0.80,
                             rep(0.0001, 7))$estimate)
    
    SSeta[i, 3, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                             wiq = diag(x = -0.3, nrow = 7, ncol = 7) + 0.3, 
                             cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                             eta = MyEta[p], zeta = 0.80,
                             rep(0.0001, 7))$estimate)
    
    SSeta[i, 4, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                             wiq = diag(x = -0.1, nrow = 7, ncol = 7) + 0.1, 
                             cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                             eta = MyEta[p], zeta = 0.80,
                             rep(0.0001, 7))$estimate)
    
    SSeta[i, 5, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                             wiq = diag(x = -0.05, nrow = 7, ncol = 7) + 0.05, 
                             cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                             eta = MyEta[p], zeta = 0.80,
                             rep(0.0001, 7))$estimate)
    
    SSeta[i, 6, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                             wiq = diag(x = 0, nrow = 7, ncol = 7), 
                             cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                             eta = MyEta[p], zeta = 0.80,
                             rep(0.0001, 7))$estimate)
  }
}

for(p in 1:5){
  for(k in 1:6){
    assign(paste0("cMyEta", p, "Brw", k), 
           data.frame(Case = c("Treff1", "Treff2",
                               "Treff3", "Treff4"),
                      SS = SSeta[,k,p],
                      Appr = paste0("Approach ", k), 
                      Alp = paste0("eta = ", MyEta[p])),
    )
  }
}

mats <- grep(x= ls(pos=1), pattern="cMyEta", value=TRUE)


Fig1SSeta <- do.call(rbind, mget(mixedsort(mats)))

Fig1SSeta$Case <- factor(Fig2SSeta$Case, 
                         labels = c(bquote(sigma["k"]^2 ~"= 0.25"),
                                    bquote(sigma["k"]^2 ~"= 0.20"),
                                    bquote(sigma["k"]^2 ~"= 0.15"),
                                    bquote(sigma["k"]^2 ~"= 0.10"))
)

m1a <- ggplot(data = Fig2SSeta, aes(x = Alp, y = SS, group = Appr)) + theme_bw() +
  geom_line(aes(color = Appr), size = 0.8) + 
  geom_point(aes(shape = Appr, color = Appr), size = 2.1) + 
  scale_x_discrete(expression(eta), labels = MyEta) + 
  scale_y_continuous(name = "Total sample size", 
                     breaks = seq(0, 350, by = 70), limits = c(0, 350)) +
  # scale_colour_discrete(labels = c("Commensurability", "Normalised weight")) +
  theme(legend.position='none', legend.title = element_blank())

fig1a <- m2a + scale_color_manual(name = element_blank(), 
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
                     values = c(13, 11, 16, 17, 0, 1)) + 
  facet_wrap(~Case, ncol = 4, labeller = label_parsed) 


#------------------------ figure 2b ------------------------#
# delta*/sigma2 = c(0.4/0.25, 0.4/0.2, 0.4/0.15, 0.4/0.1)
SSzeta = array(0, dim = c(4, 6, 5))  # dim = c(#MySig, #wiq, #MyZeta)
MySig = c(0.25, 0.2, 0.15, 0.1)
MyZeta = c(0.75, 0.80, 0.85, 0.90, 0.95)


for(p in 1:5){
  for(i in 1:4){
    SSzeta[i, 1, p] = sum(sapply(rep(MySig[i], 7), function(vari){
      NoBrwNi(vari, Ri = 1/2, eta = 0.95, zeta = MyZeta[p], targEff = 0.4)
    }))
    
    SSzeta[i, 2, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                              wiq = diag(x = -1, nrow = 7, ncol = 7) + 1, 
                              cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                              eta = 0.95, zeta = MyZeta[p],
                              rep(0.0001, 7))$estimate)
    
    SSzeta[i, 3, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                              wiq = diag(x = -0.3, nrow = 7, ncol = 7) + 0.3, 
                              cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                              eta = 0.95, zeta = MyZeta[p],
                              rep(0.0001, 7))$estimate)
    
    SSzeta[i, 4, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                              wiq = diag(x = -0.1, nrow = 7, ncol = 7) + 0.1, 
                              cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                              eta = 0.95, zeta = MyZeta[p],
                              rep(0.0001, 7))$estimate)
    
    SSzeta[i, 5, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                              wiq = diag(x = -0.05, nrow = 7, ncol = 7) + 0.05, 
                              cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                              eta = 0.95, zeta = MyZeta[p],
                              rep(0.0001, 7))$estimate)
    
    SSzeta[i, 6, p] = sum(nlm(MyBrwfun, Ri = 1/2, sig02 = rep(MySig[i], 7), s02 = 100,
                              wiq = diag(x = 0, nrow = 7, ncol = 7), 
                              cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = 0.4, 
                              eta = 0.95, zeta = MyZeta[p],
                              rep(0.0001, 7))$estimate)
  }
}

for(p in 1:5){
  for(k in 1:6){
    assign(paste0("cMyZeta", p, "Brw", k), 
           data.frame(Case = c("Treff1", "Treff2",
                               "Treff3", "Treff4"),
                      SS = SSzeta[,k,p],
                      Appr = paste0("Approach ", k), 
                      Alp = paste0("eta = ", MyZeta[p])),
    )
  }
}

mats <- grep(x= ls(pos=1), pattern="cMyZeta", value=TRUE)


Fig1SSzeta <- do.call(rbind, mget(mixedsort(mats)))

Fig1SSzeta$Case <- factor(Fig2SSzeta$Case, 
                          labels = c(bquote(sigma["k"]^2 ~"= 0.25"),
                                     bquote(sigma["k"]^2 ~"= 0.20"),
                                     bquote(sigma["k"]^2 ~"= 0.15"),
                                     bquote(sigma["k"]^2 ~"= 0.10"))
)

m1b <- ggplot(data = Fig1SSzeta, aes(x = Alp, y = SS, group = Appr)) + theme_bw() +
  geom_line(aes(color = Appr), size = 0.8) + 
  geom_point(aes(shape = Appr, color = Appr), size = 2.1) + 
  scale_x_discrete(expression(zeta), labels = MyZeta) + 
  scale_y_continuous(name = "Total sample size", 
                     breaks = seq(0, 500, by = 100), limits = c(0, 500)) +
  theme(legend.position='none', legend.title = element_blank())

fig1b <- m1b + scale_color_manual(name = element_blank(), 
                                  labels = c("No borrowing",
                                             bquote(w["qk"] ~"="~ 1),
                                             bquote(w["qk"] ~"="~ 0.3),
                                             bquote(w["qk"] ~"="~ 0.1),
                                             bquote(w["qk"] ~"="~ 0.05),
                                             bquote(w["qk"] ~"="~ 0)),
                                  values = c("lightgray", "#0072B2", "#f4a582", 
                                             "#ca0020", "#a6dba0", "#008837")) + 
  scale_shape_manual(name = element_blank(), 
                     labels = c("No borrowing",
                                bquote(w["qk"] ~"="~ 1),
                                bquote(w["qk"] ~"="~ 0.3),
                                bquote(w["qk"] ~"="~ 0.1),
                                bquote(w["qk"] ~"="~ 0.05),
                                bquote(w["qk"] ~"="~ 0)), 
                     values = c(13, 11, 16, 17, 0, 1)) + 
  facet_wrap(~Case, ncol = 4, labeller = label_parsed) 

# legend <- cowplot::get_legend(fig2b)

prow <- cowplot::plot_grid( fig1a,
                            fig1b,
                            align = 'v',
                            labels = c("(i)", "(ii)", "(iii)"),
                            hjust = -0.5,
                            ncol = 1
)

# (eps) 800*600
cowplot::plot_grid(legend, prow, ncol = 1, rel_heights = c(0.15, 1))

