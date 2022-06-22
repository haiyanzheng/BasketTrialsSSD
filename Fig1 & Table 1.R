source("subfun_and_simulators.R")


#-------------------------- Fig1: simulation scenarios -------------------------#
# Visualisation of the outcome distributions and pairwise commensurability for Scenario 1
# For reproducing other plots of Fig1
# please set SimSc = MySc2, MySc3, MySc4, MySc5, MySc6, respectively, on line 9

SimSc = MySc1

MyDatE = data.frame(mod = paste0("Subtrial ", 1:7), SimSc, 
                    upper = SimSc[,1] + 1.96*sqrt(SimSc[,2]), 
                    lower = SimSc[,1] - 1.96*sqrt(SimSc[,2]),
                    grp = "E")

MyDatC = data.frame(mod = paste0("Subtrial ", 1:7), 
                    X1 = rep(0, 7), X2 = SimSc[,2],
                    upper = 0 + 1.96*sqrt(SimSc[,2]),
                    lower = 0 - 1.96*sqrt(SimSc[,2]),
                    grp = "C")

MyDat = rbind(MyDatE, MyDatC)

levels(MyDat$grp)
MyDat$grp <- factor(MyDat$grp, levels = c("C", "E"))

trtp <- ggplot(data = MyDat, aes(x = X1, y = mod, xmin = lower, xmax = upper, colour = grp)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "darkgray") +
  geom_vline(xintercept = -0.4, linetype = "dashed", colour = "gray") +
  geom_errorbarh(height = 0.1, position = position_dodge(width=0.5)) +
  geom_point(aes(shape = grp), position = position_dodge(width=0.5)) + 
  scale_x_continuous(limits = c(-1.7, 1.25), breaks = c(-1.5, -1, -0.5, 0, 0.5, 1.0), 
                     labels = c(-1.5, -1, -0.5, 0, 0.5, 1.0), 
                     name = bquote(X[ijk])) + 
  scale_y_discrete(name = " ") + 
  theme_minimal() + theme(axis.text.y = element_blank(),
                          axis.title.x = element_text(
                            # family = "italic",
                            size = 12, angle = 0, hjust = 0.54, vjust = 9),
                          legend.position = "none", legend.title = element_blank()
  ) 

trtFig1 = trtp + guides(colour = guide_legend(reverse=TRUE)) + 
  guides(shape = guide_legend(reverse=TRUE)) + labs(subtitle = "(a) Scenario 1") 

trtFig1 

# plot of w_{qk}
HDistmat = array(0, dim = c(nrow(SimSc), nrow(SimSc)))
MyHD = HdMat(mu.hist = SimSc[, 1], var.hist = SimSc[, 2])

for(i in 1:nrow(MyHD)){
  HDistmat[MyHD[i,1], MyHD[i,2]] = MyHD[i,3]
  HDistmat[MyHD[i,2], MyHD[i,1]] = MyHD[i,3]
}


melted_HDistmat <- melt(get_lower_tri(HDistmat))
head(melted_HDistmat)

wiqFig1 = ggplot(melted_HDistmat, aes(x = Var1, y = Var2)) + 
  geom_point(aes(size = value, fill = value), shape = 21, colour = "black") + 
  scale_x_continuous(name = "", breaks = 1:7, labels = paste0("Subtrial ", 1:7)) + 
  scale_y_continuous(name = " ", breaks = 1:7, labels = paste0("Subtrial ", 1:7), 
                     limits = c(0.75, 7.25)) +
  labs(fill = bquote(w[qk])) + 
  scale_fill_gradientn(colours = c(brewer.pal(7, "Set1")[2], "White",
                                   brewer.pal(7, "Set1")[1]), na.value = NA, limits = c(0, 0.62),
                       guide = "colourbar") + 
  scale_size_area(max_size = 5, guide = FALSE) +  ## for Sc1
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.border = element_blank(), 
        # legend.position = "top", 
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

wiqFig1



# Table 1
## Proposed
MyBrwFitSc1 = nlm(MyBrwfun, Ri = 0.5, sig02 = SimSc[,2], s02 = 100, wiq = HDistmat, 
                  cctrpar = 0.05, dw = c(1.1, 1.1), br = c(54, 3), targEff = -0.4, 
                  eta = 0.95, zeta = 0.80,
                  rep(0.0001, nrow(SimSc)))

MyBrwNi1 = MyBrwFitSc1$estimate; MyBrwNi1

## No borrowing
NoBrwFitSc1 = nlm(NoBrwfun, Ri = 0.5, sig02 = SimSc[,2], targEff = -0.4, 
                  eta = 0.95, zeta = 0.80, s02 = 100, 
                  rep(0.0001, nrow(SimSc)))

NoBrwNi1 = NoBrwFitSc1$estimate; NoBrwNi1 
