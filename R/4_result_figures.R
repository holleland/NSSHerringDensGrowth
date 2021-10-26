rm(list=ls())
library(tidyverse)
library(TMB)
library(ggpubr)
library(ggExtra)
theme_set(theme_bw())
load("data/model_runs.RData")
load("data/Ndf.RData")

# The chosen model: 
model <- model.list[[3]]
df <- model$df
# ----------------------------------------------------
# -- Plot relationship between Linf and cohort size --
# ----------------------------------------------------

med <- median(N.df$Nsum)
a <- model$opt$par[1]
b <- model$opt$par[2]
df <- model$df
N.df <- filter(N.df, yearclass %in% unique(df$yearclass))
#N.df$Linf <- a + b*(log(N.df$N)-log(med))
srep <- summary(model$rep)


p1<-ggplot(N.df, aes(x = log(Nsum)))+
  #geom_rug()+
  geom_point(alpha = 0, aes(y = a))+
  geom_line(aes(y = a + b*(log(Nsum)-log(med))), lwd = .9, col = "blue")+
  geom_segment(y = a,yend=a,x=-Inf, xend = log(med), col = "blue", lty = 2) + 
  geom_segment(x = log(med),xend=log(med), y = -Inf, yend = a, lty = 2, col = "blue") +
  scale_y_continuous(name = "Asymptotic length", #expression(l[infinity]), 
                     breaks = c(seq(34,50.5,.5)),
                     labels = c(sprintf("%.1f",seq(34,50.5,.5))),
                     minor_breaks = seq(34.5,50.5,.25))+
  # geom_text(aes( x = -Inf, y = a),label = "l[0]^infinity", parse = TRUE, 
  #           hjust =-.5, vjust = -0.1,
  #           size = 5,
  #           col = "blue")+
  scale_x_continuous(name = expression(log~N[(3)]+N[(4)]),  breaks = c(2:11),
                     labels = c(2:11),
                     minor_breaks = NULL) +
  theme(axis.title.y = element_text( size = 12),
        axis.title.x = element_text(size = 12), 
        axis.text = element_text(size = 12))
(p2 <- ggExtra::ggMarginal(p1, type=c("densigram"), margins = "x", col = "blue", fill = "skyblue"))
#ggarrange(p1,p2, ncol = 1)
ggsave(p2, filename = "plots/6_Linf_vs_logN_Nsum.tiff",
       width=28, height=20, units="cm", device = "tiff", dpi = "retina")

# -----------------------------------------------------------
# --------------------- Predictions -------------------------
# -----------------------------------------------------------

# .. calculate predictions..
mean.pred <- function(x, par){
  .age = x[1]
  .logN = x[2]-log(med)
  (par[1] + par[2] * .logN) * (1-exp(-(par[3]) * (.age - par[4]) ))
}
age <- seq(0, 21.5, by = 1/12)
logN <- seq(min(log(df$Nsum)), max(log(df$Nsum)), length.out = 200)
pred.df <- expand.grid(age, logN)
names(pred.df) <- c("age", "logN")
pred.df$prediction<- apply(pred.df[,1:2], 1, mean.pred, par = model$opt$par)

df$logNfac <- factor(
  ifelse(log(df$N)<log(med),paste0("Low density (N < ",round(med,1),")"),"High density")
)
# .. figure ..
dfmeans <- df %>% group_by(Nsum, julianage) %>%
  summarize(mlength = mean(length)) %>%mutate(logN = log(Nsum))
predmean <- ggplot(pred.df, aes(x = age, y = prediction, col = logN, group = logN))+ geom_line()+
  scale_color_viridis_c(name = "Logarithmic density", breaks = seq(1,11,1))+
  scale_x_continuous(name = "Age", breaks = seq(0,26,2),expand =c(0,0), minor_breaks = seq(0,26,1))+
  scale_y_continuous(name = "Predicted length", breaks = seq(0,50,2.5), minor_breaks = seq(0,50,2.5/2))+
  guides(color = guide_colorbar(barheight = 1, barwidth = 30)) + 
  theme(legend.position = "top",
        legend.title = element_text(vjust = .9, size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(face = "bold", size = 14))+
  geom_point(data = dfmeans, aes(x=julianage, y = mlength, col = logN), size = .5)

dfmeans <- df %>% group_by(Nsum, julianage) %>%
  summarize(mlength = median(length)) %>%mutate(logN = log(Nsum))
predmedian <- ggplot(pred.df, aes(x = age, y = prediction, col = logN, group = logN))+ geom_line()+
  scale_color_viridis_c(name = "Logarithmic density", breaks = seq(1,11,1))+
  scale_x_continuous(name = "Age", breaks = seq(0,26,2),expand =c(0,0), minor_breaks = seq(0,26,1))+
  scale_y_continuous(name = "Predicted length", breaks = seq(0,50,2.5), minor_breaks = seq(0,50,2.5/2))+
  guides(color = guide_colorbar(barheight = 1, barwidth = 30)) + 
  theme(legend.position = "top",
        legend.title = element_text(vjust = .9, size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(face = "bold", size = 14))+
  geom_point(data = dfmeans, aes(x=julianage, y = mlength, col = logN), size = .5)

ggsave(predmean, filename = "plots/5a_prediction_with_mean_lengths_julianage_Nsum.tiff",
       width=25, height=18, units="cm", device = "tiff", dpi = "retina")
ggsave(predmedian, filename = "plots/5b_prediction_with_median_lengths_julianage_Nsum.tiff",
       width=25, height=18, units="cm", device = "tiff", dpi = "retina")

# -----------------------------------------------------------
# -------- Residuals vs other explanatory variables ---------
# -----------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(1)[1]

# .. gear code ..
df <- transform(df, gearcode = factor(gearcode, levels = c("Trawl", "Purse Seine","Other")))
gear1 <- ggplot(df, aes(x= julianage, y = res.norm/sigma.norm, col = gearcode)) + 
  geom_point(cex = .6) +
  facet_wrap(~gearcode, ncol = 1, strip.position = "right") + 
  scale_x_continuous(breaks = seq(0,22,2), name = "Age") + 
  scale_y_continuous(breaks =seq(-20,20,5), name = "Residuals") +
  theme(strip.background = element_rect(fill = "white"), 
        legend.position = "none")
gear2 <- ggplot(df, aes(x= gearcode, y = res.norm/sigma.norm, col = gearcode)) + 
    geom_boxplot(outlier.size = 0.1)+
    scale_y_continuous(breaks =seq(-20,20,5), name = "Residuals") + xlab("")+
    theme(legend.position = "none")
gear3 <- ggplot(df, aes(sample = res.norm/sigma.norm, col = gearcode))+
  facet_wrap( ~factor(gearcode), ncol = 1, strip.position = "right")+
  stat_qq(cex = .1) +
  geom_abline(intercept = 0, slope = 1, lwd = .5, col = "black", lty = 3)+
  scale_y_continuous(breaks =seq(-100,100,2), name = "Residual quantiles") +
  scale_x_continuous(breaks = seq(-20,20,2), name = "Theoretical Gaussian quantiles") +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none")
library(ggpubr)
gear4 <- ggarrange(gear1,ggarrange(gear2,gear3, labels = c("B","C")), ncol = 1, labels = c("A",NA))
ggsave(gear4, filename = "plots/S9_gearcode_Nsum.tiff",
       width=25, height=25, units="cm", device = "tiff", dpi = "retina")

# ---------------------------------
# -------------------------------------------------------
# ------ Summer temperature against model residuals -----
# -------------------------------------------------------
df <- model.list[[3]]$df
df$agefac <- factor(paste0("Age ", df$age))
df <- transform(df, agefac = factor(agefac, levels = paste0("Age ",1:21)))
# ..summer..
temp1 <- ggplot(df, aes(x = temp, y = res.norm/sigma.norm)) + geom_point(cex = .3)+
  facet_wrap(~agefac, ncol = 3, strip.position = "right") + 
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "", breaks = seq(-6,6,2), minor_breaks = seq(-6,6,1)) +
  theme(strip.background = element_rect(fill = "white"))
temp2 <- ggplot(df, aes(x = temp, y = res.norm/sigma.norm)) + geom_point(cex = .5)+
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "", breaks = seq(-6,6,2), minor_breaks = seq(-6,6,1))
ggtemp <- ggarrange(temp2,temp1, ncol = 2)
ggtemp <- annotate_figure(ggtemp,
                          #top = text_grob("Visualizing len", color = "red", face = "bold", size = 14),
                          bottom = text_grob("Mean summer temperature (0-min(age,3))", color = "black",
                                             vjust = -.7),
                          left = text_grob("Gaussian model residuals", color = "black", rot = 90,
                                           vjust = 2, hjust = 0.5)
                          #right = "I'm done, thanks :-)!",
                          #fig.lab = "Figure 1", 
                          #fig.lab.face = "bold"
)
ggsave(ggtemp,
       filename = "plots/7_temp_vs_Nsum_residuals.tiff",
       width=32, height=18, units="cm", device = "tiff", dpi = "retina")