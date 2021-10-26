# Estimating relationship between length and weight
rm(list=ls())
library(tidyverse)
library(TMB)
# -- load data --
load("data/official_data.RData")
load("data/Ndf.RData")
load("data/model_runs.RData")
#load("Data/prediction_dataframe.rda")
#source("R/estimate_vbgf.R")
model <- model.list[[3]]
theme_set(theme_bw())
med <- median(unlist(N.df[, "N"]))
mean.pred <- function(x, par, NonK = FALSE){
  .age = x[1]
  .logN = x[2]-log(med)
  if(!NonK){
    return((par[1] + par[2] * .logN) * (1-exp(-(par[3]) * (.age - par[4]) )))
  }else{
    return((par[1] + par[2] * .logN) * (1-exp(-(par[3] + par[4]*.logN) * (.age - par[5]) )))
  }
}

pred.df <- data.frame(length = seq(0,40,0.1))

# lw.mod <- glm(log(weight)~log(length), data = filter(df), family = gaussian(link = "identity"))
lw.mod <- glm(weight~log(length), data = filter(df, month == 1& !is.na(weight)), family = gaussian(link = "log"))

pred.df$weight <- predict(lw.mod, newdata= pred.df, type = "response")

p1 <- ggplot(filter(df,!is.na(weight), month == 1), aes(x = log(length), y = log(weight))) + geom_point() +
  geom_abline(intercept = lw.mod$coefficients[1], slope = lw.mod$coefficients[2], col = "blue", lwd = .8) + 
  scale_x_continuous(name = "Logarithmic length", breaks = seq(2,4,.25))+
  scale_y_continuous(name = "Logarithmic weight", breaks = seq(0,10,.5))
p2 <- ggplot(filter(df,!is.na(weight), month ==1), aes(x = length, y= weight)) +
  geom_point()+
  geom_line(data= pred.df, aes(x = length, y = weight),col = "blue", lwd = .8) +
  scale_x_continuous(name = "Length (cm)", breaks = seq(0,100,10)) + 
  scale_y_continuous(name = "Weight (g)", breaks = seq(0,1000,100))

pred2 <- ggpubr::ggarrange(p1,p2, nrow = 1, ncol = 2, labels = c("A","B"))
ggsave(pred2, file = "plots/S10_length_vs_weight.tiff", device = "tiff", dpi = "retina", width = 12, height = 6)

cat("Length-vs-weight model")
print(summary(lw.mod))
par(mfrow=c(2,2))
plot(lw.mod)
