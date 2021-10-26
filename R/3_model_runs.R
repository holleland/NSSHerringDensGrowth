rm(list=ls())
load("data/official_data.Rdata")
load("data/Ndf.RData")

library(TMB)
library(tidyverse)
compile("src/VBGF.cpp")
dyn.load(dynlib("src/VBGF"))


run_model <- function(column_names = c("Nsum", ""),
                                    df, 
                                    N.df, 
                                    fix_t0
                                    ){
#column_names <- c("N", "")
if(column_names[1] == "temp"){
  N1 <- df$temp - mean(df$temp)
}else if(column_names[1] == "lat"){
  N1 <- df$lat - 69
}else if(column_names[1] == ""){
  N1 <- rep(0, nrow(df))
}else{
  N1 <- log(df[,column_names[1]]) -  log(median(unlist(N.df[,column_names[1]])))
}
if(column_names[2] == "temp"){
  N2 <- df$temp - mean(df$temp)
}else if(column_names[2] == "lat"){
  N2 <- df$lat - 69
}else {
  N2 <- rep(0, nrow(df))
}


data <- list(y        = df$length,
             age      = df$julianage, 
             N1       = N1,
             N2       = N2,
             age_indx = df$age - min(df$age),
             NbyYear  = unique(N1)-mean(unique(N1)),
             weight   = rep(1, nrow(df))
)

## Parameter initial guess
parameters <- list(
  Linf        = c(37,  -0.7),
  rate        = c(.4,  0),
  t0          = 0,
  logSigma    = numeric(max(data$age_indx)   +1)
)

# -- mapping --
map <- list()
if(column_names[1] == ""){
  parameters$Linf[2] <- 0 
  map$Linf <- factor(c(1, NA))
}
if(column_names[2] == "")
  map$rate = factor(c(2,NA))
if(fix_t0)
  map$t0 = factor(NA)
map$logSigma =factor(5+1:21)

# -- Fit Gaussian model ---
obj <- MakeADFun(data, parameters, map = map, DLL = "VBGF", silent =TRUE)
timing <- system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr, 
                control = list(eval.max = 500,
                               iter.max = 500))
)
df$pred.norm <- obj$report()$mu
df$gamma.norm <- obj$report()$gamma[df$yearserialno+1]
df$res.norm <- obj$report()$residual
df$SD.norm  <- obj$report()$SD
df$sigma.norm <- exp(opt$par[names(opt$par) == "logSigma"])[as.integer(map$logSigma)[df$age]]
#df$sigma.norm <- exp(opt$par[names(opt$par) == "logSigma"])[as.integer(map$logSigma)[df$age]]
df$stdres <- df$res.norm/df$sigma.norm
rep <- sdreport(obj)

# Scenario description
scenario = paste0("L[infinity]~depends~on~", ifelse(column_names[1]=="temp", "KOLA~temperature",
                                             ifelse(column_names[1]=="N", "N[(3)]", "N[(3)]+N[(4)]")))
if(column_names[2] != "") scenario = paste0(scenario,"~and~the~rate~on~",
                                            ifelse(column_names[2]=="temp", "KOLA~temperature~",
                                                   ifelse(column_names[2]=="lat", "latitude~", column_names[2])))
if(fix_t0) scenario <- paste0(scenario, "-~t[0]~is~fixed")

model <- list(
  obj = obj, 
  opt = opt, 
  rep = rep,
  df = df,
  column_names = column_names,
  fix_t0 = fix_t0,
  scenario = scenario,
  AIC = TMBhelper::TMBAIC(opt = opt),
  timing = timing
)
return(model)
}
# ______________________________________________
#
# ----------------------------------------------
# -- FUNCTION FOR GENERATING PREDICTION PLOTS --
# ----------------------------------------------
# ______________________________________________
prediction_plot <- function(model, two, ind = 1, model1=NULL){
  theme_set(theme_bw())
  if(!is.null(model1)){
    df.model0 <- tibble(
      age = seq(0,21,0.2),
      pred = model1$opt$par[1]*(1-exp(-model1$opt$par[2]*(age-model1$opt$par[3])))
    )
  }
  if(model$column_names[2]!="") {
    
    # -----------------------------------------------------------
    # --------------------- Predictions -------------------------
    # -----------------------------------------------------------
    
    # .. calculate predictions..
    mean.pred <- function(x, par, means,ind=1){
      .age = x[1]
      .N1 = x[2]-means[1]
      .N2 = x[3] - means[2]
      (par[1] + par[2] * .N1) * (1-exp(-(par[3]+par[4]*.N2) * (.age - par[5]) ))
    }
    age <- seq(0, 21.5, by = 1/12)
    if(model$column_names[1] %in% c("N","Nsum", "TN")){
      med <- log(median(unlist(N.df[, model$column_names[1]])))
      N1 <- seq(min(log(model$df[, model$column_names[1]])), 
                  max(log(model$df[, model$column_names[1]])), 
                  length.out = 200)
    }else if(model$column_names[1] == "temp"){
      med = mean(model$df$temp)
      N1 = seq(min(model$df$temp), max(model$df$temp), length.out = 200)
    }else {N1 <- 0}
    if(model$column_names[2] == "temp"){
      N2 = seq(3,5.5,length.out = 4)
      med2 = mean(model$df$temp)
      #N2 = seq(min(model$df$temp), max(model$df$temp), length.out = 200)
    }else if(model$column_names[2]== "lat"){
      N2 = seq(58,76,length.out = 4)
      med2 = 69
      # N2 = seq(min(model$df$lat), max(model$df$lat), length.out = 200)
    }
    #tmp = seq(3,5.5,length.out = 4)
    pred.df <- expand.grid(age, N1, N2)
    names(pred.df) <- c("age", "N1", "N2")
    pred.df$prediction<- apply(pred.df[,1:3], 1, mean.pred, par = model$opt$par, means = c(med, med2))
    if(model$column_names[2] == "temp"){
      pred.df$N2fac <- factor(pred.df$N2, levels = seq(3,5.5,length.out = 4), 
                                labels = paste0(round(seq(3,5.5,length.out = 4),2), "°C"))
    }else{
      pred.df$N2fac <- factor(pred.df$N2, levels = seq(58,76,length.out = 4), 
                                labels = paste0(round(seq(58,76,length.out = 4),2), "°N"))
    }
    # .. figure ..
    pred1 <- ggplot(pred.df, aes(x = age, y = prediction, col = N1, group = N1))+ 
      geom_line()+
      scale_x_continuous(name = "Age", breaks = seq(0,26,2),
                         expand =c(0,0), minor_breaks = seq(0,26,1))+
      scale_y_continuous(name = "Predicted length", breaks = seq(0,50,5), 
                         minor_breaks = seq(0,50,5/2),
                         limits = c(0,45))+
      guides(color = guide_colorbar(barheight = 1, barwidth = 30)) + 
      facet_wrap( ~N2fac, nrow = 1)+
      theme(legend.position = "top",
            legend.title = element_text(vjust = .9, size = 12), 
            axis.text = element_text(size = 11), 
            axis.title = element_text(size = 14),
            strip.background = element_rect(fill = "white", colour = "white"),
            strip.text = element_text(size = 12))
    if(model$column_names[1] %in% c("N", "Nsum")){
      pred1 <- pred1 + scale_color_viridis_c(name = "Logarithmic density", 
                                             breaks =seq(1,11,1))
    }else{
      pred1 <- pred1 + scale_color_viridis_c(name = "Kola Temperature", 
                                             breaks =seq(1,11,.5),
                                             labels = paste0(seq(1,11,.5), "°C"))
    }
    if(!is.null(model1)){
      pred1 <- pred1 + geom_line(data = df.model0, aes(x= age, y = pred, group = NULL), col = "black", lwd = 1, 
                                 lty = 2)
    }
    #  pred1
    # ggsave(pred1, filename = "plots/5_temp_and_csize_prediction.tiff",
    #        width=25*1.5, height=10*1.5, units="cm", device = "tiff", dpi = "retina")
    # # -------------------------------------
    # -------------------------------------
    N1 = round(seq(min(N1), max(N1), length.out = 4), 2)
    N2 <- seq(min(N2), max(N2), length.out = 200)
   
    pred.df <- expand.grid(age, N1, N2)
    names(pred.df) <- c("age", "N1", "N2")
    pred.df$prediction<- apply(pred.df[,1:3], 1, mean.pred, par = model$opt$par, means = c(med, med2))
    pred.df$Linf <- (model$opt$par[1] + model$opt$par[2] * (pred.df$N1-med))
    if(model$column_names[1] == "temp"){
      pred.df$N1fac <- factor(paste0(pred.df$N1, "°C"), levels =  paste0(N1, "°C"))
    }else{
      pred.df$N1fac <- factor(paste0("Logarithmic density = ",pred.df$N1), levels = paste0("Logarithmic density = ",N1))
    }
    #pred.df$logN <- factor(pred.df$logN, levels = logN, labels = paste0("Logarithmic cohort size: ",round(logN,2)))
    # .. figure ..
    pred2 <- ggplot(pred.df, aes(x = age, y = prediction, col = N2, group = N2))+ geom_line()+
      #scale_color_viridis_c(name = "Temperature (°C)", breaks = seq(3,5.5,.5))+
      scale_x_continuous(name = "Age", breaks = seq(0,26,2),expand =c(0,0), minor_breaks = seq(0,26,1))+
      scale_y_continuous(name = "Predicted length", breaks = seq(0,50,5), minor_breaks = seq(0,50,5/2),
                         limits = c(0,45))+
      guides(color = guide_colorbar(barheight = 1, barwidth = 30)) + 
      geom_hline(aes(yintercept = Linf), lty =2, col = "darkblue")+
      facet_wrap( ~N1fac, nrow = 1)+
      theme(legend.position = "top",
            legend.title = element_text(vjust = .9, size = 12), 
            axis.text = element_text(size = 11), 
            axis.title = element_text(size = 14),
            strip.background = element_rect(fill = "white", colour = "white"),
            strip.text = element_text(size = 12))
    if(model$column_names[2]=="temp"){
      pred2 <- pred2 + scale_color_viridis_c(name = "Temperature (°C)", breaks = seq(3,5.5,.5))
    } else{
      pred2 <- pred2 + scale_color_viridis_c(name = "Latitude", breaks = seq(min(N2),max(N2),2),
                                             labels = paste0(seq(min(N2),max(N2),2), "°N"))
    }
    if(!is.null(model1)){
      pred2 <- pred2 + geom_line(data = df.model0, aes(x= age, y = pred, group = NULL), col = "black", lwd = 1, 
                                 lty = 2)
    }
    library(ggpubr)
    
    pred.gaus <- ggarrange(pred1,pred2, ncol = 1)
    
  }else{
  theme_set(theme_bw())
  
  mean.pred <- function(x, par){
    .age = x[1]
    if(model$column_names[1]!=""){
      .logN = x[2]-med
    }else{ .logN = 0}
    return((par[1] + par[2] * .logN) * (1-exp(-(par[3]) * (.age - par[4]) )))
  }
  age <- seq(0, 21.5, by = 1/12)
  if(model$column_names[1] %in% c("N","Nsum", "TN")){
    med <- log(median(unlist(N.df[, model$column_names[1]])))
    logN <- seq(min(log(model$df[, model$column_names[1]])), 
                max(log(model$df[, model$column_names[1]])), 
                length.out = 200)
  }else if(model$column_names[1] == "lat"){
    med = mean(df$lat)
    logN <- seq(min(model$df$lat), max(model$df$lat), length.out = 200)
  }else if(model$column_names[1] == "temp"){
    med = mean(df$temp)
    logN = seq(min(model$df$temp), max(model$df$temp), length.out = 200)
  }else {logN <- 0}
  pred.df <- expand.grid(age, logN)
  names(pred.df) <- c("age", "logN")
  if(!model$fix_t0 & model$column_names[1]==""){
    predpar <- c(model$opt$par[1], 0, model$opt$par[2:3])
  }else if(!model$fix_t0 & model$column_names[1]!=""){
    predpar <- model$opt$par[1:4]
  }else if(model$fix_t0 & model$column_names[1]==""){
    predpar <- c(model$opt$par[1], 0, model$opt$par[2],0)
  }else{ 
    predpar <- c(model$opt$par[1:3],0)
  }
  pred.df$prediction <- unlist(apply(pred.df[,1:2], 1, 
                                          mean.pred, 
                                          par =predpar))
  if(model$column_names[1] == "N"){
    dfmeans <- model$df %>% mutate(N = log(N)) %>% group_by(N, julianage) %>%
      summarize(mlength = mean(length)) %>% ungroup
  } else if(model$column_names[1] == "Nsum"){
    dfmeans <- model$df %>% mutate(Nsum = log(Nsum)) %>% group_by(Nsum, julianage) %>%
      summarize(mlength = mean(length)) %>% ungroup
  } else if(model$column_names[1] == "TN"){
    dfmeans <- model$df %>% mutate(TN = log(TN)) %>% group_by(TN, julianage) %>%
      summarize(mlength = mean(length)) %>% ungroup
  } else if(model$column_names[1] == "temp"){
    dfmeans <- model$df %>% group_by(temp, julianage) %>%
      summarize(mlength = mean(length)) %>% ungroup
  } else if(model$column_names[1] == "lat"){
    dfmeans <- model$df %>% group_by(lat, julianage) %>%
      summarize(mlength = mean(length)) %>% ungroup
  } else{
    dfmeans <- model$df %>% group_by(julianage) %>%
      summarize(mlength = mean(length)) %>% ungroup
  }

  if(model$column_names[1]!=""){
    dfmeans$logN <- log(unlist(dfmeans[,model$column_names[1]]))
  scalename <- ifelse(model$column_names[1] %in% c("N", "Nsum", "TN"), "Logarithmic\ndensity", 
                  ifelse(model$column_names[1] == "temp", "Kola\ntemperature\n(°C)",
                         model$column_names[1] ))
  
  pred.gaus <- ggplot(pred.df, aes_string(x = "age", y = "prediction", 
                                          col = "logN", 
                                          group = "logN"
  ))+ 
    geom_line()+
    scale_color_viridis_c(name = scalename)+#, breaks = seq(1,11,1))+
    scale_x_continuous(name = "Age", breaks = seq(0,26,2),expand =c(0,0), minor_breaks = seq(0,26,1))+
    scale_y_continuous(name = "Predicted length", breaks = seq(0,70,2.5), minor_breaks = seq(0,70,2.5/2),
                       limits = c(0,45))+
    guides(color = guide_colorbar(barheight = 20, barwidth = 1)) +
    theme(legend.position = "right",
          legend.title = element_text(vjust = .9, size = 12), 
          axis.text = element_text(size = 12), 
          axis.title = element_text(face = "bold", size = 14),
          title = element_text(size = 15, face = "bold"))+
    geom_text(label = paste0("AIC: ",
                             format(TMBhelper::TMBAIC(model$opt),
                                    digits = 2, 
                                    big.mark = " ")),
              x = Inf, y = -Inf, col = "black",#fill = "skyblue",
              vjust = -0.5, hjust = 1)+
    ggtitle(paste0("Scenario ", ind), parse(text = model$scenario)) + 
    geom_point(data = dfmeans, 
      aes_string(x="julianage", y = "mlength", 
      col = model$column_names[1]), 
      size = .5)

  }else{
    pred.gaus <- ggplot(pred.df, aes_string(x = "age", y = "prediction"))+ 
      geom_line()+
      #scale_color_viridis_c(name = "Logarithmic cohort size", breaks = seq(1,11,1))+
      scale_x_continuous(name = "Age", breaks = seq(0,26,2),expand =c(0,0), minor_breaks = seq(0,26,1))+
      scale_y_continuous(name = "Predicted length", breaks = seq(0,70,2.5), minor_breaks = seq(0,70,2.5/2),
                         limits = c(0,45))+
      guides(color = guide_colorbar(barheight = 1, barwidth = 30)) +
      theme(legend.position = "top",
            legend.title = element_text(vjust = .9, size = 12), 
            axis.text = element_text(size = 12), 
            axis.title = element_text(face = "bold", size = 14))+
      geom_text(label = paste0("AIC: ",
                               format(TMBhelper::TMBAIC(model$opt),
                                      digits = 2, 
                                      big.mark = " ")),
                x = Inf, y = -Inf, col = "black",#fill = "skyblue",
                vjust = -0.5, hjust = 1)+
      ggtitle("No density effect") + geom_point(data = dfmeans, 
                                           aes_string(x="julianage", y = "mlength"), 
                                           size = .5)
  }
  if(!is.null(model1)){
    pred.gaus <- pred.gaus + geom_line(data = df.model0, 
                                       mapping = aes_string(x= "age", y = "pred", group = NULL), 
                                       col = "black", lwd = 1, 
                                       lty = 2)
  }
  }
  #print(pred.gaus)
  beepr::beep(2)
  return(pred.gaus)
}
# ______________________________________________
#
# ----------------------------------------------
# ----- FUNCTION FOR GENERATING QQ PLOTS -------
# ----------------------------------------------
# ______________________________________________
QQplots <- function(model, ind, save = FALSE){
  theme_set(theme_bw())
  dens <- "norm"
  p1 <- ggplot(model$df, aes_string(x = "julianage", y = "stdres", 
                                    group = "length"))+ 
    geom_hline(yintercept = 0, col = "black")+
    geom_point(size = .4, col = "skyblue") + 
    scale_x_continuous(name = "Age", breaks = seq(2,20,2), minor_breaks = seq(1,21,1))+
    scale_y_continuous(name = "Standardized residuals", breaks = seq(-20,20,1))
  p2 <- ggplot(model$df, aes_string(x = "length", y = "stdres"))+
    geom_hline(yintercept = 0, col = "black")+
    geom_point(size = .4, col = "skyblue") + 
    scale_x_continuous(name = "Length", breaks = seq(-5,50,5))+
    scale_y_continuous(name = "Standardized residuals", breaks = seq(-20,20,1))+
    geom_smooth(method = "gam")
  p3 <- ggplot(model$df, aes_string(x = "length", y = "length-res.norm"))+
    geom_abline(intercept = 0,slope = 1, col = "black")+
    geom_point(size = .4, col = "skyblue") + 
    scale_x_continuous(name = "Length", breaks = seq(-5,50,5))+
    scale_y_continuous(name = "Fitted values", breaks = seq(-5,50,5))+
    geom_smooth(method = "gam")
  #g1 <- ggpubr::ggarrange(p1,p2,p3, ncol = 3)
  
  # .. by age ..
  agelabs <- paste0("Age ",1:21)
  names(agelabs) <- 1:21
  model$df$agefac <- factor(model$df$age)
  
  qq1<-ggplot(model$df, aes_string(sample = "stdres"))+
    geom_text(aes(label = paste0("Age ", agefac)), x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5,
              fontface = "bold") +
    facet_wrap( ~agefac, ncol= 7,
                labeller = labeller(agefac=agelabs))+
    stat_qq(cex = .1) + 
    geom_abline(intercept = 0, slope = 1, lwd = .8, col = "red", lty = 2)+
    scale_y_continuous(name = "Residual quantiles",
                       breaks = seq(-10,10,2),
                       minor_breaks = seq(-10,10,1))+
    scale_x_continuous(name = "Gaussian quantiles", 
                       breaks = seq(-10,10,2),
                       minor_breaks = seq(-10,10,1))+
    
    theme(strip.text =element_blank(),
          strip.background = element_blank(),
          title = element_text(size = 20)
          #axis.text = element_text(size = 12),
          # axis.title = element_text(size = 14, face = "bold")
    )+
    ggtitle(" ")#paste0("Scenario ",ind),
    #        model$scenario)
  
  # .. all in one ..
  qq2 <- ggplot(model$df, aes_string(sample = "stdres"))+
    stat_qq(cex = .1) + 
    geom_abline(intercept = 0, slope = 1, lwd = .8, col = "red", lty = 2)+
    scale_y_continuous(name = "Residual quantiles",
                       breaks = seq(-10,10,1),
                       minor_breaks = seq(-10.5,10.5,1))+
    scale_x_continuous(name = "Gaussian quantiles", 
                       breaks = seq(-10,10,1),
                       minor_breaks = seq(-10.5,10.5,1)) +
    geom_text(label = paste0("AIC: ",
                             format(TMBhelper::TMBAIC(model$opt),digits = 2, big.mark = " ")), 
              x = -Inf, y = Inf,
              color = "black", vjust = 1.2, hjust = -0.1)
  
  g1 <- ggpubr::ggarrange(qq2, p1,p2,p3, ncol = 4, labels = c("B","C","D","E"))
  qq4 <- ggpubr::ggarrange(qq1, 
                           g1, ncol = 1, labels = c("A",NA), 
                           heights = c(2, .8))
  if(save){
    ggsave(qq4, filename = paste0("plots/", 
                                  str_remove_all(str_replace_all(model$scenario," ","_"),","),
                                  ".tiff"),
           width=32, height=25, units="cm", device = "tiff", dpi = "retina")
  }
  return(qq4)
}
# ______________________________________________
#
# ----------------------------------------------
# ----------------------------------------------
# ______________________________________________

scenarios = rbind(
  expand.grid(
    N1 = c("","N", "Nsum"),#,"temp"),
    N2 = c(""),
    fix_t0 = c(FALSE)
  ),
  expand.grid(
    N1 = c("N","Nsum"),
    N2 = "",
    fix_t0 = TRUE
  ),
  expand.grid(
    N1 = c("Nsum"),#, "temp"),
    N2 = c("lat"),
    fix_t0 = c(FALSE)
  ),
  expand.grid(
    N1 = c("Nsum"),
    N2 = c("temp"),
    fix_t0 = c(FALSE)
  ))
scenarios$N1 <- as.character(scenarios$N1)
scenarios$N2 <- as.character(scenarios$N2)
model.list <-list()
t1 <- Sys.time()
for(i in 1:nrow(scenarios)){
  model.list[[i]] <- run_model(column_names = unlist(scenarios[i,c("N1","N2")]),
                                             df = df, 
                                             N.df = N.df,
                                             fix_t0 = scenarios$fix_t0[i])
  cat("-------------------------------",
      "\n--- Scenario ", format(i,width = 3), " of ", 
      format(nrow(scenarios),width =3)," ---",
      "\n-------------------------------\n")
  beepr::beep(2)
}
beepr::beep(4)
(t2 <- Sys.time())
t2-t1

cat("Did all converge?", ifelse(all(unlist(
  lapply(model.list, function(x) x$opt$convergence))==0), 
  "\n- yes", "\n- no"))
save(model.list, file = "data/model_runs.RData")

# Checking AIC
t(t(unlist(lapply(model.list, function(x){y <- x$AIC; names(y) <- x$scenario;y}))))

pred.list <- list()
for(i in c(2:nrow(scenarios))){
  pred.list[[i]] <- prediction_plot(model.list[[i]], ind = i, model1 = model.list[[1]]) 
}

library(ggpubr)
pred.panels<-ggarrange(pred.list[[2]], pred.list[[3]], pred.list[[4]], pred.list[[5]],
          ncol = 2, nrow = 2)
ggsave(plot = pred.panels, filename= "plots/S6_one_effect_predictions.tiff", device ="tiff", 
       dpi = "retina", width = 16, height = 12)
pred.panels2 <- ggarrange(pred.list[[6]],
                          pred.list[[7]], nrow = 2,ncol = 1,
                          labels = paste0("Scenario ", 6:7),
                          font.label = list(size = 20))

ggsave(plot = pred.panels2, filename= "plots/S7_two_effect_predictions.tiff", device ="tiff", 
       dpi = "retina", width = 15, height =20)

# ------------------------
# Additional filtration to exclude potential other herring species
# ------------------------
withAdditonalFiltration.model <-   run_model(
  column_names = c("Nsum", ""),
  df = filter(df, 
              !(maturity %in% 4:8 & length <27)) %>%
    filter( !(length < 28 & age >= 6),
            !(length < 29 & age >= 7),
            !(length <30  & age >= 8)),
  N.df=N.df,fix_t0 = FALSE
)
withAdditonalFiltration.model $scenario <- "With filtration"

p1 <- QQplots(model.list[[3]], save = FALSE) # plot the main model QQplots
p2 <- QQplots(withAdditonalFiltration.model,  save = FALSE)    # plot the filtration version

ggsave(p1, filename = "plots/S2.tiff", width=32, height=25, units="cm", device = "tiff", dpi = "retina" )
ggsave(p2, filename = "plots/S3.tiff", width=32, height=25, units="cm", device = "tiff", dpi = "retina" )

# ----------------------------------------
# ---- The code below can be used to look at residual plots for every scenario --- 
# ---- (commented out by default)
# ----------------------------------------
# QQ.list <- list()
# for(i in 1:nrow(scenarios)){
#   QQ.list[[i]] <- QQplots(model.list[[i]], ind = i)
# }
# for(i in 1:nrow(scenarios)){
#   QQ.panels <- ggarrange(QQ.list[[i]], labels = paste0("Scenario ", i),
#                          font.label = list(size = 20))
#   ggsave(plot = QQ.panels, filename= paste0("plots/supp_residual_plots_scenario_",i,".tiff"),
#          device ="tiff",
#          dpi = "print", width=16, height=10)
#   cat(i,"\n")
# }
# beepr::beep(4)
