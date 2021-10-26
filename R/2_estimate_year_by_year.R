# --
# -- Empirical estimation of Linf and k year by year --
# --
rm(list=ls())
library(tidyverse)
library(TMB)

# -- load data --
load("data/official_data.RData")
load("Data/Ndf.RData")

compile("src/vbgf_one_year.cpp")
dyn.load(dynlib("src/vbgf_one_year"))

par.list <- list(
  logLinf = log(41),
  lograte = log(.4),
  t0 = 0,
  logSigma = 0
)
library(tidyverse)
df <- left_join(df, 
                df %>% group_by(julianage, year) %>% count(name = "likelihoodweights"))
df <- left_join(df, 
                df %>% group_by(year) %>% count(name = "totweights"))

data.list <- list(
  y = df$length, 
  age = df$julianage,
  weights = df$likelihoodweights/df$totweights 
  # put "weights = rep(1, nrow(df))" to get unweighted version
)
obj <- MakeADFun(data.list, par.list, DLL = "vbgf_one_year")

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 300))
fixed.pars <- summary(sdreport(obj))[c(5:6, 3, 7),]
par.list <- list(
  logLinf = opt$par[1],
  lograte = opt$par[2],
  t0 = opt$par[3],
  logSigma = opt$par[3]
)

map <-list(
  lograte = factor(NA), 
  t0 = factor(NA), 
  logSigma = factor(NA)
  )

tab.year <- as.data.frame(matrix(ncol = 6+1+1+1+2+1, nrow = length(unique(df$yearclass))))
names(tab.year) <- c(
  "yearclass", "convergence","n", "Linf", "k",
  "sig",
  "seLinf", "seK",
  "seSig", 
  "minage","maxage", "maxlength"
  )
tab.year$yearclass <- unique(df$yearclass)
for(i in 1:nrow(tab.year)){
  data.list <- list(
    y = df$length[df$yearclass == tab.year$yearclass[i]], 
    age = df$julianage[df$yearclass == tab.year$yearclass[i]]#,
  )
  data.list$weights <- rep(1, length(data.list$y))
  obj <- MakeADFun(data.list, par.list, map = map,  DLL = "vbgf_one_year")
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 300))
  sumsd <- summary(sdreport(obj))
  n <- length(data.list$y)
  tab.year[i, -1] <- c(opt$convergence, n, sumsd[c(nrow(sumsd)-c(2,1,0)),],
                       min(data.list$age),
                       max(data.list$age),
                       max(data.list$y))#sumsd[c(5,6,3,7),])
}
if(!all(tab.year$convergence==0))
  tab.year[which(tab.year$convergence==1),]

tab <- left_join(tab.year, N.df, by = "yearclass")
tab$logN <- log(tab$N)
tab$logNsum <- log(tab$Nsum)

# -- Linf as function of time ---
tab.plot <- pivot_longer(tab %>% rename('Asymptotic~length'= Linf, 'log~N[(3)]' = logN), 
             cols = c("Asymptotic~length", "log~N[(3)]"))
tab.plot$lab <- ifelse(tab.plot$name == "log~N[(3)]", "A", "B")
tab.plot$name <- as.factor(tab.plot$name)
tab.plot <- transform(tab.plot, name = factor(name, 
                                              levels = c("log~N[(3)]","Asymptotic~length")))
(p1 <- ggplot(tab.plot, aes(x = yearclass, y = value)) + geom_line() +
  theme_bw() + 
  facet_wrap(~name, ncol = 1, scales = "free_y", 
             strip.position = "left", labeller = label_parsed) +
  geom_text(x = Inf, y = Inf, aes(label = lab), size = 9, hjust = 1.2, vjust = 1.2)+
  scale_x_continuous(name = "Yearclass", breaks = seq(1900,2050,5))+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 12),
        axis.title.y = element_blank()
        ))
ggsave(p1, filename = "plots/2_empirical_Linf_and_logN_by_year.tiff",
       width=28, height=20, units="cm", device = "tiff", dpi = "retina")

# -----------------------------------------------
# -------- Linf vs Log N ------------------------
# -----------------------------------------------
tab2 <- tab %>% pivot_longer(cols = c("logN", "logNsum"), names_to = "model", values_to = "logN") 
tab2 <- mutate(tab2, model = factor(ifelse(model == "logN","log~N[(3)]","log(N[(3)]+N[(4)])"),
                                  levels = c("log~N[(3)]","log(N[(3)]+N[(4)])")))
tab2 <- transform(tab2, model = factor(model,  levels = c("log~N[(3)]","log(N[(3)]+N[(4)])")))
(p2 <- ggplot(tab2, aes(x = logN, y = Linf)) + 
    facet_wrap(~model, ncol = 1, strip.position = "bottom",
               labeller = label_parsed)+
    geom_point(aes(size = n/1000, col = yearclass))+theme_bw() + 
    geom_smooth(method = "lm") + 
    geom_hline(yintercept = fixed.pars[1,1], lty = 2, col = "blue", lwd = .9) +
    scale_x_continuous(name = "", breaks = seq(1,12,1)) + 
    scale_y_continuous(name = "Asymptotic length", #expression(L[infinity]),
                       breaks = seq(30,90,1)) +
    scale_colour_viridis_c(name = "Year class")+
    scale_size(name = "Number of \nsamples \n(in thousands)")+
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "white", color = "white"),
          text = element_text(size = 18))+
    guides(color = guide_colorbar(barheight = 45))
)
ggsave(p2, filename = "plots/4_empirical_Linf_by_logN_and_logNsum.tiff",
       width=28, height=40, units="cm", device = "tiff", dpi = "retina")
