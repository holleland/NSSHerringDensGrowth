rm(list=ls())
library(tidyverse)
library(TMB)
load("data/model_runs.RData")
model <- model.list[[3]]

# -----------------------------------
# ----------- Table 1 ---------------
# -----------------------------------

ss <- summary(model$rep)
esttab <- ss[(rownames(ss) %in% c("Linf", "rate", "t0", "sigma")),]
#esttab <- head(summary(models$t$rep), n = 11)

esttab <- data.frame(esttab)
names(esttab) <- c("Estimate", "Std.Err")
esttab$`Z-score` <- esttab$Estimate/esttab$Std.Err
esttab$`p-value` <- 2*pnorm(abs(esttab$`Z-score`), lower.tail=FALSE)
esttab$` ` <- gtools::stars.pval(esttab$`p-value`)
esttab$Parameter <- rownames(esttab)
esttab <- esttab[, c(ncol(esttab), 2:ncol(esttab)-1)]

write_excel_csv(esttab, file = "tables/table1.csv", delim= ";")

# -----------------------------------
# ---------- Table S1 ---------------
# -----------------------------------
esttab <- matrix(ncol = length(model.list), nrow = 6)
esttab[1,] <- unlist(lapply(model.list, function(x){y <- x$AIC; names(y) <- x$scenario;y}))
esttab[1,] <- esttab[1,]-esttab[1,3]
esttab[-1,] <- matrix(unlist(lapply(model.list, function(x) {
  if(x$fix_t0){
    return(c(x$opt$par[1:3],NA,0))
  }else if(all(x$column_names == "")){
    return(c(x$opt$par[1],NA, x$opt$par[2], NA, x$opt$par[3]))
  }else if(all(x$column_names[1] %in% c("N","Nsum", "TN", "temp")) & x$column_names[2]==""){
    return(c(x$opt$par[1:3],NA, x$opt$par[4]))
  }else {
    return(c(x$opt$par[1:5]))
  }
}
) ), nrow = 5, byrow=F)
esttab
rownames(esttab) <- c("delta_AIC", "Linf0", "Linf1", "K0", "K1", "t0")
colnames(esttab) <- paste0("Scenario ", 1:length(model.list))
write.table(as.data.frame(esttab), file = "tables/tableS1.csv", dec= ".", sep = ";")


# .. N at age and year ..
N.df <- read.table("inputdata/N.txt",
                   header = TRUE)
names(N.df) <- str_remove(names(N.df), "X")
t(N.df)
N.df <- as.data.frame(t(N.df))
names(N.df) <- paste0("Age ", names(N.df))
N.df <- N.df %>% rownames_to_column("Year")
write_excel_csv(N.df, file = "tables/XSAM_Series.csv")
