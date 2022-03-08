# Script for structuring the data. 
# --------------------------------------------------
# -------------- LOADING PACKAGES ------------------
# --------------------------------------------------
# 
# The script will check if inputdata/kola.RData is available and if not
# it will put NAs for the df$temp column. 
# The script will assume that 
#     inputdata/N.txt - the XSAM series of number at age 
# is available. We are currently not at liberty to publish this series,
# but refer the user to Supplementary Table S3 in the paper. You can 
# download the table from there and save it as inputdata/N.txt.
# The herring individual data is downloaded automatically from NMDC
# if the file inputdata/HerringData.csv is not available. 

library(tidyverse)
library(sf)

# Download data file from 
# --------------
# Erling Kåre Stenevik (HI), Sondre Hølleland (HI), Katja Enberg (UiB), 
# Åge Høines (HI), Are Salthaug (HI), Aril Slotte (HI), Sindre Vathehol (HI), 
# Sondre Aanes (NR) (2022) Individual samples of Norwegian Spring Spawning 
# herring 1935-2019 https://doi.org/10.21335/NMDC-496562593
# ------------
if(!("HerringData.csv" %in% list.files(path = "inputdata/") )) {
  download.file(url = "https://ftp.nmdc.no/nmdc/IMR/Herring/HerringData.csv", 
                destfile = "inputdata/HerringData.csv")
}

# --------------------------------------------------
# ---------------- LOADING DATA --------------------
# --------------------------------------------------
#.. Individual data ..
ind.df <- read.csv("inputdata/HerringData.csv", 
                   sep= ",",
                     header = TRUE, 
                     na.strings = ".")
ind.df$yearclass <- ind.df$year - ind.df$age

# .. N at age and year ..
N.df <- read.table("inputdata/N.txt",
                   header = TRUE)
names(N.df) <- str_remove(names(N.df), "X")

# load temperature data
if("kola.RData" %in% list.files(path = "inputdata/")) {
  load("inputdata/kola.RData")
}

# -- Standardizing data resources --
N.df$age <- as.integer(rownames(N.df))
N.df <- gather(N.df, "year", "N", -age)
N.df$year <- as.integer(N.df$year)
N.df$yearclass <- N.df$year-N.df$age

N.df <- N.df %>% group_by(year) %>% mutate(TN = sum(N))

N.df <- filter(N.df, age %in% 3:4) %>% group_by(year) %>%
  mutate(Nsum = sum(N)) 

# -- cohort size is number of fish at age 3 --
N.df <- filter(N.df, age == 3) %>% select(N, Nsum,TN, yearclass)
N.df$year = NULL
# -- putting together the data --
df <- left_join(ind.df, N.df, by = c("yearclass"))

N.df2 <- N.df %>% mutate(yearclass = yearclass +1) %>% rename(N2 = N) %>% 
  select(yearclass, N2)
df <- left_join(df, N.df2, by = c("yearclass"))

recruit.df <- mutate(xsam.df, yearclass = year-2)  %>% 
  select(yearclass, recruit)
  
df <- left_join(df, recruit.df, by = "yearclass")

# --- filtering out restrictions from the working team ---
#.. 1935-2019
#.. maturity 1-5
#.. > 9 individuals per age per year
range(df$maturity, na.rm = T)
range(df$lon, na.rm = T)
range(df$lat, na.rm = T)
range(df$year, na.rm = T)
range(df$month[df$age > 3], na.rm = T)
range(df$month[df$age <= 3], na.rm = T)

# -- less than 9 observations removed:--
df <- left_join(df,
                df %>% group_by(age,year) %>% count(age, year) ,
                by = c("age","year"))
df <- filter(df, n >9)

# -- removing African observations: --
df2 <- filter(df, between(lat, 10, 78)|is.na(lat)) # - 4302 fish

# -- remove observations on land in Russia and Sweeden/Finland --
df2 <- filter(df2, !(between(lon, 45,46) &  between(lat, 61,62))) # Russia -1994 fish
df2 <- filter(df2, !(between(lon,12.2,15.6) & between(lat, 58.25, 62.5))) # Sweeden -174 fish
df2 <- filter(df2, !(between(lon,22,24) & between(lat, 68,69))) #Finland (1994) - 7 fish

df <- df2
rm(df2)

# -- Grouping of gear codes: --
df$gearcode <- factor(ifelse(df$gear %in% c(3700,3710,3711,3713,3714,3720), "Purse Seine", 
                      ifelse(df$gear %in% c(
                        3100,3110,3120,3130,3194,3230,3231,3235,3236,3270,3271,3293,3400,3410, 3411, 3412, 3415, 3500, 3510,  3511, 3512, 3513,3514,3516,3518,3520,3521, 3522, 3530, 3531, 3532, 3533, 3535, 3536, 3541, 3547, 3590
                      ),"Trawl","Other")))

# We only have N for years <= 2017: 
df <- filter(df, year < 2018)
df <- filter(df, !is.na(N))
df <- filter(df, !is.na(lat)) # 1 point

# -- age as double/continous number and not integer: --
df$julianage <- df$age + (df$month-1)/12

# -- matching mean summer temperature data to observations --
if("kola.RData" %in% list.files(path = "inputdata/")) {
  summer <- cbind(year = kola$year, temp = rowMeans(kola[, c("Jun","Jul","Aug")]),
                  aug = kola[, "Aug"])
  names(summer) <- c("year","temp", "aug")
  df$temp <-df$aug.temp <-  NULL
  for(i in 1:nrow(df)){
    s <- df$yearclass[i] + 0:min(df$age[i], 3)
    df$temp[i] <- mean(summer[which(summer[,"year"] %in% s), "temp"])
  }
  df <- filter(df, !is.na(temp))
}else{
  df$temp = rnorm(nrow(df), mean = 5, sd = 2)
}
save(df, file = "data/official_data.RData")
save(N.df, file = "data/Ndf.RData")
