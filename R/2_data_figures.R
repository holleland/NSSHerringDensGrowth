rm(list=ls())
# figures
library(tidyverse)
library(sf)
theme_set(theme_bw())
# -- load data --
load("data/official_data.RData")
load("data/Ndf.RData")
med <- median(N.df$N)

# --- Figure 1: map of observations ---
Norway <- st_read("shapefile/ne_10m_land.shp")
Norway <- st_crop(Norway, c(xmin = min(df$lon, na.rm=T)-2, ymin = 45,
                            xmax = max(df$lon, na.rm=T)+2, ymax = max(df$lat, na.rm=T)+2))
p1 <- ggplot(Norway) + geom_sf(fill = "grey80") +
  geom_point(data = df[order(df$N, decreasing = TRUE),],
              aes(x = lon, y = lat, col = log(N)),
              size = .8)+
  scale_x_continuous(expand = c(0,0), limits = c(min(df$lon, na.rm=T)-2,max(df$lon, na.rm=T)+2)) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(55,80,5), 69),
                     limits = c( min(df$lat, na.rm=T)-2, max(df$lat, na.rm=T)+2)) + 
  scale_color_viridis_c(name  = "Log\ndensity")+
  guides(color = guide_colorbar(barheight = 25))+
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "light blue"))+
  ggtitle("Map of observations", paste0("N = ",nrow(df)))
p1
ggsave(p1, filename = "plots/1_map_of_observations.tiff",
       width=28, height=20, units="cm", device = "tiff", dpi = "retina")

# ----------------------
# ----- Figure S1 ------
# ----------------------
ncount <- df %>% 
  filter(age <4 &  !is.na(lat)) %>% 
  mutate(lat67 = factor(ifelse(lat<69,"<69°N",">=69°N"))) %>% 
  dplyr::select(yearclass, age, lat67) %>%
  group_by(yearclass,age,lat67) %>% 
  count()
agelabs <- paste0("Age ",1:3)
names(agelabs) <- 1:3
ncount$agefac <- factor(ncount$age)
ncount <- ncount %>% group_by(yearclass, age) %>% mutate(tot = sum(n))
ncount <- transform(ncount, lat67 = factor(lat67, levels = c(">=69°N","<69°N")))

(p1per <- ggplot(ncount, aes(x = yearclass, y = n/tot, fill = lat67)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~agefac, ncol = 1,
               labeller = labeller(agefac = agelabs)) +
    scale_x_continuous(expand = c(0.01,0),
                       breaks = seq(1930,2020,10), 
                       name = "Year class")+
    scale_y_continuous(expand = c(0,0), 
                       name = "Relative number of fish in dataset (in %)", 
                       labels = scales::percent) +
    scale_fill_discrete(name = "")+
    theme(strip.background = element_rect(fill = "white"), 
          legend.position = "top",
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
)

(p2 <- ggplot(filter(N.df, between(yearclass,1932,2014)),
              aes(x = yearclass, y = N/1000)) + 
    geom_bar(stat = "identity")+
    scale_x_continuous(expand = c(0.01,0),
                       breaks = seq(1930,2020,10), 
                       name = "Year class")+
    scale_y_continuous(expand = c(0,0), 
                       name = "Cohort size (in thousands)")+
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
    )
)
noIndivid <- df %>% group_by(yearclass) %>% count()
(p3 <- ggplot(filter(noIndivid, between(yearclass,1932,2014)),
             aes(x = yearclass, y = n/1000)) + 
  geom_bar(stat = "identity") + 
  scale_x_continuous(expand = c(0.01,0),
                     breaks = seq(1930,2020,10), 
                     name = "Year class")+
  scale_y_continuous(expand = c(0,0), 
                     name = "Number of samples (in thousands)")+
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  ))
library(ggpubr)
gg1per <- ggarrange(p1per,p2,p3,ncol = 1, heights = c(3,1,1))
ggsave(gg1per, filename = "plots/S1_relcounts_NorSof69.tiff",
       width=28, height=28*5/4, units="cm", device = "tiff", dpi = "retina")

# -----------------
# --- Figure 2  ---
# -----------------
med <- median(N.df$N)
df$logNfac <- factor(
  ifelse(log(df$N)<log(med),paste0("Low density (N[3] < ",round(med,1),")"),"High density")
)
df$logNfac <- factor(
  ifelse(log(df$N)<log(quantile(N.df$N, probs = .25)),paste0("Low density (N[3] < ",round(quantile(N.df$N, probs = .25),1),")"),"High density")
)
df <- transform(df, logNfac = factor(logNfac, levels = c(paste0("Low density (N[3] < ",round(quantile(N.df$N, probs = .25),1),")"),"High density")))
boxperage <- ggplot(df, aes(x=factor(age), y = length, col = logNfac)) + geom_boxplot() + 
  xlab("Age") + ylab("Length")+
  scale_color_discrete(labels = c(as.expression(bquote("Low density (" ~ N[(3)] ~"<"~.(round(quantile(N.df$N, probs = .25),1))~")")),
                                  "High density"))+
  theme(legend.title = element_blank(),
        legend.position = "top")
# -- Nsum --
df$logNfac2 <- factor(
  ifelse(log(df$Nsum)<log(quantile(N.df$Nsum, probs = .25)),paste0("Low density (N[(3)]+N[(4)] < ",round(quantile(N.df$Nsum, probs = .25),1),")"),"High density")
)
df <- transform(df, logNfac2 = factor(logNfac2, levels = c(paste0("Low density (N[(3)]+N[(4)] < ",round(quantile(N.df$Nsum, probs = .25),1),")"),"High density")))
boxperage2 <- ggplot(df, aes(x=factor(age), y = length, col = logNfac2)) + geom_boxplot() + 
  xlab("Age") + ylab("Length")+
  scale_color_discrete(labels = c(as.expression(bquote("Low density (" ~ N[(3)]+N[(4)] ~"<"~.(round(quantile(N.df$Nsum, probs = .25),1))~")")),
                                  "High density"))+
  theme(legend.title = element_blank(),
        legend.position = "top")
ggsave(ggarrange(boxperage,boxperage2, ncol = 1, nrow = 2), 
       filename = "plots/3_boxplots_length_per_age_Nsum_combined.tiff",
       width=20, height=1.7*14, units="cm", device = "tiff", dpi = "retina")

# -------------- 
# -- Latitude -- 
# -------------- 
df$latfac <- factor(
  ifelse(df$lat<=69, "South of 69°N","North of 69°N")
)
df <- transform(df, latfac = factor(latfac, levels = c("South of 69°N","North of 69°N")))
boxperage <- ggplot(df, aes(x=factor(age), y = length, col = latfac)) + geom_boxplot() + 
  xlab("Age") + ylab("Length")+
  theme(legend.title = element_blank(),
        legend.position = "top")
ggsave(boxperage, filename = "plots/S5_latitude_boxplots_length_per_age.tiff",
       width=20, height=14, units="cm", device = "tiff", dpi = "retina")

# --------------- 
# -- Gear code -- 
# --------------- 
df <- transform(df, gearcode = factor(gearcode, levels = c("Trawl","Purse Seine","Other")))
boxperage <- ggplot(df, aes(x=factor(age), y = length, col = gearcode)) + geom_boxplot() + 
  xlab("Age") + ylab("Length")+
  theme(legend.title = element_blank(),
        legend.position = "top")
ggsave(boxperage, filename = "plots/S8_gearcode_boxplots_length_per_age.tiff",
       width=20, height=14, units="cm", device = "tiff", dpi = "retina")

# -------------------
# --- Temperature ---
# -------------------
df$tempfac <- factor(
  ifelse(df$temp<round(mean(df$temp, na.rm =T),2), "Cold (< 4.35°C)","Warm (>= 4.35°C)")
)
df <- transform(df, tempfac = factor(tempfac, levels = c("Cold (< 4.35°C)","Warm (>= 4.35°C)")))
boxperage <- ggplot(filter(df,!is.na(temp)), aes(x=factor(age), y = length, col = tempfac)) + geom_boxplot() + 
  xlab("Age") + ylab("Length")+
  theme(legend.title = element_blank(),
        legend.position = "top")
ggsave(boxperage, filename = "plots/S4_boxplots_temp_length_per_age.tiff",
       width=20, height=14, units="cm", device = "tiff", dpi = "retina")
