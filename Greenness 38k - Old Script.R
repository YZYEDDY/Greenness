library(raster)
library(sf)
library(dplyr)
library(rgdal)
library(foreach)
library(parallel)
library(doParallel)
library(openair)
library(lmerTest)
library(splines)
library(lubridate)
library(qwraps2)
library(lmtest)
length <- readRDS("length")
a <- merge(length0701, length, by.x = "sid", by.y = "sid", all = T)
length0701$length1500 <- a$road.length
a <- count(length0701, sid)
length07011 <- length07011[order(length07011$sid, length07011$date),]
#stack####
list <- list.files(path = "Greenness of GD", pattern = '.NDVI', all.files = T, full.names = T)
allndvi <- stack(list)
allndvi <- allndvi/10000
#address
address <- dplyr::select(semen, "sid", "date", "longitude.wgs", "latitude.wgs")
address <- st_as_sf(address, coords = c("longitude.wgs", "latitude.wgs"), crs = crs(allndvi))
#project
allndvi.pr <- projectRaster(allndvi, crs = "+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
address.pr <- st_transform(address, "+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#date
rastertime <- gsub("Greenness of GD/MOD13Q1.A", "", list) %>%
  gsub("\\.mosaic\\.006\\.\\d*\\.psmcrpgs_\\d*\\.250m_16_days_NDVI-250m_16_days_NDVI.tif", "", .) 
rastertime <- as.Date(rastertime, format = "%Y%j")
names(allndvi.pr) <- as.character(rastertime)
#get NDVI####
brick(positivendvi.pr, "pndvi")
allndvi.pr <- brick("allndvi.tif")
#Convert ndvi ("<" converts the values selected to 1/TRUE)
allndvi.pr[allndvi.pr == -0.3] <- NA
positivendvi.pr <- allndvi.pr
positivendvi.pr[positivendvi.pr < 0.2] <- 0
#get NDVI by time and address (change buffer and output object)
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach (i = 1:nrow(address.pr), .packages = c("raster", "sf"), .combine = "rbind") %dopar% {
    extract(allndvi.pr, address.pr[i,], buffer = 400, fun = mean)
    }
  a <- data.frame(t(a))
  a <- foreach(i = 1:ncol(a), .packages = c("openair"), .combine = "c") %dopar% {
  mean(selectByDate(data.frame(date = rastertime, a[,i]), start = period1819[i,1], end = period1819[i,2])[,2])
  }
  stopCluster(cl)
   b <- a
})
a <- extract(allndvi.pr, address.pr[2,], buffer = 400, fun = mean)
a <- data.frame(t(a))
mean(selectByDate(data.frame(date = rastertime, a), start = start[2], end = end[2])[,2])

#Education
levels(length0701$education)
new.levels <- c("中专高中", "中专高中", "大专", "本科以上", "本科以上", "本科以上", "本科以上")
length0701$edu <- factor(new.levels[length0701$education])
#Season 6-11
length0701$season <- 'cold'
length0701[month(length0701$date) == 10,]$season <- 'warm'
length0701$season <- factor(length0701$season)
summary(length07011$season)
#Road length
length07011[,38:45][is.na(length07011[,38:45])] <- 0
length07011[,63][is.na(length07011[,63])] <- 0
library(dvmisc)
length0701$length1500_cut <- quant_groups(length0701$age, 2)
#For the first specimen
a <- length0701 %>%
  group_by(sid) %>%
  .[order(length0701$date),] %>%
  filter(row_number() == 1)
length07011$length1500 <- a$length1500
#Analysis total effect####
length0701%>%
  lmer(sperm.concentration ~ pndvi400 + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid),.) %>%
  summary() 
I(pndvi400/IQR(length0701$pndvi400))
length07011 %>%
  lm(normal.forms ~ pndvi400.cut + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), data = .) %>%
  summary()
#Basic characteristics
library(qwraps2)
#Likelihood ratio test for random term
rand(lmer(semen.volume ~ ndvi400.cut + education + race + age_cut + pm2.5 + (1|sid), length0701[length0701$confidence >= 70 | length0701$comprehension >= 20,]))
#Correlation
cor.test(length0701$pm2.5, length0701$length3000, method = c("pearson", "kendall", "spearman"))
library(ellipse)
a <- cor(select(length0701, "ndvi500", "ndvi1000","pndvi400", "pndvi500", "pndvi1000", "pm2.5", "rh", "tmp", "length100", "length300", "length400", "length500", "length1000", "length2000","length3000", "length5000", "semen.volume", "sperm.concentration", "pr", "normal.forms"), use = "complete.obs")
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab") 
plotcorr(round(a, digits = 2), col=rgb(colorfun((a+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1))
#Generalized additive mixed model (to be completed)####
library(mgcv)
length0701 %>%
  gamm(normal.forms ~ s(pndvi400) + edu + han + age_cut + s(tmp) + s(rh) + s(length1000) + as.character(year(date)), data = ., random=list(sid=~1)) %>%
  summary() 



#Results display####
#Estimates for all observations
system.time({
  y <- c(17, 18, 15)
  x <- c(49)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- length0701 %>%
        lmer(length0701[,y[j]] ~ length0701$ndvi400.cut + length0701[,x[i]] + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + as.character(year(date)) + (1|sid), .) %>%
        summary(.) %>%
        .[[10]] %>%
        .[c(2:4), 1]
      b <- length0701 %>%
        lmer(length0701[,y[j]] ~ length0701$ndvi400.cut + length0701[,x[i]] + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + as.character(year(date)) + (1|sid), .) %>%
        confint(.) %>%
        .[c(4:6),]
      data.frame(a, b, row.names = paste(i, j, c("cut2", "cut3", "cut4")))
    }}
  stopCluster(cl)
})
write.table(round(a, 2), file = "a", sep = ",", row.names = T, col.names = T)
#p trend
a <- length0701 %>%
  group_by(ndvi400.cut) %>%
  summarize(ndvi400.cut.median = median(ndvi400))
b <- merge(length0701, a, by.x = "ndvi400.cut", by.y = "ndvi400.cut")
summary(factor(b$ndvi400.cut.median))
summary(length0701$pndvi1500.cut)
a <- length07011 %>%
  group_by(ndvi400.cut) %>%
  summarize(ndvi400.cut.median = median(ndvi400))
b <- merge(length07011, a, by.x = "ndvi400.cut", by.y = "ndvi400.cut")
summary(factor(b$ndvi400.cut.median))
summary(length07011$pndvi1500.cut)
#table 2 
mean_sd(length0701[length0701$pndvi400.cut == "[0,0.13]",]$length400, na_rm = T)
mean_sd(length0701[length0701$pndvi400.cut == "(0.13,0.211]",]$length400, na_rm = T)
mean_sd(length0701[length0701$pndvi400.cut == "(0.211,0.3]",]$length400, na_rm = T)
mean_sd(length0701[length0701$pndvi400.cut == "(0.3,0.787]",]$length400, na_rm = T)

median_iqr(length07011[length07011$pndvi400.cut == "[0,0.129]",]$normal.forms, na_rm = T)
median_iqr(length07011[length07011$pndvi400.cut == "(0.129,0.21]",]$normal.forms, na_rm = T)
median_iqr(length07011[length07011$pndvi400.cut == "(0.21,0.306]",]$normal.forms, na_rm = T)
median_iqr(length07011[length07011$pndvi400.cut == "(0.306,0.763]",]$normal.forms, na_rm = T)
#Stratification####
library(ggplot2)
ggplot(length0701, aes(x=date, y = pndvi400)) +
  geom_point()

length0701$season2 <- "cold"
length0701[month(length0701$date)> 6 & month(length0701$date) < 13,]$season2 <- "warm"
length0701$season2 <- factor(length0701$season2)

I(pndvi400/IQR(length0701$pndvi400))
#Age
length0701%>%
  lmer(pr ~ pndvi400.cut.median + han + abs + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid), data = ., subset = age_bi == "low") %>%
  summary() 
length07011 %>%
  lm(normal.forms ~ pndvi400.cut + edu + han + abs + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), data = ., subset = age_bi == "high") %>%
  confint()
summary(lm(normal.forms ~ pndvi400*age_bi + han + abs + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011))
lrtest(lm(normal.forms ~ pndvi400 + age_bi + han + abs + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011),
       lm(normal.forms ~ pndvi400 + edu + han + abs + age + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011)
)
lrtest(lmer(semen.volume ~ pndvi400*age_bi + han + abs + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid), length0701),
       lmer(semen.volume ~ pndvi400 + age_bi + edu + han + abs + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid), length0701)
)

system.time({
  y <- c(16, 17, 14)
  x <- c(48)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- length0701 %>%
        lmer(length0701[,y[j]] ~ length0701$pndvi400.cut + length0701[,x[i]] + edu + han + abs + ns(tmp, df = 3) + ns(rh, df = 3) + as.character(year(date)) + (1|sid), ., subset = age_bi == "high") %>%
        summary(.) %>%
        .[[10]] %>%
        .[c(2:4), 1]
      b <- length0701 %>%
        lmer(length0701[,y[j]] ~ length0701$pndvi400.cut + length0701[,x[i]] + edu + han + abs + ns(tmp, df = 3) + ns(rh, df = 3) + as.character(year(date)) + (1|sid), ., subset = age_bi == "high") %>%
        confint(.) %>%
        .[c(4:6),]
      data.frame(a, b, row.names = paste(i, j, c("cut2", "cut3", "cut4")))
    }}
  stopCluster(cl)
})
write.table(round(a, 2), file = "a", sep = ",", row.names = T, col.names = T)
#Season
summary(factor(month(length07011[length07011$season == 'cold',]$date)))
length0701%>%
  lmer(pr ~ I(pndvi400/IQR(length0701$pndvi400)) + han + abs + edu + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid), data = ., subset = season1 == 'cold') %>%
  summary() 
length07011 %>%
  lm(normal.forms ~ pndvi400.cut + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) +  length400_cut + as.character(year(date)), data = ., subset = season1 == 'warm') %>%
  summary()
summary(lm(normal.forms ~ pndvi400*season + han + abs + edu + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011))
lrtest(lm(normal.forms ~ pndvi400.cut*season1 + han + abs + edu + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011),
       lm(normal.forms ~ pndvi400.cut + season1 + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011)
)
lrtest(lmer(semen.volume ~ pndvi400*season1 + han + abs + edu + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid), length0701),
       lmer(semen.volume ~ pndvi400 + season1 + han + abs + edu + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)) + (1|sid), length0701)
)

system.time({
  y <- c(17, 18, 15)
  x <- c(49)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- length0701 %>%
        lmer(length0701[,y[j]] ~ length0701$pndvi400.cut + length0701[,x[i]] + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + as.character(year(date)) + (1|sid), ., subset = season1 == "warm") %>%
        summary(.) %>%
        .[[10]] %>%
        .[c(2:4), 1]
      b <- length0701 %>%
        lmer(length0701[,y[j]] ~ length0701$pndvi400.cut + length0701[,x[i]] + edu + han + abs + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + as.character(year(date)) + (1|sid), ., subset = season1 == "warm") %>%
        confint(.) %>%
        .[c(4:6),]
      data.frame(a, b, row.names = paste(i, j, c("cut2", "cut3", "cut4")))
    }}
  stopCluster(cl)
})
write.table(round(a, 2), file = "a", sep = ",", row.names = T, col.names = T)


