library(splines)
library(lubridate)
library(qwraps2)
library(lmtest)
library(dvmisc)
library(ggplot2)
library(viridis)
library(ggthemes)
library(RColorBrewer)
library(ggspatial)
library(ggpubr)
#Temporary files
showTmpFiles()
rasterTmpFile()
rasterOptions(tmpdir="E:/R Temp Files")
removeTmpFiles(h=0)
work.semen <- work.semen[order(work.semen$sid, work.semen$date),]
work.semen$date0 <- work.semen$date - 90
a <- merge(work.semen, work.baseline, by.x = "sid", by.y = "sid")
#GD border
a <- st_read("D:/Research Projects/地图信息等/中国地图/sheng/CN-sheng-A.shp")
sheng <- a[a$SHENG_ID == 44,]
plot(sheng[1,7])
shi <- st_read("D:/Research Projects/地图信息等/中国地图/shi/CN-shi-A.shp")
shi <- shi[shi$sheng == 440000,]
shi <- st_transform(shi, crs = crs)

#Geocode
address <- st_as_sf(address, coords = c('longitude.wgs', 'latitude.wgs'), crs = crs(ndvi))
address <- st_transform(address, crs = crs)
#stack####
list <- list.files(path = "Greenness of GD", pattern = '.NDVI', all.files = T, full.names = T)
ndvi <- stack(list)
writeRaster(pndvi, "pndvi", format = 'GTiff', overwrite = T)
#date
rastertime <- gsub("Greenness of GD/MOD13Q1.A", "", list) %>%
  gsub("\\.mosaic\\.006\\.\\d*\\.psmcrpgs_\\d*\\.250m_16_days_NDVI-250m_16_days_NDVI.tif", "", .) 
rastertime <- as.Date(rastertime, format = "%Y%j")
names(ndvi) <- as.character(rastertime)
#NDVI values
writeRaster(ndvi, "ndvi", format = 'GTiff', overwrite = T)
writeRaster(pndvi, "pndvi", format = 'GTiff', overwrite = T)
pndvi <- brick("pndvi.tif")
ndvi <- calc(a, fun = function(x) ifelse(x == -3000, NA, x*0.0001))
pndvi <- ndvi
pndvi <- calc(pndvi, fun = function(x) ifelse(x < 0.2, 0, x))
#Sheng border
a <- st_read("C:/Users/qq502/Desktop/地图信息等/中国地图/sheng/CN-sheng-A.shp")
sheng <- a[a$SHENG_ID == 44,]
sheng <- st_transform(sheng, crs =  "+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#CRS
crs = "+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ndvi <- projectRaster(ndvi, crs = crs)
plot(ndvi[[1]])
plot(address[,3], add = T)
plot(sheng[1,7], add = T)
#get NDVI by time and address (change buffer and output name)(crs transformed here)
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach (i = 1:nrow(address1), .packages = c("raster", "sf"), .combine = "rbind") %dopar% {
    a <- extract(pndvi, address1[i,], buffer = 400, fun = mean, na.rm = T)
  }
  a <- data.frame(t(a))
  a <- foreach(i = 1:ncol(a), .packages = c("openair"), .combine = "c") %dopar% {
    mean(selectByDate(data.frame(date = rastertime, a[,i]), start = sample$date0[i], end = sample$date[i])[,2], na.rm = T)
  }
  stopCluster(cl)
  pndvi400 <- a
})
sample0914$mndvi400 <- mndvi400
sample0914$mndvi400.cut <- quant_groups(sample0914$mndvi400, 4)
#Alternative
names(pndvi) <- format(rastertime, "X%Y%m%d")
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach (i = 1:nrow(address1), .packages = c("raster", "sf"), .combine = "rbind") %dopar% {
    a <- extract(pndvi, address1[i,], buffer = 400, fun = mean, na.rm = T)
  }
  stopCluster(cl)
})
c(paste0("X", format(seq(from = as.Date("2015-09-14"), to = as.Date("2015-11-17"), by = 16), "%Y%m%d")))
pndvi[[c(paste0("X", format(seq(from = as.Date("2015-09-14"), to = as.Date("2015-11-17"), by = 16), "%Y%m%d")))]]
#Total motility
sample$tm <- sample$pr +sample$np
#Semen volume
sample <- semen[!is.na(semen$semen.volume),]
#sperm number
sample$sperm.number <- sample$semen.volume*sample$sperm.concentration
#Total Motile Sperm Count
sample$tmsc <- sample$tm*sample$sperm.number
#Normal forms
sample1 <- sample1[!is.na(sample1$normal.forms),]
#Edu
semen$education <- factor(semen$education, levels = c("中专","初中","高中","大专","本科","硕士","博士","博士后"))
levels(semen$education)
new.levels <- c("hsl","hsl","hsl","diploma","bh","bh","bh","bh" )
semen$edu <- factor(new.levels[semen$education])
#Race-Han
semen$race <- factor(a$race)
semen$han <- 'han'
semen[semen$race != "汉族" & semen$race != "汉", ]$han <- 'no'
semen$han <- factor(semen$han)
summary(semen$han)
#Season 6-11
semen$season <- 'cold'
semen[month(semen$date) == 6,]$season <- 'warm'
semen$season <- factor(semen$season)
summary(semen$season)
#Age
a <- as.Date(semen$date) - as.Date(semen$dob)
semen$age <- as.numeric(floor(a/365))
sample$age.bi <- quant_groups(sample$age, 2)
#Abs
semen1$abs.cut <- "(2,3)"
semen1[semen1$abs %in% c(4,5),]$abs.cut <- "(4,5)"
semen1[semen1$abs %in% c(6,7),]$abs.cut <- "(6,7)"
semen1$abs.cut <- factor(semen1$abs.cut)
#Road length (G,S,X)
road.all = st_read("C:/Users/qq502/Desktop/地图信息等/china-latest-free.shp/gis_osm_roads_free_1.shp")
road <- road.all[road.all$fclass %in% c("motorway", "trunk", "primary", "secondary", "motorway_link", "trunk_link", "primary_link", "secondary_link"),]
bbox <- st_bbox(c(xmin = 109, xmax = 117.3, ymax = 25.7, ymin = 20), crs = 4326) %>%
  st_as_sfc() 
road.GD0 <- st_transform(st_intersection(road, bbox), crs= crs) 
system.time({
  road.length.tmp <- st_intersection(road.GD0[, 11], st_buffer(address1[,1], dist = 400))
  road.length.tmp$road.length <- st_length(road.length.tmp)
  length400.1 <- aggregate(road.length ~ sid, data = road.length.tmp, function(x) as.numeric(sum(x))/1000)
})
a <- merge(semen1, length, by.x = "sid", by.y = "sid", all = T)
sample1[is.na(sample1$length400),]$length400 <- 0
sample0914$length400.cut <- quant_groups(sample0914$length400, 4)
#For the first specimen
address2 <- address1 %>%
  group_by(sid) %>%
  .[order(sample$date),] %>%
  filter(row_number() == 1)
address2 <- st_as_sf(data.frame(address2))
semen1$pndvi500 <- a
#Analysis effect####volume/concentration/number/tm/pr/tmsc/forms
sample0914 %>%
  lmer(tmsc ~ mndvi400.cut + edu + han + abs.cut + age.cut + length400.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid),.) %>%
  summary() 
I(pndvi400/IQR(sample0914$pndvi400))

sample0914%>%
  lm(sperm.concentration ~ I(mndvi400/IQR(sample0914$mndvi400)) + edu + han + abs.cut + age.cut + length1500.cut + ns(tmp, df = 3) + ns(rh, df = 3), data = .) %>%
  summary()

library(car)
vif(a)
#Results display####
#Table 3
system.time({
  y <- c(39,40,41,42,35,44)
  x <- c(32,26,29)
  z <- c(21,17,19)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + sample0914[,z[i]] + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[[10]] %>%
        .[2, 1]
      b <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + sample0914[,z[i]] + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        confint(.) %>%
        .[4,]
      data.frame(a, b[1], b[2], row.names = paste(colnames(sample0914)[x[i]], colnames(sample0914)[y[j]]))
    }}
  stopCluster(cl)
})
write.table(round(a, 3), file = "a", sep = ",", row.names = F, col.names = T)
system.time({
  y <- c(39,40,41,42,35,44)
  x <- c(33,27,30)
  z <- c(21,17,19)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ sample0914[,x[i]] + sample0914[,z[i]] + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[[10]] %>%
        .[c(2:4), 1]
      b <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ sample0914[,x[i]] + sample0914[,z[i]] + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        confint(.) %>%
        .[c(4:6),]
      data.frame(a, b, row.names = paste(colnames(sample0914)[x[i]], colnames(sample0914)[y[j]], c("cut2", "cut3", "cut4")))
    }}
  stopCluster(cl)
})
write.table(round(a, 2), file = "a", sep = ",", row.names = T, col.names = T)
#p trend
a <- sample0914 %>%
  group_by(pndvi1500.cut) %>%
  summarize(pndvi1500.cut.median = median(pndvi1500))
b <- merge(b, a, by.x = "pndvi1500.cut", by.y = "pndvi1500.cut")
summary(factor(b$pndvi400.cut.median))
summary(sample0914$pndvi400.cut)
a <- arrange(b,sid,date)
b <- a[,order(colnames(a))]
a <- b %>% relocate(sid, date)
sample0914 <- a


system.time({
  y <- c(39,40,41,42,35,44)
  x <- c(34,28,31)
  z <- c(21,17,19)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ sample0914[,x[i]] + sample0914[,z[i]] + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[[10]] %>%
        .[2, 5]
    }}
  stopCluster(cl)
})
write.table(round(a, 3), file = "a", sep = ",", row.names = F, col.names = T)
#Table 1
summary(semen1[!is.na(semen1$normal.forms),]$pndvi400.cut)

p <- function(i,j) {n_perc0(sample0914[sample0914$pndvi400.cut == levels(sample0914$pndvi400.cut)[i],]$abs.cut == levels(sample0914$abs.cut)[j], digits = 1)}
a <- t(sapply(c(1:4), function(i) sapply(c(1:4), p , i)))
write.table(a, file = "a", quote = F, sep = ";", row.names = T, col.names = T)
sapply(c(1:3), function(j) p(2, j))

round(summary(sample0914$abs.cut)/387.54, 2)
mean_sd(sample0914[sample0914$pndvi400.cut == levels(sample0914$pndvi400.cut)[4],]$rh)
n_perc(factor(year(a[a$pndvi400.cut == levels(a$pndvi400.cut)[4],]$date)) == 2019)

a <- sample0914 %>%
  count(sid)
summary(a$n == 1) 
a[a$n ==38,]
#Table 2
a <- data.frame(
  mean = round(sapply(c(39,40,41,42,35,44,24), function(x) mean(sample0914[,x], na.rm = T)), 3),
  sd = round(sapply(c(39,40,41,42,35,44,24), function(x) sd(sample0914[,x], na.rm = T)), 3),
  t(round(sapply(c(39,40,41,42,35,44,24), function(x) quantile(sample0914[,x], probs = c(0.01,0.25,0.5,0.75,0.99), na.rm = T)), 3))
)
write.table(a, file = "a", quote = F, sep = ";", row.names = T, col.names = T)
#Table 4 Stratified
library(ggplot2)
ggplot(sample, aes(x=date, y = rh)) +
  geom_point()

I(pndvi400/IQR(sample$pndvi400, na.rm = T))
#Age
sample0914%>%
  lmer(semen.volume ~ I(pndvi400/IQR(sample0914$pndvi400)) + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + fertility + (1|sid), data = ., subset = age.bi == levels(sample0914$age.bi)[2]) %>%
  summary() 
sample0914 %>%
  lm(normal.forms ~ I(pndvi400/IQR(sample0914$pndvi400)) + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + fertility, data = ., subset = age.bi == levels(sample0914$age.bi)[1]) %>%
  summary()
sample0914%>%
  lmer(tm ~ I(pndvi400/IQR(sample0914$pndvi400)) + han + abs.cut + age.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + (1|sid), data = ., subset = fertility == levels(sample0914$fertility)[2]) %>%
  summary() 
sample0914 %>%
  lm(sperm.concentration ~ I(pndvi400/IQR(sample0914$pndvi400)) + edu + han + age.cut + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut, data = ., subset = fertility == levels(sample0914$fertility)[1]) %>%
  summary()
with(sample0914, prop.table(table(age.bi, fertility)))
library(car)
vif(sample0914 %>%
      lm(normal.forms ~ I(pndvi400/IQR(sample0914$pndvi400)) + edu + han + age.bi + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + marital + fertility, data = .))
y <- c(39,40,41,42,35,44)
lrtest(lm(normal.forms ~ I(pndvi400/IQR(sample0914$pndvi400))*age.bi + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut, sample0914),
       lm(normal.forms ~ I(pndvi400/IQR(sample0914$pndvi400)) + age.bi + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut, sample0914)
)
lrtest(lmer(semen.volume ~ I(pndvi400/IQR(sample0914$pndvi400))*age.bi + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + (1|sid), data = sample0914),
       lmer(semen.volume ~ I(pndvi400/IQR(sample0914$pndvi400)) + age.bi + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + (1|sid), data = sample0914)
)
a <- function(x) {lrtest(lmer(sample0914[,y[x]] ~ I(pndvi400/IQR(sample0914$pndvi400))*age.bi + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + (1|sid), data = sample0914),
       lmer(sample0914[,y[x]] ~ I(pndvi400/IQR(sample0914$pndvi400)) + age.bi + han + abs.cut + edu + ns(tmp, df = 3) + ns(rh, df = 3) + length400.cut + (1|sid), data = sample0914)
)}
a(6)

system.time({
  y <- c(39,40,41,42,35,44)
  x <- c(32,26,29)
  z <- c(21,17,19)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + sample0914[,z[i]] + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), ., subset = age.bi == levels(sample0914$age.bi)[2]) %>%
        summary(.) %>%
        .[[10]] %>%
        .[2, 1]
      b <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + sample0914[,z[i]] + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), ., subset = age.bi == levels(sample0914$age.bi)[2]) %>%
        confint(.) %>%
        .[4,]
      data.frame(a, b[1], b[2], row.names = paste(colnames(sample0914)[x[i]], colnames(sample0914)[y[j]]))
    }}
  stopCluster(cl)
})
write.table(round(a, 2), file = "a", sep = ",", row.names = F, col.names = T)
#Season
summary(factor(month(length07011[length07011$season == 'cold',]$date)))
sample%>%
  lmer(pr ~ I(pndvi400/IQR(sample$pndvi400, na.rm = T)) + han + abs.cut + edu + age.cut + ns(rh, df = 3) + length400.cut + (1|sid), data = ., subset = season == 'cold') %>%
  summary() 
sample1 %>%
  lm(normal.forms ~ I(pndvi400/IQR(sample1$pndvi400, na.rm = T)) + edu + han + abs.cut + age.cut + ns(rh, df = 3) + length400.cut, data = ., subset = season == 'warm') %>%
  summary()
summary(lm(normal.forms ~ pndvi400*season + han + abs + edu + age_cut + ns(tmp, df = 3) + ns(rh, df = 3) + length400_cut + as.character(year(date)), length07011))

lrtest(lm(normal.forms ~ I(pndvi400/IQR(sample1$pndvi400, na.rm = T))*season + han + abs.cut + edu + age.cut + ns(rh, df = 3) + length400.cut, sample1),
       lm(normal.forms ~ I(pndvi400/IQR(sample1$pndvi400, na.rm = T)) + season + han + abs.cut + edu + age.cut + ns(rh, df = 3) + length400.cut, sample1)
)
lrtest(lmer(sperm.concentration ~ I(pndvi400/IQR(sample$pndvi400, na.rm = T))*season + han + abs.cut + age.cut + edu + ns(rh, df = 3) + length400.cut + (1|sid), data = sample),
       lmer(sperm.concentration ~ I(pndvi400/IQR(sample$pndvi400, na.rm = T)) + season + han + abs.cut + age.cut + edu + ns(rh, df = 3) + length400.cut + (1|sid), data = sample)
)
a <- function(x){lrtest(lmer(x ~ I(pndvi400/IQR(sample$pndvi400, na.rm = T))*season + han + abs.cut + age.cut + edu + ns(rh, df = 3) + length400.cut + (1|sid), data = sample),
       lmer(x ~ I(pndvi400/IQR(sample$pndvi400, na.rm = T)) + season + han + abs.cut + age.cut + edu + ns(rh, df = 3) + length400.cut + (1|sid), data = sample)
)}
a(sample[,y[5]])

system.time({
  y <- c(42,43,44,45,38,47)
  x <- c(35,29,32,24)
  z <- c(21,17,19,21)
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      a <- sample %>%
        lmer(sample[,y[j]] ~ I(sample[,x[i]]/IQR(sample[,x[i]], na.rm = T)) + sample[,z[i]] + edu + han + abs.cut + age.cut + ns(rh, df = 3) + (1|sid), ., subset = season == 'warm') %>%
        summary(.) %>%
        .[[10]] %>%
        .[2, 1]
      b <- sample %>%
        lmer(sample[,y[j]] ~ I(sample[,x[i]]/IQR(sample[,x[i]], na.rm = T)) + sample[,z[i]] + edu + han + abs.cut + age.cut + ns(rh, df = 3) + (1|sid), ., subset = season == 'warm') %>%
        confint(.) %>%
        .[4,]
      data.frame(a, b[1], b[2], row.names = paste(colnames(sample)[x[i]], colnames(sample)[y[j]]))
    }}
  stopCluster(cl)
})
write.table(round(a, 2), file = "a", sep = ",", row.names = T, col.names = T)

#Exposure plot
#Mean NDVI in 2016-2019, in GD
ndvi.mean <- calc(pndvi, fun=mean, na.rm = T)
writeRaster(ndvi.mean, "ndvimean", format = 'GTiff', overwrite = T)
ndvi.mean <- raster("ndvimean.tif")
sheng.sp <- as_Spatial(sheng$geometry)
ndvi.gd <- mask(ndvi.mean, sheng.sp)
writeRaster(ndvi.gd, "ndvimeangd", format = 'GTiff', overwrite = T)
a <- as(ndvi.gd, "SpatialPixelsDataFrame")
a <- as.data.frame(a)
colnames(a) <- c("value", "x", "y")
ndvi.gd.df <- a
ndvi.gd.df$cut <- cut(ndvi.gd.df$value, breaks = c(0,0.2,0.3,0.5,0.7, 1), labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8",">=0.8"), right = F, ordered_result = T)
#Mean NDVI in 400m for each subject
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach (i = 1:nrow(address2), .packages = c("raster", "sf"), .combine = "c") %dopar% {
    a <- extract(ndvi.mean, address2[i,], buffer = 400, fun = mean, na.rm = T)
  }
  stopCluster(cl)
  address2$ndvi.mean <- a
})
address2$ndvi.mean.cut <- cut(address2$ndvi.mean, breaks = c(0,0.2,0.3,0.5,0.7,1), labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8",">=0.8"), right = F, ordered_result = T)
#Transfer sf points to raster
address2.sp <- as_Spatial(address2)
address2.sp <- SpatialPoints(address2.sp)
address2.r <- rasterize(address2.sp, address2.r, fun = "count")
address2.r <- aggregate(address2.r, fact = 2)
address2.df <- data.frame(as(address2.r, "SpatialPixelsDataFrame"))
address2.df$cut <- cut(address2.df$layer, breaks = c(1,11,31,51,71,101,151,201,301,1000), labels = c("1-10","11-30","31-50","51-70","71-100","101-150","151-200","201-300",">300"), right = F, ordered_result = T)
#Plot 1
a <- ggplot() +
  #geom_sf(data = map.hubei.city, color = NA, fill = "gray98", size = 0.2) +
  geom_sf(data = shi, color = NA, fill = "gray99") +
  #geom_tile(data = ndvi.gd.df, aes(x = x, y = y, fill = cut), alpha = 0.8) +
  geom_tile(data = address2.df, aes(x=x, y=y, fill = cut), alpha = 0.7) +
  #geom_sf(data = hamm.map, stroke = 0, color = "#993333", size = 0.6, alpha = 0.3) +
  #geom_sf(data = address2, color = "#045a8d", size = 1.3, alpha = 0.7, shape = 1) +  
  geom_sf(data = shi, color = "gray33", fill = NA, size = 0.1) +
  geom_sf(data = shi.gd, color = "black", fill = NA, size = 0.3) +
  coord_sf(datum = 4326) +
  annotation_scale(location = "br", 
                   width_hint = 0.2, 
                   line_width = 0.75, 
                   line_col = "gray22",
                   height = unit(0.1, "cm"), 
                   pad_x = unit(3.5, "cm"), 
                   pad_y = unit(0.3, "cm"), 
                   text_cex = 0.4, 
                   text_col = "gray22", 
                   bar_cols = c("white", "gray22")) +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(0.4, "cm"), width = unit(0.4, "cm"),
                         pad_x = unit(3.8, "cm"), pad_y = unit(0.75, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 7, text_col = "gray22", line_col = "gray22", fill = c("white", "gray22"))) +
  scale_fill_manual(values = rev(brewer.pal(9, "RdYlBu"))) +
  #scale_fill_manual(values = c("#fee391","#c7e9c0","#74c476","#238b45","#00441b"))+
  #scale_fill_manual(values = c("#08306b","#2171b5","#6baed6","#a1d99b","#fee391","#fcae91","#fb6a4a","#de2d26","#a50f15")) +
  #labs(fill = "NDVI values: ")+
  labs(fill = "No. of subjects:") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "gray88", size = 0.1),
        #axis.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(color = "gray33", size = 5),
        axis.ticks = element_blank(),
        legend.position = c(0.85,0.30),
        #legend.position = "none",
        #legend.direction = "horizontal",
        legend.key.width = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.title = element_text(size = 4, color = "gray22"),
        legend.text = element_text(size = 4, color = "gray22"),
        legend.margin = margin(t = 0.1, unit = "lines"))
#Plot 2
b <- ggplot() +
  #geom_sf(data = map.hubei.city, color = NA, fill = "gray98", size = 0.2) +
  geom_sf(data = shi, color = NA, fill = "gray99") +
  geom_tile(data = ndvi.gd.df, aes(x = x, y = y, fill = cut), alpha = 0.8) +
  #geom_tile(data = address2.df, aes(x=x, y=y, fill = cut), alpha = 0.7) +
  #geom_sf(data = hamm.map, stroke = 0, color = "#993333", size = 0.6, alpha = 0.3) +
  #geom_sf(data = address2, color = "#045a8d", size = 1.3, alpha = 0.7, shape = 1) +  
  geom_sf(data = shi, color = "gray33", fill = NA, size = 0.1) +
  geom_sf(data = shi.gd, color = "black", fill = NA, size = 0.3) +
  coord_sf(datum = 4326) +
  annotation_scale(location = "br", 
                   width_hint = 0.2, 
                   line_width = 0.75, 
                   line_col = "gray22",
                   height = unit(0.1, "cm"), 
                   pad_x = unit(3.5, "cm"), 
                   pad_y = unit(0.3, "cm"), 
                   text_cex = 0.4, 
                   text_col = "gray22", 
                   bar_cols = c("white", "gray22")) +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(0.4, "cm"), width = unit(0.4, "cm"),
                         pad_x = unit(3.8, "cm"), pad_y = unit(0.75, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 7, text_col = "gray22", line_col = "gray22", fill = c("white", "gray22"))) +
  #scale_fill_manual(values = rev(brewer.pal(10, "RdYlBu"))) +
  scale_fill_manual(values = c("#fee391","#c7e9c0","#74c476","#238b45","#00441b"))+
  #scale_fill_manual(values = c("#08306b","#2171b5","#6baed6","#a1d99b","#fee391","#fcae91","#fb6a4a","#de2d26","#a50f15")) +
  labs(fill = "NDVI values: ")+
  #labs(fill = "No. of subjects:") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "gray88", size = 0.1),
        #axis.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(color = "gray33", size = 5),
        axis.ticks = element_blank(),
        legend.position = c(0.85,0.30),
        #legend.position = "none",
        #legend.direction = "horizontal",
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.30, "cm"),
        legend.title = element_text(size = 5, color = "gray22"),
        legend.text = element_text(size = 5, color = "gray22"),
        legend.margin = margin(t = 0.1, unit = "lines"))
ggarrange(a, b)
