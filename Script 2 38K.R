library(splines) # Regression Spline Functions and Classes
library(lubridate) # Make Dealing with Dates a Little Easier
library(qwraps2) # Quick Wraps 2
library(lmtest) # Testing Linear Regression Models
library(lmerTest) # lmer
library(dvmisc) # Convenience Functions, Moving Window Statistics, and Graphics
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(viridis) # Default Color Maps from 'matplotlib'
library(ggthemes) # Extra Themes, Scales and Geoms for 'ggplot2'
library(RColorBrewer) # ColorBrewer Palettes
library(ggspatial) # Spatial Data Framework for ggplot2 
library(ggpubr) # 'ggplot2' Based Publication Ready Plots
library(parallel) # makeCluster detectCores stopCluster
library(doParallel) # registerDoParallel
library(sf) # st_read st_transform st_as_sf %>% st_sf st_intersection st_buffer st_distance st_bbox st_as_sfc st_length as_Spatial

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
writeRaster(ndvi, "ndvi", format = 'GTiff', overwrite = T)
#date
rastertime <- gsub("Greenness of GD/MOD13Q1.A", "", list) %>%
  gsub("\\.mosaic\\.006\\.\\d*\\.psmcrpgs_\\d*\\.250m_16_days_NDVI-250m_16_days_NDVI.tif", "", .) 
rastertime <- as.Date(rastertime, format = "%Y%j")
names(ndvi) <- as.character(rastertime)

#NDVI values
writeRaster(ndvi, "ndvi", format = 'GTiff', overwrite = T)
writeRaster(pppndvi, "pppndvi", format = 'GTiff', overwrite = T)
pppndvi <- brick("pppndvi.tif")
ndvi <- calc(ndvi, fun = function(x) ifelse(x == -3000, NA, x*0.0001))
pppndvi <- calc(ndvi, fun = function(x) ifelse(x < 0, NA, x))

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
    a <- extract(pppndvi, address1[i,], buffer = 200, fun = mean, na.rm = T)
  }
  a <- data.frame(t(a))
  pppndvi200 <- foreach(i = 1:ncol(a), .packages = c("openair"), .combine = "c") %dopar% {
    mean(selectByDate(data.frame(date = rastertime, a[,i]), start = sample$date0[i], end = sample$date[i])[,2], na.rm = T)
  }
  stopCluster(cl)
})
sample0914$pppndvi200 <- pppndvi200
sample0914$pppndvi200.cut <- quant_groups(sample0914$pppndvi200, 4)

#Median for p trend
a <- data.table(sample0914)
a[, pppndvi200.cut.median := median(pppndvi200), pppndvi200.cut]
a[, pppndvi1000.cut.median := median(pppndvi1000), pppndvi1000.cut]
a[, pppndvi1500.cut.median := median(pppndvi1500), pppndvi1500.cut]
#Reorder columns
sample0914 <- a %>% 
  relocate(pppndvi400.cut, pppndvi400.cut.median,
           pppndvi1000, pppndvi1000.cut, pppndvi1000.cut.median,
           pppndvi1500, pppndvi1500.cut, pppndvi1500.cut.median,
           .after = pppndvi400) 

#PM2.5
#load("mc.rda")#Four files included
system.time({
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  pm2.5 <- foreach(i = 1:nrow(address1), .packages = c("sf", "data.table", "lubridate"), .combine = "rbind") %dopar%{
    stations <- apc.study.area[apc.study.area$date %in% seq(from = address1$date[i] - 89, to = address1$date[i], by = 1), c('station.id', 'date', 'pm2.5')]
    stations <- st_sf(merge(stations, station.point[,c('station.id', 'geometry')]))
    stations <- st_intersection(stations, st_buffer(address1[i,], dist = 50000))
    if (nrow(stations) > 0) {
      dist <- data.table(data.table(date = stations$date, pm2.5 = stations$pm2.5), dist = as.numeric(st_distance(address1[i,], stations)))
      dist[, weight := 1/(dist**2)]
      pm2.5 <- dist[, sum(pm2.5*weight, na.rm = T)/sum(weight, na.rm = T), date]
      pm2.5 <- sum(pm2.5$V1)/90
      return(data.table(data.table(address1[i, c('sid', 'date')]), pm2.5))
    }
  }
})
a <- merge(sample0914, pm2.5[,-"geometry"], all = T)

#Alternative NDVI
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
#Season 6-11 (+2)
ggplot(sample0914, aes(x = date, y = tmp)) +
  geom_point()
setDT(sample0914)
sample0914$season <- 'w'
sample0914[month(date) %in% c(2:7), season := 'c']
sample0914$season <- factor(as.character(sample0914$season))

#Season 5-10 (+2)
setDT(sample0914)
sample0914$season1 <- 'w'
sample0914[month(date) %in% c(1:6), season1 := 'c']
sample0914[,season1 := factor(season1)]

#Season 5-10
setDT(sample0914)
sample0914$season2 <- 'c'
sample0914[month(date) %in% c(5:10), season2 := 'w']
sample0914[,season2 := factor(season2)]

#Season 6-11
setDT(sample0914)
sample0914$season3 <- 'c'
sample0914[month(date) %in% c(6:11), season3 := 'w']
sample0914[,season3 := factor(season3)]

#Age
a <- as.Date(semen$date) - as.Date(semen$dob)
semen$age <- as.numeric(floor(a/365))
sample$age.bi <- quant_groups(sample$age, 2)
levels(sample0914$age.bi) <- c('l', 'h')
#Abs
semen1$abs.cut <- "(2,3)"
semen1[semen1$abs %in% c(4,5),]$abs.cut <- "(4,5)"
semen1[semen1$abs %in% c(6,7),]$abs.cut <- "(6,7)"
semen1$abs.cut <- factor(semen1$abs.cut)
#Road length (G,S,X)
road = st_read("D:/Research Projects/地图信息等/china-latest-free.shp/gis_osm_roads_free_1.shp")
road <- road[road$fclass %in% c("motorway", "trunk", "primary", "secondary", "motorway_link", "trunk_link", "primary_link", "secondary_link"),]
road <- st_bbox(c(xmin = 109, xmax = 117.3, ymax = 25.7, ymin = 20), crs = 4326) %>%
  st_as_sfc() %>% 
  st_intersection(road, .) %>% 
  st_transform(., crs = crs)
system.time({
  road.length.tmp <- st_intersection(road[, 11], st_buffer(address1[,1], dist = 200))
  road.length.tmp$road.length <- st_length(road.length.tmp)
  length200 <- aggregate(road.length ~ sid, data = road.length.tmp, function(x) as.numeric(sum(x))/1000)
})
a <- merge(sample0914, length200, all = T)
a[is.na(length200) == T, length200 := 0]
a[is.na(sample0914$length200),]$length200 <- 0
a$length200.cut <- quant_groups(a$length200, 4)
#For the first specimen
a <- sample0914 %>%
  group_by(sid) %>%
  .[order(.$date),] %>%
  filter(row_number() == 1)
address2 <- st_as_sf(data.frame(address2))
semen1$pndvi500 <- a

#Results display####
#Table 3model <- sample0914 %>%
  lmer(semen.volume ~ I(pppndvi400/IQR(.$pppndvi400)) + length400.cut + pm2.5 + season + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid),.) %>% 
  summary %>% 
  .[["coefficients"]]
I(pppndvi400/IQR(.$pppndvi400))
I(pndvi400/IQR(.$pndvi400))

model <- sample0914 %>%
  lm(normal.forms ~ I(pppndvi1500/IQR(.$pppndvi1500)) + length1500.cut + pm2.5 + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3), data = .) %>%
  summary()%>%
  .[["coefficients"]] 
estimate <- model[2, 1]
error <- model[2, 2]
estimate <- model[c(2:4), 1]#for cut
error <- model[c(2:4), 2]#for cut
(output <- paste0(round(estimate, 3), ' (', round(estimate - 1.96*error, 3), ', ', round(estimate + 1.96*error, 3), ')') )
I(ppndvi400/IQR(.$ppndvi400))

library(car) # Companion to Applied Regression # Companion to Applied Regression
vif(a)

#continuous
system.time({
  y <- c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr')
  x <- c('pppndvi400', 'pppndvi1000', 'pppndvi1500')
  z <- c('length400.cut', 'length1000.cut', 'length1500.cut')
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]], na.rm = T)) + sample0914[,z[i]] + pm2.5 + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .,) %>%
        summary(.) %>%
        .[["coefficients"]] 
      estimate <- model[2, 1]
      error <- model[2, 2]
      output <- paste0(round(estimate, 3), ' (', round(estimate - 1.96*error, 3), ', ', round(estimate + 1.96*error, 3), ')') 
      data.frame(output, row.names = paste(x[i], y[j]))
    }}
  stopCluster(cl)
})
write.table(a, file = "a", quote = F, sep = ";", row.names = T, col.names = T)
#cut
system.time({
  y <- match(c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr'), names(sample0914))
  x <- match(c('pppndvi400.cut', 'pppndvi1000.cut', 'pppndvi1500.cut'), names(sample0914))
  z <- match(c('length400.cut', 'length1000.cut', 'length1500.cut'), names(sample0914))
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a.cut <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ sample0914[,x[i]] + sample0914[,z[i]] + pm2.5 + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[["coefficients"]]
        estimate <- model[c(2:4), 1]
        error <- model[c(2:4), 2]
        output <- paste0(round(estimate, 3), ' (', round(estimate - 1.96*error, 3), ', ', round(estimate + 1.96*error, 3), ')') 
        output <- t(data.frame(output))
    }}
  a.cut <- data.frame(a.cut, row.names = paste(rep(colnames(sample0914)[x], each = length(y)), rep(colnames(sample0914)[y], times = length(x))))
  stopCluster(cl)
})
write.table(a.cut, file = "a", sep = ";", quote = F, row.names = T, col.names = T)
#linear trend
system.time({
  y <- match(c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr'), names(sample0914))
  x <- match(c('pppndvi400.cut.median', 'pppndvi1000.cut.median', 'pppndvi1500.cut.median'), names(sample0914))
  z <- match(c('length400.cut', 'length1000.cut', 'length1500.cut'), names(sample0914))
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + sample0914[,z[i]] + pm2.5 + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[["coefficients"]] 
      output <- model[2, 5]
      data.frame(output, row.names = paste(colnames(sample0914)[x[i]], colnames(sample0914)[y[j]]))
    }}
  stopCluster(cl)
})
write.table(round(a, 3), file = "a", quote = F, sep = ";", row.names = T, col.names = T)

#Table 4 Stratified
model <- sample0914 %>%
  lm(normal.forms ~ I(pppndvi1000/IQR(pppndvi1000)) + pm2.5 + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3), .,
       subset = season3 == 'w') %>%
  summary(.) %>%
  .[["coefficients"]] 
estimate <- model[2, 1]
error <- model[2, 2]
(output <- paste0(round(estimate, 3), ' (', round(estimate - 1.96*error, 3), ', ', round(estimate + 1.96*error, 3), ')') )

lrtest(lmer(pr ~ I(pppndvi400/IQR(pppndvi400)) + pm2.5 + fertility + edu + han + age.cut + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), sample0914),
       lmer(pr ~ I(pppndvi400/IQR(pppndvi400)) * season - season + pm2.5 + fertility + edu + han + age.cut + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), sample0914))

lrtest(lm(normal.forms ~ I(pppndvi400/IQR(pppndvi400)) + age.cut + pm2.5 + fertility + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3), sample0914),
       lm(normal.forms ~ I(pppndvi400/IQR(pppndvi400)) * season3 - season3 + age.cut + pm2.5 + fertility + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3), sample0914))
#continuous
system.time({
  y <- c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr')
  x <- c('pppndvi400', 'pppndvi1000', 'pppndvi1500')
  z <- c('length400.cut', 'length1000.cut', 'length1500.cut')
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + pm2.5 + fertility + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .,
             subset = season3 == 'w') %>%
        summary(.) %>%
        .[["coefficients"]] 
      estimate <- model[2, 1]
      error <- model[2, 2]
      output <- paste0(round(estimate, 3), ' (', round(estimate - 1.96*error, 3), ', ', round(estimate + 1.96*error, 3), ')') 
      data.frame(output, row.names = paste(x[i], y[j]))
    }}
  stopCluster(cl)
})
write.table(a, file = "a", quote = F, sep = ";", row.names = T, col.names = T)
#Effect modification
system.time({
  y <- c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr')
  x <- c('pppndvi400', 'pppndvi1000', 'pppndvi1500')
  z <- c('length400.cut', 'length1000.cut', 'length1500.cut')
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    foreach(j = 1:length(y), .packages = c("lmtest", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- lrtest(lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + age.cut + pm2.5 + fertility + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), sample0914,),
                      lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) * season3 - season3 + age.cut + pm2.5 + fertility + edu + han + abs.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), sample0914,))
      output <- model[2, 5]
      data.frame(output, row.names = paste(x[i], y[j]))
    }}
  stopCluster(cl)
})
write.table(round(a, 3), file = "a", quote = F, sep = ";", row.names = T, col.names = T)

#cut
system.time({
  y <- match(c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr'), names(sample0914))
  x <- match(c('pppndvi400.cut', 'pppndvi1000.cut', 'pppndvi1500.cut'), names(sample0914))
  z <- match(c('length400.cut', 'length1000.cut', 'length1500.cut'), names(sample0914))
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a.cut <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    a <- foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ sample0914[,x[i]] + sample0914[,z[i]] + pm2.5 + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[[10]]
      estimate <- model[c(2:4), 1]
      error <- model[c(2:4), 2]
      output <- paste0(round(estimate, 3), ' (', round(estimate - 1.96*error, 3), ', ', round(estimate + 1.96*error, 3), ')') 
      output <- t(data.frame(output))
    }}
  a.cut <- data.frame(a.cut, row.names = paste(rep(colnames(sample0914)[x], each = length(y)), rep(colnames(sample0914)[y], times = length(x))))
  stopCluster(cl)
})
write.table(a.cut, file = "a", sep = ";", quote = F, row.names = T, col.names = T)
#linear trend
system.time({
  y <- match(c('semen.volume', 'sperm.concentration', 'sperm.number', 'tm', 'pr'), names(sample0914))
  x <- match(c('pppndvi400.cut.median', 'pppndvi1000.cut.median', 'pppndvi1500.cut.median'), names(sample0914))
  z <- match(c('length400.cut', 'length1000.cut', 'length1500.cut'), names(sample0914))
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  a <- foreach(i = 1:length(x), .packages = c("parallel", "doParallel"), .combine = "rbind") %dopar%{
    foreach(j = 1:length(y), .packages = c("dplyr", "splines", "lubridate", "lmerTest"), .combine = "rbind") %dopar%{
      model <- sample0914 %>%
        lmer(sample0914[,y[j]] ~ I(sample0914[,x[i]]/IQR(sample0914[,x[i]])) + sample0914[,z[i]] + pm2.5 + edu + han + abs.cut + age.cut + ns(tmp, df = 3) + ns(rh, df = 3) + (1|sid), .) %>%
        summary(.) %>%
        .[[10]] 
      output <- model[2, 5]
      data.frame(output, row.names = paste(colnames(sample0914)[x[i]], colnames(sample0914)[y[j]]))
    }}
  stopCluster(cl)
})
write.table(round(a, 3), file = "a", quote = F, sep = ";", row.names = T, col.names = T)

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

