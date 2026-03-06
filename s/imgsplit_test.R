library(imager)
library(stringr)
library(data.table)
# setwd('c:/users/habcam.local/Documents/')
imagelist1 <- read.csv('../data/raw/groundtruth2224.csv')
names(imagelist1) <- c("X1",'imagename')
imagelist1$year <- as.numeric(substr(imagelist1$imagename,8,11))
imagelist <- subset(imagelist1,year> 2021)
# setwd('d:/')
for (i in 1:length(imagelist$imagename)){
  thisimage <- load.image(paste0('images2224/',imagelist$imagename[i]))
  splitimagelist <- imsplit(thisimage,'x',2)
  ifelse(imagelist$year==2024,side<- 2,side<- 1) #side=1 left, side=2 right
  save.image(splitimagelist[[side]],paste0('split2224v2/',imagelist$imagename[i]))
}