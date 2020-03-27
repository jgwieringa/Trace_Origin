library(spatstat)
library(sp)
library(adegenet)
library(raster)

setwd("~/TE")

labo_te=read.csv("labo-ICPMS-finalconc.csv")
labo_loc=labo_te[,-4:-27]
labo_loc=labo_loc[,-1]


labo_te_pc=labo_te[,-1:-3]


te.graph=chooseCN(labo_loc,type=5,d1=4,d2=20)

te.spca1=spca(labo_te_pc,xy=labo_loc,cn=te.graph,scannf=FALSE)
te.spca1

labo_pc1=cbind(labo_loc,te.spca1$ls[1])
labo_pc1

odds=3

data=cbind(labo_loc, te.spca1$ls[,1])
###Average of all sites
te.sd=1.98 ####Calculated in arcmap for within site variation
###Highest site
#te.sd=3.89
#te.sd=2.7


basemap.cca <- read.csv("basemap.csv")
names(basemap.cca) <- c("Longitude", "Latitude","sPCA") # rename columns

data <- read.csv('te.pc1.csv', header = TRUE)

assign <- list() # create an empty list; after the assignment, the list will contain one dataframe for each individual
for(i in 1:nrow(data)){    
  prob <- numeric()
  for(j in 1:nrow(basemap.cca)){  # loops through each row of the basemap
    prob[j] <- dnorm(data$te.spca1.ls...1.[i], mean = basemap.cca$sPCA[j], te.sd) # probability that indv i originated from basemap cell j
  }
  
  rel.prob <- prob/max(prob) # calculate relative probability of originating in each raster cell
  occ <- ifelse(rel.prob < 0.5 , 0,1) # convert relative probs to "likely" or "unlikely"
  assign[[i]] <- data.frame(id = basemap.cca$sPCA, x = basemap.cca$Longitude, y = basemap.cca$Latitude, z = rel.prob, o = occ) # create data from with raster xy, rel prob and occ. for indv. i
}


###Creating figures
ext_pts=labo_loc

ext_prob=matrix(nrow=nrow(data),ncol=1)
ext_prob2=ext_prob


for(i in 1:nrow(data)){
  ####Creating assignment maps
  ind=assign[[i]]  ####Pulling out individual
  indb=ind[,-1]  ####Removing the ID column
  ind_prob=indb[,-4]
  ind_yn=indb[,-3]
  
  ###creating probability raster
  ra=rasterFromXYZ(ind_prob)
  fn=paste(i, '_prob.tif', sep = '')
  writeRaster(ra,fn,format='GTiff',overwrite=TRUE)
  
  point=matrix(c(data[i,]$Long,y=data[i,]$Lat),nrow=1,ncol=2)
  ext_prob[i]=extract(ra, point)
  
  fn=paste(i, '_prob.pdf', sep = '')
  pdf(file=fn)
  plot(ra)
  points(x=data[i,]$Long,y=data[i,]$Lat, col='black',pch=16)
  dev.off()
  
  ###Creaing yes/no raster
  ra=rasterFromXYZ(ind_yn)
  fn=paste(i, '_yn.tif', sep = '')
  writeRaster(ra,fn,format='GTiff',overwrite=TRUE)
  
  point=matrix(c(data[i,]$Long,y=data[i,]$Lat),nrow=1,ncol=2)
  ext_prob2[i]=extract(ra, point)
  
  fn=paste(i, '_yn.pdf', sep = '')
  pdf(file=fn)
  plot(ra)
  points(x=data[i,]$Long,y=data[i,]$Lat, col='black',pch=16)
  dev.off()
  
}

ext_pts_prob=cbind(ext_pts,ext_prob)
ext_pts_prob2=cbind(ext_pts,ext_prob2)
sum(ext_pts_prob2$ext_prob2)



###Testing this

test=read.csv("labo_test_assignment.csv")
pc_test=matrix(nrow=nrow(test),ncol=4)
pc_test[,1]=test[,1]
pc_test[,2]=test[,2]
pc_test[,3]=test[,3]

colnames(pc_test)=c('ID','Long','Lat','PC')
assign=list()

for(i in 1:nrow(test)){
  labo_te_pc2=rbind(labo_te_pc,test[i,4:27])
  labo_loc2=rbind(labo_loc,test[i,2:3])
  
  te.graph2=chooseCN(labo_loc2,type=5,d1=4,d2=20)
  
  te.spca2=spca(labo_te_pc2,xy=labo_loc2,cn=te.graph2,scannf=FALSE)
  
  pc_test[i,4]=te.spca2$ls[128,1]
  
}


for(i in 1:nrow(pc_test)){    
  prob <- numeric()
  for(j in 1:nrow(basemap.cca)){  # loops through each row of the basemap
    prob[j] <- dnorm(pc_test[i,4], mean = basemap.cca$sPCA[j], te.sd) # probability that indv i originated from basemap cell j
  }
  
  rel.prob <- prob/max(prob) # calculate relative probability of originating in each raster cell
  occ <- ifelse(rel.prob < 0.3 , 0,1) # convert relative probs to "likely" or "unlikely"
  assign[[i]] <- data.frame(id = basemap.cca$sPCA, x = basemap.cca$Longitude, y = basemap.cca$Latitude, z = rel.prob, o = occ) # create data from with raster xy, rel prob and occ. for indv. i
}


ext_pts=pc_test[,2:3]

ext_prob=matrix(nrow=nrow(pc_test),ncol=1)
ext_prob2=ext_prob


for(i in 1:nrow(pc_test)){
  ID=pc_test[,1]
  ####Creating assignment maps
  ind=assign[[i]]  ####Pulling out individual
  indb=ind[,-1]  ####Removing the ID column
  ind_prob=indb[,-4]
  ind_yn=indb[,-3]
  
  ###creating probability raster
  ra=rasterFromXYZ(ind_prob)
  fn=paste(ID[i], '_prob.tif', sep = '')
  writeRaster(ra,fn,format='GTiff',overwrite=TRUE)
  
  point=matrix(c(pc_test[i,2],y=pc_test[i,3]),nrow=1,ncol=2)
  ext_prob[i]=extract(ra, point)
  
  fn=paste(ID[i], '_prob.pdf', sep = '')
  pdf(file=fn)
  plot(ra)
  points(x=pc_test[i,2],y=pc_test[i,3], col='black',pch=16)
  dev.off()
  
  ###Creaing yes/no raster
  ra=rasterFromXYZ(ind_yn)
  fn=paste(ID[i], '_yn.tif', sep = '')
  writeRaster(ra,fn,format='GTiff',overwrite=TRUE)
  
  point=matrix(c(pc_test[i,2],y=pc_test[i,3]),nrow=1,ncol=2)
  ext_prob2[i]=extract(ra, point)
  
  fn=paste(ID[i], '_yn.pdf', sep = '')
  pdf(file=fn)
  plot(ra)
  points(x=pc_test[i,2],y=pc_test[i,3], col='black',pch=16)
  dev.off()
  
}


ext_pts_prob=cbind(ext_pts,ext_prob)
ext_pts_prob2=cbind(ext_pts,ext_prob2)
sum(ext_pts_prob2[,3])
