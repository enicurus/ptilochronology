###Ptilochronology data manipulation and analysis
require(raster)
require(rgdal)
require(geosphere)
require(reshape)
require(ggplot2)
require(maptools)
require(maps)
require(mapdata)
require(dismo)

###########
####SAVE ALL MODELS AS TABLES #####

#######


###Age effects

feath<-read.csv("/Users/ryanterrill/Dropbox/Ptilochronology/Feathergrowthdatabase.csv")

feath<-feath[feath$growth!=0,]
plot(feath$growth[feath$species=="capensis"])
age.g<-feath[feath$species%in%c("capensis","gujanensis","major"),]
age.g<-age.g[age.g$age%in%c("a","j"),]
age.g$speciesName<-paste(age.g$Genus,age.g$species,sep=" ")


pdf("~/Dropbox/Ptilochronology/AgeGrowthPlots,pdf")
ggplot(age.g,aes(age,growth))+geom_boxplot(outlier.shape=NA)+facet_wrap(~speciesName)+theme_bw()+geom_jitter(alpha=.2,width=.1)
dev.off()

zc.age<-t.test(feath$growth[feath$species=="capensis"&feath$age=="a"],feath$growth[feath$species=="capensis"&feath$age=="j"])


### Age and sex for all 4 species



tm.age<-t.test(feath$growth[feath$species=="major"&feath$age=="a"],feath$growth[feath$species=="major"&feath$age=="j"])

cg.age<-t.test(feath$growth[feath$species=="gujanensis"&feath$age=="a"],feath$growth[feath$species=="gujanensis"&feath$age=="j"])



####
##Geolocation

alt<-raster("~/Desktop/alt.NewWorld.30arcsec.grd")



pdf("~/Dropbox/Ptilochronology/Piaya_cayana_map.pdf")
plot(alt,col=gray.colors(100),colNA="lightblue")
points(as.numeric(as.character(feath$longitude[feath$species=="cayana"])),as.numeric(as.character(feath$latitude[feath$species=="cayana"])),pch=16,cex=.6,col="red")
title("Piaya cayana   n=454")
dev.off()

pdf("~/Dropbox/Ptilochronology/Zonotrichia_map.pdf")
plot(alt,col=gray.colors(100),colNA="lightblue")
points(as.numeric(as.character(feath$longitude[feath$species=="capensis"])),as.numeric(as.character(feath$latitude[feath$species=="capensis"])),pch=16,cex=.6,col="red")
title("Zonotrichia capensis   n=712")
dev.off()

pdf("~/Dropbox/Ptilochronology/Taraba_map.pdf")
plot(alt,col=gray.colors(100),colNA="lightblue")
points(as.numeric(as.character(feath$longitude[feath$species=="major"])),as.numeric(as.character(feath$latitude[feath$species=="major"])),pch=16,cex=.6,col="red")
title("Taraba major  n=186")
dev.off()


pdf("~/Dropbox/Ptilochronology/Cyclarhis_map.pdf")
plot(alt,col=gray.colors(100),colNA="lightblue")
points(as.numeric(as.character(feath$longitude[feath$species=="gujanensis"])),as.numeric(as.character(feath$latitude[feath$species=="gujanensis"])),pch=16,cex=.6,col="red")
title("Cyclarhis gujanensis  n=408")
dev.off()




pdf("~/Dropbox/Ptilochronology/GrowthLocalities.pdf")
plot(alt,col=gray.colors(100),colNA="lightblue")
points(as.numeric(as.character(feath$longitude)),as.numeric(as.character(feath$latitude)),pch=16,cex=.4,col="red")

title("Feather Growth Sampling Localitities")
dev.off()

ex<- c(-130,-30,-60,40)
alt.crop<-crop(alt,extent(ex))



#######
#Make data frame of records with lat long and growth rate

feath<-feath[feath$age!="j",]
spatialGrowth<-data.frame(paste(feath$Genus,feath$species,sep="_"),feath$growth,feath$latitude,feath$longitude)
colnames(spatialGrowth)<-c("species","growth","latitude","longitude")
spatialGrowth<-spatialGrowth[complete.cases(spatialGrowth),]
 write.csv(spatialGrowth,"~/Dropbox/Ptilochronology/spatialGrowth.csv")

sg<-read.csv("~/Dropbox/Ptilochronology/spatialGrowth.csv")
sg<-sg[sg$growth!=0,]
sg<-sg[sg$species!="Onychorhynchus_mexicanus",]
sg<-sg[sg$species!="Furnarius_leucopus",]

pdf("~/Dropbox/Ptilochronology/species_histograms.pdf")
ggplot(sg,aes(growth))+geom_histogram(bins=20)+facet_wrap(~species,scales="free")
dev.off()


#plot growth by latitude by species

sg<-sg[sg$species%in%c("Cyclarhis_gujanensis","Piaya_cayana","Taraba_major","Zonotrichia_capensis"),]

sg.plot<-sg
sg.plot$species<-gsub("_"," ",sg$species)



pdf("~/Dropbox/Ptilochronology/latitude_plot.pdf")


ggplot(sg.plot,aes(y=growth,x=abs(as.numeric(as.character(latitude)))))+geom_point(fill="black",alpha=.6)+geom_smooth(method="lm")+facet_wrap(~species,scales="free")+xlab("Absolute Latitude")+ylab("Feather Growth (mm/day)")+theme_bw()+ theme(strip.background = element_rect(fill="transparent"))

dev.off()

lat.lm.pc<-lm(sg$growth[sg$species=="Piaya_cayana"]~sg$latitude[sg$species=="Piaya_cayana"])
lat.lm.tm<-lm(sg$growth[sg$species=="Taraba_major"]~sg$latitude[sg$species=="Taraba_major"])
lat.lm.cg<-lm(sg$growth[sg$species=="Cyclarhis_gujanensis"]~sg$latitude[sg$species=="Cyclarhis_gujanensis"])
lat.lm.zc<-lm(sg$growth[sg$species=="Zonotrichia_capensis"]~sg$latitude[sg$species=="Zonotrichia_capensis"])

summary(lat.lm.pc)
summary(lat.lm.tm)
summary(lat.lm.cg)
summary(lat.lm.zc)


####test significance of growth*latitude by species

cg.lm<-lm(sg$growth[sg$species=="Cyclarhis_gujanensis"]~abs(as.numeric(as.character(sg$latitude[sg$species=="Cyclarhis_gujanensis"]))))
pc.lm<-lm(sg$growth[sg$species=="Piaya_cayana"]~abs(as.numeric(as.character(sg$latitude[sg$species=="Piaya_cayana"]))))
tm.lm<-lm(sg$growth[sg$species=="Taraba_major"]~abs(as.numeric(as.character(sg$latitude[sg$species=="Taraba_major"]))))
zc.lm<-lm(sg$growth[sg$species=="Zonotrichia_capensis"]~abs(as.numeric(as.character(sg$latitude[sg$species=="Zonotrichia_capensis"]))))

latGrowMat<-matrix(ncol=2,nrow=4)
colnames(latGrowMat)<-c("adj_R2","p")
rownames(latGrowMat)<-c("Cyclarhis","Piaya","Taraba","Zonotrichia")

latGrowMat[1,1]<-summary(cg.lm)$adj.r.squared
latGrowMat[1,2]<-summary(cg.lm)$coefficients[,4][2]

latGrowMat[2,1]<-summary(pc.lm)$adj.r.squared
latGrowMat[2,2]<-summary(pc.lm)$coefficients[,4][2]

latGrowMat[3,1]<-summary(tm.lm)$adj.r.squared
latGrowMat[3,2]<-summary(tm.lm)$coefficients[,4][2]

latGrowMat[4,1]<-summary(zc.lm)$adj.r.squared
latGrowMat[4,2]<-summary(zc.lm)$coefficients[,4][2]

write.csv(latGrowMat,"~/Dropbox/Ptilochronology/Latitiude_Growth_Stats.csv")



###make some maps for sampling distributions

landsea<-readShapeSpatial("/Users/ryanterrill/Downloads/landsea_mask/landsea_mask.shp")
landsea.crop<-crop(landsea,extent(ex))
pdf("~/Dropbox/Ptilochronology/Cyclarhis_sampling.pdf")
plot(landsea.crop, main="Cyclarhis gujanensis Sampling Localities")
plot(cg.range,col="lightblue",add=TRUE)
points(cg.locs,pch=16)
dev.off()

pdf("~/Dropbox/Ptilochronology/Piaya_sampling.pdf")
plot(landsea.crop, main="Piaya cayana Sampling Localities")
plot(pc.range,col="lightblue",add=TRUE)
points(pc.locs,pch=16)
dev.off()

pdf("~/Dropbox/Ptilochronology/Taraba_sampling.pdf")
plot(landsea.crop, main="Taraba major Sampling Localities")
plot(tm.range,col="lightblue",add=TRUE)
points(tm.locs,pch=16)
dev.off()

pdf("~/Dropbox/Ptilochronology/Zonotrichia_sampling.pdf")
plot(landsea.crop, main="Zonotrichia capensis  Sampling Localities")
plot(zc.range,col="lightblue",add=TRUE)
points(zc.locs,pch=16)
dev.off()






cg.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Cyclarhis"])),as.numeric(as.character(feath$latitude[feath$Genus=="Cyclarhis"])))
pc.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Piaya"])),as.numeric(as.character(feath$latitude[feath$Genus=="Piaya"])))
tm.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Taraba"])),as.numeric(as.character(feath$latitude[feath$Genus=="Taraba"])))
zc.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Zonotrichia"])),as.numeric(as.character(feath$latitude[feath$Genus=="Zonotrichia"])))


colnames(cg.locs)<-colnames(pc.locs)<-colnames(tm.locs)<-colnames(zc.locs)<-c("latitude","longitude")
cg.locs<-cg.locs[complete.cases(cg.locs),]
pc.locs<-pc.locs[complete.cases(pc.locs),]
tm.locs<-tm.locs[complete.cases(tm.locs),]
zc.locs<-zc.locs[complete.cases(zc.locs),]

###test distribution of each species for normality

cg.norm<-shapiro.test(sg$growth[sg$species=="Cyclarhis_gujanensis"])
fr.norm<-shapiro.test(sg$growth[sg$species=="Furnarius_rufus"])
pc.norm<-shapiro.test(sg$growth[sg$species=="Piaya_cayana"])
ps.norm<-shapiro.test(sg$growth[sg$species=="Pitangus_sulphuratus"])
tm.norm<-shapiro.test(sg$growth[sg$species=="Taraba_major"])
zc.norm<-shapiro.test(sg$growth[sg$species=="Zonotrichia_capensis"&sg$growth<3.5])

#columns to add:

#mean mass
##redo this - calculate mean mass from the data

cgtb<-read.csv("/Users/ryanterrill/Dropbox/Ptilochronology/bodymass_missing.csv",na.strings="",stringsAsFactors=FALSE)
class(cgtb$Mass)="numeric"


#calculate mean mass for all 4 species


cg.meanMass<-mean(cgtb[cgtb$genus=="Cyclarhis"&is.na(cgtb$mass)==FALSE,]$mass)
tm.meanMass<-mean(cgtb[cgtb$genus=="Taraba"&is.na(cgtb$mass)==FALSE,]$mass)

pc.meanMass<-mean(feath[feath$Genus=="Piaya"&is.na(feath$mass)==FALSE,]$mass)
zc.meanMass<-mean(feath[feath$Genus=="Zonotrichia"&is.na(feath$mass)==FALSE,]$mass)

feath$meanMass[feath$species=="gujanensis"]=cg.meanMass
feath$meanMass[feath$species=="rufus"]=45
feath$meanMass[feath$species=="cayana"]=pc.meanMass
feath$meanMass[feath$species=="sulphuratus"]=61.92
feath$meanMass[feath$species=="major"]=tm.meanMass
feath$meanMass[feath$species=="capensis"]=zc.meanMass
feath$meanMass[feath$species=="coronatus"]=15.72
feath$meanMass[feath$species=="cunicularius"]=28.2
feath$meanMass[feath$species=="malancholicus"]=37.48

feath<-feath[feath$growth!=0,]
feath4<-feath[feath$Genus=="Taraba"|feath$Genus=="Cyclarhis"|feath$Genus=="Zonotrichia"|feath$Genus=="Piaya",]

pdf("~/Dropbox/Ptilochronology/meanMass_growth.pdf")
ggplot(feath4,aes(log(meanMass),log(growth)))+geom_point(alpha=.1,size=2)+geom_smooth(method="lm")+theme_bw()
dev.off()

meanMass.lm<-lm(log(feath4$growth)~log(feath4$meanMass))
summary(meanMass.lm)


###Mass

#make a new data frame with masses and growth rates





###Match Mass measurements from SC to feath database




feathT<-feath

for (i in cgtb$number){
feath[which(feath$number==i),5]=cgtb[which(cgtb$number==i),6]
}

length(feath$mass[!is.na(feath$mass)])

feath.bm<-data.frame(feath[!is.na(feath$mass),])
feath.bm<-feath.bm[feath.bm$growth>0,]
feath.bm<-feath.bm[feath.bm$Genus%in%c("Piaya","Cyclarhis","Taraba","Zonotrichia"),]

feath.bm$name=NA
feath.bm$name[feath.bm$Genus=="Piaya"]="Piaya cayana"
feath.bm$name[feath.bm$Genus=="Taraba"]="Taraba major"
feath.bm$name[feath.bm$Genus=="Cyclarhis"]="Cyclarhis gujanensis"
feath.bm$name[feath.bm$Genus=="Zonotrichia"]="Zonotrichia capensis"

feath.bm$mass<-as.numeric(gsub("not recorded","NA",feath.bm$mass))

PC.masstest<-lm(feath.bm$growth[feath.bm$Genus=="Piaya"]~feath.bm$mass[feath.bm$Genus=="Piaya"])
TM.masstest<-lm(feath.bm$growth[feath.bm$Genus=="Taraba"]~feath.bm$mass[feath.bm$Genus=="Taraba"])
CG.masstest<-lm(feath.bm$growth[feath.bm$Genus=="Cyclarhis"]~feath.bm$mass[feath.bm$Genus=="Cyclarhis"])
ZC.masstest<-lm(feath.bm$growth[feath.bm$Genus=="Zonotrichia"]~feath.bm$mass[feath.bm$Genus=="Zonotrichia"])


summary(PC.masstest)
summary(TM.masstest)
summary(CG.masstest)
summary(ZC.masstest)

pdf("~/Dropbox/Ptilochronology/Mass~growth_4sp.pdf")
ggplot(feath.bm,aes(y=log(growth),x=log(mass)))+geom_point(fill="black",alpha=.6)+geom_smooth(method="lm")+facet_wrap(~name,scales="free")+theme_bw()
dev.off()



####Model test each


nrow(feath.bm[feath.bm$Genus=="Piaya"&feath.bm$growth>0,])
nrow(feath.bm[feath.bm$Genus=="Zonotrichia",])

pc.growth.norm<-shapiro.test(feath.bm$growth[feath.bm$Genus=="Piaya"])
pc.mass.norm<-shapiro.test(feath.bm$growth[feath.bm$Genus=="Piaya"])


zc.growth.norm<-shapiro.test(feath.bm$growth[feath.bm$Genus=="Zonotrichia"])
zc.mass.norm<-shapiro.test(feath.bm$growth[feath.bm$Genus=="Zonotrichia"])


zc.gm<-lm(feath.bm$growth[feath.bm$Genus=="Zonotrichia"]~feath.bm$mass[feath.bm$Genus=="Zonotrichia"])

pc.gm<-lm(feath.bm$growth[feath.bm$Genus=="Piaya"]~feath.bm$mass[feath.bm$Genus=="Piaya"])







#distance from edge of range
#read in range maps

cg.ns<-readShapeSpatial("/Users/ryanterrill/Dropbox/Ptilochronology/Cyclarhis_gujanensis.shp") 
pc.ns<-readShapeSpatial("/Users/ryanterrill/Dropbox/Ptilochronology/Piaya_cayana.shp")
tm.ns<-readShapeSpatial("/Users/ryanterrill/Dropbox/Ptilochronology/Taraba_major.shp")
zc.ns<-readShapeSpatial("/Users/ryanterrill/Dropbox/Ptilochronology/Zonotrichia_capensis.shp")

cg.growth.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Cyclarhis"])),as.numeric(as.character(feath$latitude[feath$Genus=="Cyclarhis"])),feath$growth[feath$Genus=="Cyclarhis"])
pc.growth.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Piaya"])),as.numeric(as.character(feath$latitude[feath$Genus=="Piaya"])),feath$growth[feath$Genus=="Piaya"])
tm.growth.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Taraba"])),as.numeric(as.character(feath$latitude[feath$Genus=="Taraba"])),feath$growth[feath$Genus=="Taraba"])
zc.growth.locs<-data.frame(as.numeric(as.character(feath$longitude[feath$Genus=="Zonotrichia"])),as.numeric(as.character(feath$latitude[feath$Genus=="Zonotrichia"])),feath$growth[feath$Genus=="Zonotrichia"])

colnames(cg.growth.locs)<-colnames(pc.growth.locs)<-colnames(tm.growth.locs)<-colnames(zc.growth.locs)<-c("latitude","longitude","growth")
cg.growth.locs<-cg.growth.locs[complete.cases(cg.growth.locs),]
pc.growth.locs<-pc.growth.locs[complete.cases(pc.growth.locs),]
tm.growth.locs<-tm.growth.locs[complete.cases(tm.growth.locs),]
zc.growth.locs<-zc.growth.locs[complete.cases(zc.growth.locs),]

cg.locs<-data.frame(cg.growth.locs[,1],cg.growth.locs[,2])
pc.locs<-data.frame(pc.growth.locs[,1],pc.growth.locs[,2])
tm.locs<-data.frame(tm.growth.locs[,1],tm.growth.locs[,2])
zc.locs<-data.frame(zc.growth.locs[,1],zc.growth.locs[,2])


colnames(cg.locs)<-colnames(pc.locs)<-colnames(tm.locs)<-colnames(zc.locs)<-c("latitude","longitude")
cg.locs<-cg.locs[complete.cases(cg.locs),]
pc.locs<-pc.locs[complete.cases(pc.locs),]
tm.locs<-tm.locs[complete.cases(tm.locs),]
zc.locs<-zc.locs[complete.cases(zc.locs),]

############################################################################
#Raw distance from edge doesn't really make sense - take this out of the study ####
############################################################################

###Calculate dist from edge of range for each locality by species
### These are grayed out because they take forever - load the points fromt the csv files
#pc.dists<-dist2Line(pc.locs,pc.ns)
#cg.dists<-dist2Line(cg.locs,cg.ns)
#tm.dists<-dist2Line(tm.locs,tm.ns)
#zc.dists<-dist2Line(zc.locs,zc.ns)


write.csv(pc.dists,"~/Dropbox/Ptilochronology/Piaya_dits2edge.csv")
write.csv(cg.dists,"~/Dropbox/Ptilochronology/Cyclarhis_dits2edge.csv")
write.csv(tm.dists,"~/Dropbox/Ptilochronology/Taraba_dits2edge.csv")
write.csv(zc.dists,"~/Dropbox/Ptilochronology/Zonotrichia_dits2edge.csv")


###Compare dist to growth

cg.gd<-data.frame(cg.dists[,1],cg.growth.locs$growth)
colnames(cg.gd)<-c("dist2edge","growth")
cg.gd<-cg.gd[cg.gd$growth!=0,]
cg.gd.lm<-lm(cg.gd[,1]~cg.gd[,2])
summary(cg.gd.lm)

pdf("~/Dropbox/Ptilochronology/Cyclarhis_growth_x_dist2edge.pdf")
ggplot(cg.gd,aes(growth,dist2edge))+geom_point()+geom_smooth(method="lm")+theme_bw()
dev.off()

pc.gd<-data.frame(pc.dists[,1],pc.growth.locs$growth)
colnames(pc.gd)<-c("dist2edge","growth")
pc.gd<-pc.gd[pc.gd$growth!=0,]
pc.gd.lm<-lm(pc.gd[,1]~pc.gd[,2])
summary(pc.gd.lm)




pdf("~/Dropbox/Ptilochronology/Piaya_growth_x_dist2edge.pdf")
ggplot(pc.gd,aes(growth,dist2edge))+geom_point()+geom_smooth(method="lm")+theme_bw()
dev.off()

tm.gd<-data.frame(tm.dists[,1],tm.growth.locs$growth)
colnames(tm.gd)<-c("dist2edge","growth")
tm.gd<-tm.gd[tm.gd$growth!=0,]
tm.gd.lm<-lm(tm.gd[,1]~tm.gd[,2])
summary(tm.gd.lm)

pdf("~/Dropbox/Ptilochronology/Taraba_growth_x_dist2edge.pdf")
ggplot(tm.gd,aes(growth,dist2edge))+geom_point()+geom_smooth(method="lm")+theme_bw()
dev.off()


zc.gd<-data.frame(zc.dists[,1],zc.growth.locs$growth)
colnames(zc.gd)<-c("dist2edge","growth")
zc.gd<-zc.gd[zc.gd$growth!=0,]
zc.gd.lm<-lm(zc.gd[,1]~zc.gd[,2])
summary(zc.gd.lm)

pdf("~/Dropbox/Ptilochronology/Zonotrichia_growth_x_dist2edge.pdf")
ggplot(zc.gd,aes(growth,dist2edge))+geom_point()+geom_smooth(method="lm")+theme_bw()
dev.off()

###compare SD in Corrientes to whole data set for CG and ZC

cg.corrientes<-feath[feath$species=="gujanensis",]
cg.corrientes<-cg.corrientes$growth[1:66]
sd(cg.corrientes)

cg.all<-feath$growth[feath$species=="gujanensis"]

zc.corrientes<-feath[feath$species=="capensis",]
zc.corrientes<-zc.corrientes$growth[4:78]
sd(zc.corrientes)

zc.all<-feath$growth[feath$species=="capensis"]



#elevation

cg.sp<-SpatialPointsDataFrame(data=cg.locs,coords=cg.growth.locs[,c("latitude","longitude")])
cg.el<-extract(alt.crop,cg.sp)

pc.sp<-SpatialPointsDataFrame(data=pc.locs,coords=pc.growth.locs[,c("latitude","longitude")])
pc.el<-extract(alt.crop,pc.sp)

tm.sp<-SpatialPointsDataFrame(data=tm.locs,coords=tm.growth.locs[,c("latitude","longitude")])
tm.el<-extract(alt.crop,tm.sp)

zc.sp<-SpatialPointsDataFrame(data=zc.locs,coords=zc.growth.locs[,c("latitude","longitude")])
zc.el<-extract(alt.crop,zc.sp)

cg.el.g<-data.frame(cg.el,cg.growth.locs$growth);cg.el.g<-cg.el.g[cg.el.g[,2]!=0,]
cg.el.g.lm<-lm(cg.el.g[,1]~cg.el.g[,2]);summary(cg.el.g.lm)

pdf("~/Dropbox/Ptilochronology/Cyclarhis_elevation_Growth.pdf")
plot(cg.el.g,main="Cyclarhis gujanensis elevation x growth rate")
dev.off()

pc.el.g<-data.frame(pc.el,pc.growth.locs$growth);pc.el.g<-pc.el.g[pc.el.g[,2]!=0,]
pc.el.g.lm<-lm(pc.el.g[,1]~pc.el.g[,2]);summary(pc.el.g.lm)


pdf("~/Dropbox/Ptilochronology/Piaya_elevation_Growth.pdf")
plot(pc.el.g,main="Piaya cayana elevation x growth rate")
dev.off()

tm.el.g<-data.frame(tm.el,tm.growth.locs$growth);tm.el.g<-tm.el.g[tm.el.g[,2]!=0,]
tm.el.g.lm<-lm(tm.el.g[,1]~tm.el.g[,2]);summary(tm.el.g.lm)

pdf("~/Dropbox/Ptilochronology/Taraba major_Growth.pdf")
plot(tm.el.g,main="Taraba major elevation x growth rate")
dev.off()

zc.el.g<-data.frame(zc.el,zc.growth.locs$growth);zc.el.g<-zc.el.g[zc.el.g[,2]!=0,]
zc.el.g.lm<-lm(zc.el.g[,1]~zc.el.g[,2]);summary(zc.el.g.lm)

pdf("~/Dropbox/Ptilochronology/Zontrichia_Growth.pdf")
plot(zc.el.g,main="Zonotrichia capensis elevation x growth rate")
dev.off()


#precip, temp from worldclim

BClim<-getData("worldclim",var="bio",res=2.5);BClim<-crop(BClim,extent(ex))
tmin<-getData("worldclim",var="tmin",res=2.5);tmin<-crop(tmin,extent(ex))
tmax<-getData("worldclim",var="tmax",res=2.5);tmax<-crop(tmax,extent(ex))
tavg<-getData("worldclim",var="tmean",res=2.5);tavg<-crop(tavg,extent(ex))


#BioClim vars - Cyclarhis
cg.bc<-extract(BClim,cg.sp)
cg.bc.g<-cbind(cg.bc,cg.growth.locs$growth)
cg.bc.g<-data.frame(cg.bc.g)
colnames(cg.bc.g)[20]<-"growth"
cg.bc.g<-cg.bc.g[cg.bc.g$growth!=0,]
cg.bc.melt<-melt(cg.bc.g,id="growth")

pdf("~/Dropbox/Ptilochronology/Cyclarhis_bioclim_models.pdf")
ggplot(cg.bc.melt,aes(y=growth,x=value))+geom_point()+facet_wrap(~variable,scales="free")+theme_bw()+ggtitle("Cyclarhis gujanensis growth vs Bioclim variables")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

summary(lm(cg.bc.g$growth~cg.bc.g[,1]))$coefficients[2,4]
summary(lm(cg.bc.g$growth~cg.bc.g[,1]))$r.squared

cg.bio.mat<-matrix(nrow=19,ncol=2)
for(i in 1:19){
	cg.bio.mat[i,1]<-summary(lm(cg.bc.g$growth~cg.bc.g[,i]))$coefficients[2,4]
	cg.bio.mat[i,2]<-summary(lm(cg.bc.g$growth~cg.bc.g[,i]))$r.squared
	}
	
colnames(cg.bio.mat)<-c("p_value","R2")	
	
write.csv(cg.bio.mat,"~/Dropbox/Ptilochronology/CycGuj_bioCim_models.csv")	
	
#BioClim vars - Piaya
pc.bc<-extract(BClim,pc.sp)
pc.bc.g<-cbind(pc.bc,pc.growth.locs$growth)
pc.bc.g<-data.frame(pc.bc.g)
colnames(pc.bc.g)[20]<-"growth"
pc.bc.g<-pc.bc.g[pc.bc.g$growth!=0,]
pc.bc.melt<-melt(pc.bc.g,id="growth")

pdf("~/Dropbox/Ptilochronology/Piaya_bioclim_models.pdf")
ggplot(pc.bc.melt,aes(y=growth,x=value))+geom_point()+facet_wrap(~variable,scales="free")+theme_bw()+ggtitle("Cyclarhis gujanensis growth vs Bioclim variables")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

summary(lm(pc.bc.g$growth~pc.bc.g[,1]))$coefficients[2,4]
summary(lm(pc.bc.g$growth~pc.bc.g[,1]))$r.squared

pc.bio.mat<-matrix(nrow=19,ncol=2)
for(i in 1:19){
	pc.bio.mat[i,1]<-summary(lm(pc.bc.g$growth~pc.bc.g[,i]))$coefficients[2,4]
	pc.bio.mat[i,2]<-summary(lm(pc.bc.g$growth~pc.bc.g[,i]))$r.squared
	}
pc.bio.mat
	
colnames(pc.bio.mat)<-c("p_value","R2")	
write.csv(pc.bio.mat,"~/Dropbox/Ptilochronology/Piaya_bioCim_models.csv")	
	
	
#BioClim vars - Taraba
tm.bc<-extract(BClim,tm.sp)
tm.bc.g<-cbind(tm.bc,tm.growth.locs$growth)
tm.bc.g<-data.frame(tm.bc.g)
colnames(tm.bc.g)[20]<-"growth"
tm.bc.g<-tm.bc.g[tm.bc.g$growth!=0,]
tm.bc.melt<-melt(tm.bc.g,id="growth")

pdf("~/Dropbox/Ptilochronology/Taraba_bioclim_models.pdf")
ggplot(tm.bc.melt,aes(y=growth,x=value))+geom_point()+facet_wrap(~variable,scales="free")+theme_bw()+ggtitle("Cyclarhis gujanensis growth vs Bioclim variables")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

summary(lm(tm.bc.g$growth~tm.bc.g[,1]))$coefficients[2,4]
summary(lm(tm.bc.g$growth~tm.bc.g[,1]))$r.squared

tm.bio.mat<-matrix(nrow=19,ncol=2)
for(i in 1:19){
	tm.bio.mat[i,1]<-summary(lm(tm.bc.g$growth~tm.bc.g[,i]))$coefficients[2,4]
	tm.bio.mat[i,2]<-summary(lm(tm.bc.g$growth~tm.bc.g[,i]))$r.squared
	}
tm.bio.mat
	
colnames(tm.bio.mat)<-c("p_value","R2")	
write.csv(tm.bio.mat,"~/Dropbox/Ptilochronology/Taraba_bioCim_models.csv")	
	
	
#Bioclimvars - Zonotrichia
zc.bc<-extract(BClim,zc.sp)
zc.bc.g<-cbind(zc.bc,zc.growth.locs$growth)
zc.bc.g<-data.frame(zc.bc.g)
colnames(zc.bc.g)[20]<-"growth"
zc.bc.g<-zc.bc.g[zc.bc.g$growth!=0,]
zc.bc.melt<-melt(zc.bc.g,id="growth")

pdf("~/Dropbox/Ptilochronology/Zonotrichia_bioclim_models.pdf")
ggplot(zc.bc.melt,aes(y=growth,x=value))+geom_point()+facet_wrap(~variable,scales="free")+theme_bw()+ggtitle("Cyclarhis gujanensis growth vs Bioclim variables")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

summary(lm(zc.bc.g$growth~zc.bc.g[,1]))$coefficients[2,4]
summary(lm(zc.bc.g$growth~zc.bc.g[,1]))$r.squared

zc.bio.mat<-matrix(nrow=19,ncol=2)
for(i in 1:19){
	zc.bio.mat[i,1]<-summary(lm(zc.bc.g$growth~zc.bc.g[,i]))$coefficients[2,4]
	zc.bio.mat[i,2]<-summary(lm(zc.bc.g$growth~zc.bc.g[,i]))$r.squared
	}
zc.bio.mat
	
colnames(zc.bio.mat)<-c("p_value","R2")	
write.csv(zc.bio.mat,"~/Dropbox/Ptilochronology/Zonotrichia_bioCim_models.csv")	
	
		
	
		

#habitat suitability


##Species distribution models - code for making models in /Dropbox/Ptilochronology/Ptilo_SDM.R 



#####################################################################################

#Cyclarhis
#load ebird data


##save the model as a raster
cg.pred_me<-raster("~/Dropbox/Ptilochronology/Cyclarhis_sdm")


#plot the model
pdf("~/Dropbox/Ptilochronology/Cyclarhis_sdm.pdf")
plot(cg.pred_me,main="Cyclarhis gujanensis Distribution Model")
dev.off()

#plot the model with sampling localities
pdf("~/Dropbox/Ptilochronology/Cyclarhis_sdm_sampling_locs.pdf")
plot(cg.pred_me,main="Cyclarhis gujanensis Distribution Model with \n Feather Growth Sampling Localities")
points(cg.locs,pch=16)
dev.off()


#extract climate suitability at specimen localities

cg.g.sdm<-extract(cg.pred_me,cg.locs)

cg.sdm.growth<-data.frame(cg.g.sdm,cg.growth.locs$growth);colnames(cg.sdm.growth)<-c("climSuit","growth")
cg.sdm.growth<-cg.sdm.growth[cg.sdm.growth$growth!=0,]

lm.cg.sdm.growth<-lm(cg.sdm.growth$climSuit~cg.sdm.growth$growth)
summary(lm.cg.sdm.growth)


pdf("~/Dropbox/Ptilochronology/Cyclarhis_growth_climSuitability.pdf")
ggplot(cg.sdm.growth,aes(y=growth,x=climSuit))+geom_point()+theme_bw()+geom_smooth(method="lm")+ggtitle("Cyclarhis gujanensis - Feather Growth Rate vs Climate Suitability")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

#####################################################################################

#####################################################################################

#Piaya
#load ebird data


##read in the SDM
pc.pred_me<-raster("~/Dropbox/Ptilochronology/Piaya_sdm")


#plot the model
pdf("~/Dropbox/Ptilochronology/Piaya_sdm.pdf")
plot(pc.pred_me,main="Piaya cayana Distribution Model")
dev.off()

#plot the model with samplinf localities
pdf("~/Dropbox/Ptilochronology/Piaya_sdm_sampling_locs.pdf")
plot(pc.pred_me,main="Piaya cayana Distribution Model with \n Feather Growth Sampling Localities")
points(pc.locs,pch=16)
dev.off()


#extract climate suitability at specimen localities

pc.g.sdm<-extract(pc.pred_me,pc.locs)

pc.sdm.growth<-data.frame(pc.g.sdm,pc.growth.locs$growth);colnames(pc.sdm.growth)<-c("climSuit","growth")
pc.sdm.growth<-pc.sdm.growth[pc.sdm.growth$growth!=0,]

lm.pc.sdm.growth<-lm(pc.sdm.growth$climSuit~pc.sdm.growth$growth)
summary(lm.pc.sdm.growth)

pdf("~/Dropbox/Ptilochronology/Piaya_growth_climSuitability.pdf")
ggplot(pc.sdm.growth,aes(y=growth,x=climSuit))+geom_point()+theme_bw()+geom_smooth(method="lm")+ggtitle("Piaya gujanensis - Feather Growth Rate vs Climate Suitability")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

#####################################################################################


#####################################################################################

#Taraba
#load ebird data


##read the SDM
tm.pred_me<-raster("~/Dropbox/Ptilochronology/Taraba_sdm")


#plot the model
pdf("~/Dropbox/Ptilochronology/Taraba_sdm.pdf")
plot(tm.pred_me,main="Taraba major Distribution Model")
dev.off()

#plot the model with samplinf localities
pdf("~/Dropbox/Ptilochronology/Taraba_sdm_sampling_locs.pdf")
plot(tm.pred_me,main="Taraba major Distribution Model with \n Feather Growth Sampling Localities")
points(tm.locs,pch=16)
dev.off()


#extract climate suitability at specimen localities

tm.g.sdm<-extract(tm.pred_me,tm.locs)

tm.sdm.growth<-data.frame(tm.g.sdm,tm.growth.locs$growth);colnames(tm.sdm.growth)<-c("climSuit","growth")
tm.sdm.growth<-tm.sdm.growth[tm.sdm.growth$growth!=0,]

lm.tm.sdm.growth<-lm(tm.sdm.growth$climSuit~tm.sdm.growth$growth)
summary(lm.tm.sdm.growth)

pdf("~/Dropbox/Ptilochronology/Taraba_growth_climSuitability.pdf")
ggplot(tm.sdm.growth,aes(y=growth,x=climSuit))+geom_point()+theme_bw()+geom_smooth(method="lm")+ggtitle("Taraba major - Feather Growth Rate vs Climate Suitability")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

#####################################################################################

#####################################################################################

#Zonotrichia
#load SDM

zc.pred_me<-raster("~/Dropbox/Ptilochronology/Zonotrichia_sdm")


#plot the model
pdf("~/Dropbox/Ptilochronology/Zonotrichia_sdm.pdf")
plot(zc.pred_me,main="Zonotrichia capensis Distribution Model")
dev.off()

#plot the model with sampling localities
pdf("~/Dropbox/Ptilochronology/Zonotrichia_sdm_sampling_locs.pdf")
plot(zc.pred_me,main="Zonotrichia capensis Distribution Model with \n Feather Growth Sampling Localities")
points(zc.locs,pch=16)
dev.off()


#extract climate suitability at specimen localities

zc.g.sdm<-extract(zc.pred_me,zc.locs)

zc.sdm.growth<-data.frame(zc.g.sdm,zc.growth.locs$growth);colnames(zc.sdm.growth)<-c("climSuit","growth")
zc.sdm.growth<-zc.sdm.growth[zc.sdm.growth$growth!=0,]

lm.zc.sdm.growth<-lm(zc.sdm.growth$climSuit~zc.sdm.growth$growth)
summary(lm.zc.sdm.growth)


pdf("~/Dropbox/Ptilochronology/Zonotrichia_growth_climSuitability.pdf")
ggplot(zc.sdm.growth,aes(y=growth,x=climSuit))+geom_point()+theme_bw()+geom_smooth(method="lm")+ggtitle("Zonotrichia capensis - Feather Growth Rate vs Climate Suitability")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

#####################################################################################





grow<-feath$growth[feath$growth!=0]
