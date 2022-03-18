## checking distributions of outputs with HydroLAKES

library(foreign); library(raster)
#read HydroLAKES table
HL <- read.dbf("input/HydroLAKES_polys_v10.dbf")

#read maps: Depth
depth <- raster("output/Depth_avg.tif")
depth_median <- raster("Cosas_Rafa/Results_median/Depth_avg.tif")

#number of pixels
length(depth[][!is.na(depth[])])

#plot
pdf("output/depth_density.pdf")
plot(density(log10(HL$Depth_avg)), ylim=c(0,2))
lines(density(log10(depth[][!is.na(depth[])])), col="blue")
lines(density(log10(depth_median[][!is.na(depth_median[])])), col="red")
dev.off()

#read maps: Area and plots
area <- raster("output/Lake_area.tif")
area_median <- raster("./Cosas_Rafa/Results_median/Lake_area.tif")
area_rafa2<-  raster("./Cosas_Rafa/Results/Lake_area.tif")

#plot
pdf("output/area_density.pdf")
plot(density(log10(HL$Lake_area)), ylim=c(0,2))
lines(density(log10(area[][!is.na(area[])])), col="blue")
lines(density(log10(area_median[][!is.na(area_median[])])), col="red")
lines(density(log10(area_rafa2[][!is.na(area_rafa2[])])), col="green")

dev.off()

#read maps: lake type and plots
ltype <- raster("outputs/Lake_type.tif")
ltype_data <- ltype[][!is.na(ltype[])]
paste0("Natural lakes numbers: from ",round(length(which(HL$Lake_type==1))/length(HL$Lake_type)*100,2),
      "% to ", round(length(which(ltype_data==1))/length(ltype_data)*100,2), "%")
paste0("Reservoirs numbers: from ",round(length(which(HL$Lake_type==2))/length(HL$Lake_type)*100,2),
       "% to ", round(length(which(ltype_data==2))/length(ltype_data)*100,2), "%")
paste0("Natural-Reservoirs numbers: from ",round(length(which(HL$Lake_type==3))/length(HL$Lake_type)*100,2),
       "% to ", round(length(which(ltype_data==3))/length(ltype_data)*100,2), "%")

#number of big lakes (>0.5 degrees)
HLid_biglakes <- raster("output/HLid_biglakes.tif")
length(HLid_biglakes[][!is.na(HLid_biglakes[])])