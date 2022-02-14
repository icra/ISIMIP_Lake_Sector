## checking distributions of outputs with HydroLAKES

library(foreign); library(raster)
#read HydroLAKES table
HL <- read.dbf("/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.dbf")

#read maps: Depth
depth <- raster("/home/ry4902/ISIMIP_Lake_Sector/Cosas_Rafa/Results/Depth_avg.tif")

#number of pixels
length(depth[][!is.na(depth[])])

#plot
plot(density(log10(HL$Depth_avg)), ylim=c(0,2))
lines(density(log10(depth[][!is.na(depth[])])), col="blue")

#read maps: Area and plots
area <- raster("/home/ry4902/ISIMIP_Lake_Sector/Cosas_Rafa/Results/Lake_area.tif")
plot(density(log10(HL$Lake_area)), ylim=c(0,2))
lines(density(log10(area[][!is.na(area[])])), col="blue")

#read maps: lake type and plots
ltype <- raster("/home/ry4902/ISIMIP_Lake_Sector/Cosas_Rafa/Results/Lake_type.tif")
ltype_data <- ltype[][!is.na(ltype[])]
paste0("Natural lakes numbers: from ",round(length(which(HL$Lake_type==1))/length(HL$Lake_type)*100,2),
      "% to ", round(length(which(ltype_data==1))/length(ltype_data)*100,2), "%")
paste0("Reservoirs numbers: from ",round(length(which(HL$Lake_type==2))/length(HL$Lake_type)*100,2),
       "% to ", round(length(which(ltype_data==2))/length(ltype_data)*100,2), "%")
paste0("Natural-Reservoirs numbers: from ",round(length(which(HL$Lake_type==3))/length(HL$Lake_type)*100,2),
       "% to ", round(length(which(ltype_data==3))/length(ltype_data)*100,2), "%")




