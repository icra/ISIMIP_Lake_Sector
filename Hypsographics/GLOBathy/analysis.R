library(foreign)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

#load HydroLakes
data<-read.dbf("/home/rmarce/Cloud/a. WATExR/ISIMIP/Area_depth February 2022/Hydrolakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.dbf")
#load Dmax from Khazaei et al. (2022)
Dmax <-read.csv("/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv",header=T)

#checking order. sum must be zero 
sum(data$Hylak_id-Dmax$Hylak_id)

#adding Dmax to Hydrolakes
data$Dmax_Khazaei <- Dmax$Dmax_use_m
data$Dmax_Khazaei_cone <- Dmax$Dmax_cone_m


#cheking ratios Dmax vs Dmean
ratio_mm <-data$Dmax_Khazaei/data$Depth_avg 
summary(ratio_mm)
hist(log10(ratio_mm))
boxplot(ratio_mm)
abline(0,1)

dif_mm <- data$Dmax_Khazaei-data$Depth_avg
summary(dif_mm)
hist(dif_mm)
boxplot(dif_mm)

dif_mm_cone <- data$Dmax_Khazaei_cone-data$Depth_avg
summary(dif_mm_cone)
hist(dif_mm_cone)

length(which(dif_mm<0))/length(dif_mm)*100
hist(log10(data$Lake_area[which(dif_mm<0)]))
summary(data$Lake_area[which(dif_mm<0)])

which(data$Lake_area==16717.890) 
data[12,]

#checking outlier
which(ratio_mm>200)
data[8278,]

hist(log10(ratio_mm))
plot(data$Dmax_Khazaei,data$Depth_avg,log="xy")

data$Vd  <- 3*data$Depth_avg/data$Dmax_Khazaei
hist(log10(data$Vd))
hist((data$Vd),xlim=c(0,3),breaks=400)

log10(3)
length(which(data$Vd>3))/length(data$Vd)*100
quantile(data$Vd,0)
boxplot((data$Vd))
median(data$Vd)


r <- hist(log10(data$Vd))
plot(r)
# or alternatively:
barplot(r$counts, log="y", col="white", names.arg=r$breaks[-1])
summary(data$Vd)
#volume of a cone
pi*R^2*h/3
#volume of a truncated cone
pi*h*(r^2+r*R+R^2)/3

#volume
Area*mean_depth
pi*R^2*mean_depth

pi*R^2*mean_depth
pi*R^2*h/3

mean_depth
h/3

3*mean_depth/h





############ Testing hypsografics
nc_data <- nc_open('/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_hAV_relationships.nc')
# Save the print(nc) dump to a text file
{
  sink('gimms3g_ndvi_1982-2012_metadata.txt')
  print(nc_data)
  sink()
}
id <- ncvar_get(nc_data, "h")
(id[1:11,1:11])
ncatt_get

Dmax_nc <- ncvar_get(nc_data, "lake_attributes")
Dmax_nc[1:4,1:3]


V <- ncvar_get(nc_data, "V")
V[1:11,1:3]
A <- ncvar_get(nc_data, "A")
A[1:11,1:3]

V[11,1:3]/A[11,1:3]*1000

A_Khazaei <- A[11,]
V_Khazaei <- V[11,]*1000
Dmax_Khazaei <- Dmax_nc[1,]
Dmean_Khazaei <-Dmax_nc[2,] 


plot(data$Vol_total,V_Khazaei,log="xy")
abline(1,1)
plot(data$Lake_area,A_Khazaei,log="xy")
abline(1,1)
plot(data$Depth_avg,Dmean_Khazaei,log="xy")
abline(1,1)
sum(data$Vol_total)
sum(V_Khazaei)
sum(data$Lake_area)
sum(A_Khazaei)


dif_area <- (data$Lake_area-A_Khazaei)/data$Lake_area*100
hist(dif_area)
which(dif_area==0)

dif_volume <- (data$Vol_total-V_Khazaei)/data$Vol_total*100
hist(dif_volume)
which(dif_volume==0)

dif_Dmean <- (data$Depth_avg-Dmean_Khazaei)/data$Depth_avg*100
hist(dif_Dmean)
which(dif_Dmean==0)

hist((V_Khazaei/A_Khazaei-Dmean_Khazaei)/Dmean_Khazaei*100)
head(Dmean_Khazaei)

Vd_Khazaei  <- 3*Dmean_Khazaei/Dmax_Khazaei
hist((Vd_Khazaei))
summary(Vd_Khazaei)

#save.image(".RData")

