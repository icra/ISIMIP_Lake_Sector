library(foreign)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(sf)
library(dplyr)
library(stars)
library(viridis)


#load HydroLakes polygons
data<-read.dbf("/home/rmarce/Cloud/a. WATExR/ISIMIP/Area_depth February 2022/Hydrolakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.dbf")
#load maximum depth (Dmax) from Khazaei et al. (2022)
Dmax <-read.csv("/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv",header=T)

#checking order to make sure we do not mess with different lakes, sum must be zero 
sum(data$Hylak_id-Dmax$Hylak_id)

#adding Dmax to Hydrolakes
data$Dmax_Khazaei <- Dmax$Dmax_use_m #this is the one that paper reccomends
data$Dmax_Khazaei_cone <- Dmax$Dmax_cone_m #assuming a cone, provided in that same paper


#comparing Dmax_Khazaei vs mean depth in Hydrolakes. Some impossible situations.
#with ratios
ratio_mm <-data$Dmax_Khazaei/data$Depth_avg 
summary(ratio_mm)
hist(log10(ratio_mm))

#with differences
dif_mm <- data$Dmax_Khazaei-data$Depth_avg
summary(dif_mm)
hist(dif_mm)
boxplot(dif_mm)

#now with the cone assumption, better behaviour
dif_mm_cone <- data$Dmax_Khazaei_cone-data$Depth_avg
summary(dif_mm_cone)
hist(dif_mm_cone)

#percent of impossible values with Dmax_Khazaei
length(which(dif_mm<0))/length(dif_mm)*100

#cheking the volume development using Dmax_Khazaei (should between 0 and 3)
data$Vd  <- 3*data$Depth_avg/data$Dmax_Khazaei
summary(data$Vd)
hist(log10(data$Vd)) 
hist((data$Vd),xlim=c(0,3),breaks=400)
#nonsense values (>3) and median <1, not very realistic distribution
# 
# #volume of a cone
# pi*R^2*h/3
# #volume of a truncated cone
# pi*h*(r^2+r*R+R^2)/3
# 
# #volume
# Area*mean_depth
# pi*R^2*mean_depth
# # pi*R^2*mean_depth
# pi*R^2*h/3
# 
# mean_depth
# h/3
# # 3*mean_depth/h
#




############ Testing hypsografics in Khazaei
nc_data <- nc_open('/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_hAV_relationships.nc')
# Save the print(nc) dump to a text file
{
  sink('gimms3g_ndvi_1982-2012_metadata.txt')
  print(nc_data)
  sink()
}
#load the depth index of the hypsografics (11 depths per)
id <- ncvar_get(nc_data, "h")
(id[1:11,1:11])

#load the lake attributes:max depth, mean depth, surface area, total volume (m, m, km2, km3)
attri_nc <- ncvar_get(nc_data, "lake_attributes")
attri_nc[1:4,1:3]

#Collecting maximum volume and area, max and mean depth
A_Khazaei <- attri_nc[3,] 
V_Khazaei <- attri_nc[4,]
Dmax_Khazaei <- attri_nc[1,] 
Dmean_Khazaei <-attri_nc[2,] 

#comparing Hydrolakes and Khazaei volumes. Volume in hydrolakes is in 0.001 km3
plot(data$Vol_total*0.001,V_Khazaei,log="xy")
abline(0,1)
#comparing Hydrolakes and Khazaei areas
plot(data$Lake_area,A_Khazaei,log="xy")
abline(0,1)
#comparing Hydrolakes and Khazaei mean depth
plot(data$Depth_avg,Dmean_Khazaei,log="xy")
abline(0,1)

#difference in total volum and area, in percent 
(sum(V_Khazaei)-sum(data$Vol_total)*0.001)/(sum(data$Vol_total)*0.001)*100
(sum(A_Khazaei)-sum(data$Lake_area))/(sum(data$Lake_area))*100

#distribution of individual differences, in percent
dif_area <- (data$Lake_area-A_Khazaei)/data$Lake_area*100
hist(dif_area)

dif_volume <- (data$Vol_total*0.001-V_Khazaei)/(data$Vol_total*0.001)*100
hist(dif_volume)

dif_Dmean <- (data$Depth_avg-Dmean_Khazaei)/data$Depth_avg*100
hist(dif_Dmean)

# volume development using data from Khazaei exclusively
Vd_Khazaei  <- 3*Dmean_Khazaei/Dmax_Khazaei
hist((Vd_Khazaei))
summary(Vd_Khazaei)


#ANalyzing the results for the selected representative lakes
#recovering lakes ID from Khazaei databse
nc_lake_ID <- ncvar_get(nc_data, "lake_id")
dim(nc_lake_ID)

#loading the selected representative lakes for ISIMIP3
data_selected<-read.dbf("/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/HL_selected.dbf")
head(data_selected$Hylak_id)

#attributes of the selected representative lakes  
attri_nc_selected <- attri_nc[,data_selected$Hylak_id]
dim(attri_nc_selected)
attri_nc_selected[,1:10]

#Atttributes to variables
A_Khazaei_selected <- attri_nc_selected[3,]
V_Khazaei_selected <- attri_nc_selected[4,]
Dmax_Khazaei_selected <- attri_nc_selected[1,]
Dmean_Khazaei_selected <-attri_nc_selected[2,] 


#Volume development for selected lakes
Vd_Khazaei_selected  <- 3*Dmean_Khazaei_selected/Dmax_Khazaei_selected
hist((Vd_Khazaei_selected))

#comparing with whole-database results
summary(Vd_Khazaei_selected)
summary(Vd_Khazaei)

plot(density(Vd_Khazaei), col="red")
lines(density(Vd_Khazaei_selected))
plot(quantile(Vd_Khazaei,p=seq(0,1,0.01)),type="l")
lines(quantile(Vd_Khazaei_selected,p=seq(0,1,0.01)), col="red")

#comparing volumes and areas with Hydrolakes for selected representative lakes
plot(data$Vol_total[data_selected$Hylak_id]*0.001,V_Khazaei_selected,log="xy")
abline(0,1)
plot(data$Lake_area[data_selected$Hylak_id],A_Khazaei_selected,log="xy")
abline(0,1)

#save.image(".RData")

plot(Dmax_Khazaei_selected,Vd_Khazaei_selected,log="x")
plot(A_Khazaei_selected,Vd_Khazaei_selected,log="x")


#saving rasters for ISIMIP3

#loading lake_ID of representative lakes for reference to assign values to rasters
raster_id <- raster("/home/rmarce/ISIMIP_Lake_Sector/output/Hylak_id.tif")
raster_id_asmatrix <- as.matrix(raster_id)#to speed up

#initializing matrices that will eventually become the rasters
matrix_vd <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix
matrix_A <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix
matrix_V <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix
matrix_Dmax <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix
matrix_Dmean <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix

#loop for assigning values in the matrices 
for (i in 1:360){ 
  for (j in 1:720){ 
    #print(i);print(j)
    if (is.na(raster_id_asmatrix[i,j])) {#to avoid the loop crashing
      }else{
        #identifying which lake goes to the current position
      position <- which(data_selected$Hylak_id == raster_id_asmatrix[i,j])
      matrix_vd[i,j] <-  Vd_Khazaei_selected[position]
      matrix_A[i,j] <-  A_Khazaei_selected[position]
      matrix_V[i,j] <-  V_Khazaei_selected[position]
      matrix_Dmax[i,j] <-  Dmax_Khazaei_selected[position]
      matrix_Dmean[i,j] <-  Dmean_Khazaei_selected[position]
    }
  }
}

## converting to raster and writing data
Vd_raster <- raster(matrix_vd)
A_raster <- raster(matrix_A)
V_raster <- raster(matrix_V)
Dmax_raster <- raster(matrix_Dmax)
Dmean_raster <- raster(matrix_Dmean)

extent(Vd_raster) <- extent(c(-180,180,-90,90))
extent(A_raster) <- extent(c(-180,180,-90,90))
extent(V_raster) <- extent(c(-180,180,-90,90))
extent(Dmax_raster) <- extent(c(-180,180,-90,90))
extent(Dmean_raster) <- extent(c(-180,180,-90,90))

crs(Vd_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(A_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(V_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(Dmax_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(Dmean_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

writeRaster(Vd_raster,"./Results/Vd_raster.tif", overwrite=T)
writeRaster(A_raster,"./Results/A_raster.tif", overwrite=T)
writeRaster(V_raster,"./Results/V_raster.tif", overwrite=T)
writeRaster(Dmax_raster,"./Results/Dmax_raster.tif", overwrite=T)
writeRaster(Dmean_raster,"./Results/Dmean_raster.tif", overwrite=T)


plot(Vd_raster,col=viridis(2))
plot(log10(Dmax_raster),col=viridis(2))

aa<-(V_raster/(A_raster*Dmean_raster*0.001))
summary(aa)
hist(aa)
#looks like rounding errors

bb<-(Vd_raster/(3*Dmean_raster/Dmax_raster))
summary(bb)


#Collating hypsografics


V <- ncvar_get(nc_data, "V")
V[1:11,1:3]
A <- ncvar_get(nc_data, "A")
A[1:11,1:3]