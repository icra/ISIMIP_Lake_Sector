library(raster)
#present area per pixel map
frac <- raster("/home/dmercado/ISIMIP_Lake_Sector/Water_area_per_pixel/frac_areas.tif")
#all years with dam construction since 1850 (included)
years <- as.numeric(unlist(strsplit(unlist(strsplit(list.files(path="previous/",pattern="NA"), "frac_areas_NA_")), ".tif")))

for (c in 1:length(years)){
  print(paste("year:", years[c]))
  frac_temp <- raster(paste0("previous/frac_areas_", years[1],".tif"))
  if (c>1){
    for (y in years[2:c]){
      frac_temp <- frac_temp + raster(paste0("previous/frac_areas_", y,".tif"))
    }
    rest_temp <- frac - frac_temp
    m <- c(0, 0, NA,1,2,1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    rest_temp_NA <- reclassify(rest_temp, rclmat, include.lowest=TRUE)
    writeRaster(rest_temp_NA,paste0("final/frac_areas_", years[c],".tif"),overwrite=T)
  }else{
    rest_temp <- frac - frac_temp
    m <- c(0, 0, NA,1,2,1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    rest_temp_NA <- reclassify(rest_temp, rclmat, include.lowest=TRUE)
    writeRaster(rest_temp_NA,paste0("final/frac_areas_", years[c],".tif"),overwrite=T)
  }
  
}
