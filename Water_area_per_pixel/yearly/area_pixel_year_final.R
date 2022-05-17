library(raster)
#present day area per pixel map
frac <- raster("/home/dmercado/ISIMIP_Lake_Sector/Water_area_per_pixel/frac_areas.tif")
#all years with dam construction since 1850 (included)
years <- as.numeric(unlist(strsplit(unlist(strsplit(list.files(path="previous/",pattern="NA"), "frac_areas_NA_")), ".tif")))

for (c in 1:length(years)){
  print(paste("year:", years[c]))
  rest_temp <- frac
  if(c==length(years)){
    writeRaster(frac, paste0("final/frac_areas_NA_", years[c],".tif"), overwrite=T)
  }else{
    for (y in years[(c+1):length(years)]){
      rest_temp <- rest_temp - raster(paste0("previous/frac_areas_", y,".tif"))
    }
    m <- c(0, 0, NA,1,2,1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    rest_temp_NA <- reclassify(rest_temp, rclmat, include.lowest=TRUE)
    writeRaster(rest_temp_NA,paste0("final/frac_areas_NA_", years[c],".tif"), overwrite=T)
  }
}

