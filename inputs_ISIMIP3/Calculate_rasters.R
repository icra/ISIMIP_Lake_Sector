library(rgdal);library(raster)

HL_id <- raster("Hylak_id.tif")
G_areas <- read.csv("/home/rmarce/ISIMIP_Lake_Sector/Hypsographics/GLOBathy/Results/GLOBATHY_hypso_areas_representative.csv", header=F)
G_level <- read.csv("/home/rmarce/ISIMIP_Lake_Sector/Hypsographics/GLOBathy/Results/GLOBATHY_hypso_levels_representative.csv", header=F)
G_volume <- read.csv("/home/rmarce/ISIMIP_Lake_Sector/Hypsographics/GLOBathy/Results/GLOBATHY_hypso_volumes_representative.csv", header=F)
G_id <- read.csv("/home/rmarce/ISIMIP_Lake_Sector/Hypsographics/GLOBathy/Results/GLOBATHY_lakes_ID_representative.csv", header=F)

for (level in 1:11){
  area_level <- HL_id
  level_level <- HL_id
  volume_level <- HL_id
  
  area_level_index <- HL_id
  level_level_index <- HL_id
  volume_level_index <- HL_id
  
  for (pixel in 1:dim(G_id)[1]){
    area_level[][(area_level_index[]==G_id[pixel,1])] <- G_areas[level,pixel]
    level_level[][(level_level_index[]==G_id[pixel,1])] <- G_level[level,pixel]
    volume_level[][(volume_level_index[]==G_id[pixel,1])] <- G_volume[level,pixel]
  }
  writeRaster(area_level, paste0("rasters_hypsographic/area_level",level,".tif"), overwrite=T)
  writeRaster(level_level, paste0("rasters_hypsographic/level_level",level,".tif"), overwrite=T)
  writeRaster(volume_level, paste0("rasters_hypsographic/volume_level",level,".tif"), overwrite=T)
  print(paste0("termina level", level))
}


G_fitareas <- read.csv("/home/rmarce/ISIMIP_Lake_Sector/Hypsographics/GLOBathy/Results/GLOBATHY_hypso_fitArea_representative.csv", header=F)
G_fitvolume <- read.csv("/home/rmarce/ISIMIP_Lake_Sector/Hypsographics/GLOBathy/Results/GLOBATHY_hypso_fitVolume_representative.csv", header=F)

for (level in 1:3){
  fitarea_level <- HL_id
  fitvolume_level <- HL_id
  
  fitarea_level_index <- HL_id
  fitvolume_level_index <- HL_id
  
  
  for (pixel in 1:dim(G_id)[1]){
    fitarea_level[][(fitarea_level_index[]==G_id[pixel,1])] <- G_fitareas[level,pixel]
    fitvolume_level[][(fitvolume_level_index[]==G_id[pixel,1])] <- G_fitvolume[level,pixel]
  }
  writeRaster(area_level, paste0("rasters_hypsographic/fitarea_level",level,".tif"))
  writeRaster(volume_level, paste0("rasters_hypsographic/fitvolume_level",level,".tif"))
  print(paste0("termina fit level",	level))
}

