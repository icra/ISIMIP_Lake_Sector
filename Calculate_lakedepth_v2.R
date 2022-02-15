## ---------------------------
##
## Script name: Calculate_lakedepth
##
## Purpose of script: To calculate lake depth from HydroLAKES v1.0 database
##
## Authors: Daniel Mercado-Bettín and Rafael Marcé
## Institute: Catalan Institute for Water Research (ICRA)
## 
## Date Created: 2022-02-11
## 
## Email: dmercado@icra.cat, rmarce@icra.cat
##
## ---------------------------
##
## Notes: 
## more details about this code in GitHub: https://github.com/icra/ISIMIP_Lake_Sector  
## 
## ---------------------------

library(sf);library(raster);library(matrixStats)

## Function to calculate real areas of each pixel according to latitude degree (Earth distortion)
#This is needed to select the big lakes for each pixel
areakm2lat <- function(latitude, resolution=0.5, R=6371007){
  # R=6371007m is the authalic earth radius at equator
  
  height <- resolution*pi/180*R #height of the cells, same value for the whole grid
  width <- (sin((latitude+resolution/2)*pi/180)-sin((latitude-resolution/2)*pi/180))*R #cells width
  a_km2 <- width*height/1e6 #cells area depending on latitude
  return(a_km2)
}

## opening HydroLAKES
HL <- shapefile("input/HydroLAKES_polys_v10.shp")

## calculate "pseudocentroid" inside polygon  
HL_sf <- st_as_sf(HL) 
HL_cent <- st_point_on_surface(HL_sf) 

## get attribute table and add coordinates:
HL_df <- data.frame(st_drop_geometry(HL_cent), st_coordinates(HL_cent))

## filling new matrix for each longitude and latitude with the weighted median values
#empty matrix to save all weighted median map for depth
wm_matrix <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix to save all weighted median map for depth

#empty matrix to save all weighted median map for the rest of the data in HydroLAKES
wm_list <- list();c<-0; var_vector <- c(1,6:13,15:19)
for (v in 1:length(var_vector)){ #vector of data to be saved from HydroLAKES database (HL_df)
  c <- c+1
  wm_list[[c]] <- matrix(data = NA, nrow=360, ncol = 720)
} 

#empty matrix to save all position of the big lakes > 0.5 degrees
wm_biglakes <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix to save all weighted median map for depth

c_lon <- 0; 
for (lon in seq(-180, 179.5, 0.5)){ #loop in longitude
  c_lon <- c_lon+1; c_lat <- 0
  for (lat in seq(90, -89.5, -0.5)){ #loop in latitude
    c_lat <- c_lat+1
    pos_lakes <- which( (HL_df$X>lon & HL_df$X<=(lon+0.5)) & (HL_df$Y<lat & HL_df$Y>=(lat-0.5)) ) 
    if (length(pos_lakes)>0){  #only doing something when having a lake
      print(paste("start lon:", lon, "and lat:", lat))  
      max_val <- max(HL_df$Lake_area[pos_lakes]) #greater area of the lakes contained in the pixel
      if(max_val>areakm2lat(lat)){ #looking for lake areas > 0.5 degrees
        posbig_max <- which(HL_df$Lake_area[pos_lakes]==max_val)
        wm_biglakes[c_lat,c_lon] <- HL_df$Hylak_id[pos_lakes[posbig_max]] #save id of lakes with area > 0.5 degrees
      }
      wm_temp <- weightedMedian(HL_df$Depth_avg[pos_lakes],HL_df$Lake_area[pos_lakes],interpolate=F, ties="max") #weighted median calculation
      wm_matrix[c_lat,c_lon] <- wm_temp #add to final matrix
      position_subset <- which(HL_df$Depth_avg[pos_lakes]==wm_temp) #find position to write more variables
      if (length(position_subset)>1){position_subset <- which(HL_df$Lake_area[pos_lakes[position_subset]] == max(HL_df$Lake_area[pos_lakes[position_subset]]))}
      if (length(position_subset)>1){position_subset <- 1}
      for (v in 1:length(var_vector)){wm_list[[v]][c_lat,c_lon] <- HL_df[pos_lakes[position_subset],var_vector[v]]}
    }
  }
}

## converting to raster and writing data
#depth
wm_raster <- raster(wm_matrix)
extent(wm_raster) <- extent(c(-180,180,-90,90))
crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(wm_raster,"output/Depth_avg.tif", overwrite=T)
#rest of variables
for (v in 1:length(var_vector)){
  HL_name <- names(HL_df)[var_vector[v]]
  wm_raster <- raster(wm_list[[v]])
  extent(wm_raster) <- extent(c(-180,180,-90,90))
  crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  writeRaster(wm_raster, paste0("output/", HL_name,".tif"), overwrite=T)
}
#save HL id of big lakes
wm_raster <- raster(wm_biglakes)
extent(wm_raster) <- extent(c(-180,180,-90,90))
crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(wm_raster,"output/HLid_biglakes.tif", overwrite=T)

## save shapefile with the final selected lakes
HL_id <- wm_list[[1]][!is.na(wm_list[[1]])]
HL_sf_subset <- subset(HL_sf, HL_sf$Hylak_id %in% HL_id)
st_write(HL_sf_subset, "output/HL_selected.shp", append=FALSE)

## save shapefile with the initial inside-polygon-centroids of all lakes
st_write(HL_cent, "output/HL_cent.shp", append=FALSE)