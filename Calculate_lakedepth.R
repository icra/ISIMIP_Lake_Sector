## ---------------------------
##
## Script name: Calculate_lakedepth
##
## Purpose of script: To calculate lake depth from HydroLAKES v1.0 database
##
## Author: Daniel Mercado-Bettín and Rafael Marcé
## Institute: Catalan Institute for Water Research (ICRA)
## 
## Date Created: 2022-02-10
## 
## Email: dmercado@icra.cat
##
## ---------------------------
##
## Notes: 
## more details about this code in GitHub: https://github.com/icra/ISIMIP_Lake_Sector  
##
## ---------------------------

library(sf);library(raster)

## weighted median function, modified from https://github.com/spatstat/spatstat.geom/blob/main/R/weightedStats.R
weighted_median <- function(x, w){
  probs <- 0.5 #to calculate median, one of the quantiles
  
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x_o <- x[oo]
  w_o <- w[oo]
  Fx <- cumsum(w_o)/sum(w_o)

  #'
  #this method is apply when having more than 2 points
  out <- approx(Fx, x_o, xout=probs, ties="ordered", rule=2,
                method="linear")
  
  closest_pos <- which(abs(out$y-x)==min(abs(out$y-x))) #position of the closest value
  if (length(closest_pos)>1){
    closest_pos_area <- which(w[closest_pos]==max(w[closest_pos])) #lake with greater area selected, in case there are more than one depths with the same value
    closest_pos <- closest_pos[closest_pos_area][1]
  }
  result <- c(x[closest_pos], closest_pos) #saving weighted median depth and its position in initial 'x' vector
  return(result)
}

## opening HydroLAKES
HL <- shapefile("/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/test_weightedMean/HL_test.shp")

## calculate "pseudocentroid" inside polygon  
HL_sf <- st_as_sf(HL) 
HL_cent <- st_point_on_surface(HL_sf) 

## get attribute table and add coordinates:
HL_df <- data.frame(st_drop_geometry(HL_cent), st_coordinates(HL_cent))

## filling new matrix for each longitude and latitude with the weighted median values
wm_matrix <- matrix(data = NA, nrow=360, ncol = 720) #empty matrix to save all weighted median map for depth

#empty matrix to save all weighted median map for the rest of the data in HydroLAKES
wm_list <- list();c<-0; var_vector <- c(1,6:19)
for (v in 1:length(var_vector)){ #vector of data to be saved from HydroLAKES database (HL_df)
  c <- c+1
  wm_list[[c]] <- matrix(data = NA, nrow=360, ncol = 720)
} 
c_lon <- 0; 
for (lon in seq(-180, 179.5, 0.5)){ #loop in longitude
  c_lon <- c_lon+1; c_lat <- 0
  for (lat in seq(90, -89.5, -0.5)){ #loop in latitude
    c_lat <- c_lat+1
    pos_lakes <- which( (HL_df$X>lon & HL_df$X<=(lon+0.5)) & (HL_df$Y<lat & HL_df$Y>=(lat-0.5)) ) 
    #print(pos_lakes)
    #columns to save 1:21
    if (length(pos_lakes)>0){
      max_val <- max(HL_df$Lake_area[pos_lakes])
      if(max_val>3080){
        if (length(pos_lakes)>1){
          print("hay dos lagos grandes OJO")
          stop("hay dos lagos grandes OJO")
        }
        wm_matrix[c_lat,c_lon] <- HL_df$Depth_avg[pos_lakes]
        for (v in 1:length(var_vector)){wm_list[[v]][c_lat,c_lon] <- HL_df[pos_lakes,var_vector[v]]}
      }else{
        #print(pos_lakes)
        if (length(pos_lakes)==1){
          wm_temp <- HL_df$Depth_avg[pos_lakes]
          wm_matrix[c_lat,c_lon] <- wm_temp
          for (v in 1:length(var_vector)){wm_list[[v]][c_lat,c_lon] <- HL_df[pos_lakes,var_vector[v]]}
          #print(paste("1111111111111111ésta es la wm:", wm_temp))
          #print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
          #print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
        }else if(length(pos_lakes)==2){
          #print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
          #print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
          pos2_max <- which(HL_df$Lake_area[pos_lakes]==max(HL_df$Lake_area[pos_lakes]))
          pos_lakes <- pos_lakes[pos2_max]
          wm_temp <- HL_df$Depth_avg[pos_lakes]
          wm_matrix[c_lat,c_lon] <- wm_temp
          for (v in 1:length(var_vector)){wm_list[[v]][c_lat,c_lon] <- HL_df[pos_lakes,var_vector[v]]}
          #print(paste("222222222222ésta es la wm:", wm_temp))
        }else{
          wm_temp <- weighted_median(x=HL_df$Depth_avg[pos_lakes], w=HL_df$Lake_area[pos_lakes])
          wm_matrix[c_lat,c_lon] <- wm_temp[1]
          for (v in 1:length(var_vector)){wm_list[[v]][c_lat,c_lon] <- HL_df[pos_lakes[wm_temp[2]],var_vector[v]]}
          #print(paste("ésta es la wm:", wm_temp[1], "con esta posición", wm_temp[2]))
          #print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
          #print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
        }
      }
    }
  }
}

## converting to raster and writing data
#depth
wm_raster <- raster(wm_matrix)
extent(wm_raster) <- extent(c(-180,180,-90,90))
crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(wm_raster,"/home/ry4902/ISIMIP_Lake_Sector/test.tif", overwrite=T)
#rest of variables
for (v in 1:length(var_vector)){
  HL_name <- names(HL_df)[var_vector[v]]
  wm_raster <- raster(wm_list[[v]])
  extent(wm_raster) <- extent(c(-180,180,-90,90))
  crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  writeRaster(wm_raster, paste0("/home/ry4902/ISIMIP_Lake_Sector/test",HL_name,".tif"), overwrite=T)
}

## save shapefile with the final selected lakes
HL_id <- wm_list[[1]][!is.na(wm_list[[1]])]
field_id_pts <- subset(HL_sf, HL_sf$Hylak_id %in% HL_id)
st_write(field_id_pts, "/home/ry4902/ISIMIP_Lake_Sector/HL_selected.shp")
