library(sf);library(raster)

#Opening HydroLAKES
HL <- shapefile("/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/test_weightedMean/HL_test.shp")

#Calculate "pseudocentroid" inside polygon  
HL_sf <- st_as_sf(HL) 
HL_cent <- st_point_on_surface(HL_sf) 

#get attribute table and add coordinates:
HL_df <- data.frame(st_drop_geometry(HL_cent), st_coordinates(HL_cent))

#weighted median function, modified from (source:https://stackoverflow.com/questions/2748725/is-there-a-weighted-median-function)
weighted_median <- function(x, w) {
  w_o <- w[order(x)]
  x_o <- x[order(x)]
  
  prob <- cumsum(w_o)/sum(w_o)
  ps <- which(abs(prob - .5) == min(abs(prob - .5)))
  
  #get position
  pos <- which(x==x_o[ps] & w==w_o[ps])
  
  return(c(x_o[ps], pos))
}
##3QUITAR BIG LAKES
#AGREGAR IGUAL
wm_matrix <- matrix(data = NA, nrow=360, ncol = 720)
c_lon <- 0; 
for (lon in seq(-180, 180, 0.5)){ #loop in longitude
  c_lon <- c_lon+1; c_lat <- 0
  for (lat in seq(90, -90, -0.5)){ #loop in latitude
    c_lat <- c_lat+1
    pos_lakes <- which( (HL_df$X>lon & HL_df$X<(lon+0.5)) & (HL_df$Y>lat&HL_df$Y<(lat+0.5))) 
    #columns to save 1:21
    if (length(pos_lakes)>0){
      #print(pos_lakes)
      wm_temp <- weighted_median(x=HL_df$Depth_avg[pos_lakes], w=HL_df$Lake_area[pos_lakes])[1]
      wm_matrix[c_lat-1,c_lon] <- wm_temp
      print(paste("ésta es la wm:", wm_temp[1], "con esta posición", wm_temp[2]))
      print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
      print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
    }
  }
}

#write.table(wm_matrix, file = "/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/test_weightedMean/wm_matrix.ascii", col.names = F, row.names = F)
wm_raster <- raster(wm_matrix)
extent(wm_raster) <- extent(c(-180,180,-90,90))
#extent(wm_raster) <- extent(c(-180,180,-90,90))
res(wm_raster)
crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

writeRaster(wm_raster,"/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/test_weightedMean/test.tif", overwrite=T)

#Borrar
###basura que puede ser util:
writeOGR(out, "./home/ry4902/Documents/ISIMIP/AnalysisDepthArea", "test", 
         driver = "ESRI Shapefile")

t<-(data.frame(unlist(matrix(out$geometry))))
df1 <- data.frame()
for (i in 1:61){
  print(unlist(t[i,]))
  df1 <- rbind(df1, unlist(t[i,]))
}

write.csv(df1, file = "/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/centroids_Big3080km2_inside.csv",
          row.names = F)

length(HL_cent$geometry) #total number of lakes
unlist(HL_cent$geometry[l]) #Longitude and latitude of l lake

write.csv(df1, file = "/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/centroids_Big3080km2_inside.csv",
          row.names = F)


