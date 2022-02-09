library(sf);library(raster)

#weighted median function, modified from https://github.com/spatstat/spatstat.geom/blob/main/R/weightedStats.R
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
  #This method is apply when having more than 2 points
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

#Opening HydroLAKES
HL <- shapefile("/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/test_weightedMean/HL_test.shp")

#Calculate "pseudocentroid" inside polygon  
HL_sf <- st_as_sf(HL) 
HL_cent <- st_point_on_surface(HL_sf) 

#get attribute table and add coordinates:
HL_df <- data.frame(st_drop_geometry(HL_cent), st_coordinates(HL_cent))

#Filling new matrix for each longitude and latitude with the weighthed mean values
wm_matrix <- matrix(data = NA, nrow=360, ncol = 720) #matrix to save all weighted mean map
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
        print(paste("lagos mayores de 1degree estaban en estas posiciones:", pos_lakes))
        if (length(pos_lakes)>1){
          print("hay dos lagos grandes OJO")
          stop("hay dos lagos grandes OJO")
        }
        wm_matrix[c_lat,c_lon] <- HL_df$Depth_avg[pos_lakes]
      }else{
        #print(pos_lakes)
        if (length(pos_lakes)==1){
          wm_temp <- HL_df$Depth_avg[pos_lakes]
          wm_matrix[c_lat,c_lon] <- wm_temp
          print(paste("1111111111111111ésta es la wm:", wm_temp))
          print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
          print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
        }else if(length(pos_lakes)==2){
          pos2_max <- which(HL_df$Lake_area[pos_lakes]==max(HL_df$Lake_area[pos_lakes]))
          pos_lakes <- pos_lakes[pos2_max]
          wm_temp <- HL_df$Depth_avg[pos_lakes]
          wm_matrix[c_lat,c_lon] <- wm_temp
          print(paste("222222222222ésta es la wm:", wm_temp))
          print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
          print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
          
        }else{
          wm_temp <- weighted_median(x=HL_df$Depth_avg[pos_lakes], w=HL_df$Lake_area[pos_lakes])
          wm_matrix[c_lat,c_lon] <- wm_temp[1]
          print(paste("ésta es la wm:", wm_temp[1], "con esta posición", wm_temp[2]))
          print(paste("a partir de estas depths:", paste(HL_df$Depth_avg[pos_lakes], collapse = ' ')))
          print(paste("y estas areas:", paste(HL_df$Lake_area[pos_lakes], collapse = ' ')))
        }
      }
    }
  }
}

#write.table(wm_matrix, file = "/home/ry4902/Documents/ISIMIP/AnalysisDepthArea/test_weightedMean/wm_matrix.ascii", col.names = F, row.names = F)
wm_raster <- raster(wm_matrix)
extent(wm_raster) <- extent(c(-180,180,-90,90))
crs(wm_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

writeRaster(wm_raster,"/home/ry4902/ISIMIP_Lake_Sector/test.tif", overwrite=T)

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

#weighted median function, modified from (source:https://stackoverflow.com/questions/2748725/is-there-a-weighted-median-function)
weighted_median <- function(x, w) {
  w_o <- w[order(x)]
  x_o <- x[order(x)]
  
  prob <- cumsum(w_o)/sum(w_o)
  ps <- which(abs(prob - .5) == max(abs(prob - .5)))
  
  #get position
  pos <- which(x==x_o[ps] & w==w_o[ps])
  
  return(c(x_o[ps], pos))
}



#weighted median function, modified from https://github.com/spatstat/spatstat.geom/blob/main/R/weightedStats.R
weighted_median <- function(x, w, na.rm=TRUE) {
  probs <- 0.5
  
  #if type=1, the next largest value is returned (this is the right-continuous inverse of the left-continuous cumulative distribution function);
  
  #if type=2, the average of the two surrounding values is returned (the average of the right-continuous and left-continuous inverses);
  
  #if type=4, linear interpolation is performed
  
  #type <- 4
  
  #x <- as.numeric(as.vector(x))
  #w <- as.numeric(as.vector(w))
  #if(is.na(m <- match(type, c(1,2,4))))
  #  stop("Argument 'type' must equal 1, 2 or 4", call.=FALSE)
  #type <- c(1,2,4)[m]
  #if(anyNA(x) || anyNA(w)) {
  #  ok <- !(is.na(x) | is.na(w))
  #  x <- x[ok]
  #  w <- w[ok]
  #  }
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x_o <- x[oo]
  w_o <- w[oo]
  Fx <- cumsum(w_o)/sum(w_o)
  #' 
  if(anyDuplicated(x)) {
    dup <- rev(duplicated(rev(x_o)))
    x_o <- x_o[!dup]
    Fx <- Fx[!dup]
  }
  #'
  out <- approx(Fx, x_o, xout=probs, ties="ordered", rule=2,
                method="linear")
  #out <- switch(as.character(type),
  #              "1" = approx(Fx, x, xout=probs, ties="ordered", rule=2,
  #                           method="constant", f=1),
  #              "2" = approx(Fx, x, xout=probs, ties="ordered", rule=2,
  #                           method="constant", f=1/2),
  #              "4" = approx(Fx, x, xout=probs, ties="ordered", rule=2,
  #                           method="linear"))
  
  closest_pos <- which(abs(out$y-x)==min(abs(out$y-x))) #posición del valor más cercano
  random_closest_pos <- sample(closest_pos, 1) #randomly selection in case there are more than one depths with the same value
  result <- c(x[random_closest_pos], random_closest_pos) #saving weighted median depth and its position in initial 'x' vector
  #names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}


out <- switch(as.character(type),
              "1" = approx(Fx, x, ties="ordered", rule=2,
                           method="constant", f=1),
              "2" = approx(Fx, x, ties="ordered", rule=2,
                           method="constant", f=1/2),
              "4" = approx(Fx, x, ties="ordered", rule=2,
                           method="linear"))

