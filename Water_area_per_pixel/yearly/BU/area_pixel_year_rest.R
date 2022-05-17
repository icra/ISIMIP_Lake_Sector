Sys.time()
library(sf); library(raster); library(lwgeom); library(rgdal)

#open grand database
grand <- st_read("/home/dmercado/ISIMIP_Lake_Sector_inputs/GRanD_dams_v1_1.shp")
#select id with year of construction greater or equal than 1850 (historical scenario first year)
grand_id <- grand$GRAND_ID[grand$YEAR>=1850]
#open hydrolakes database
HL <- st_read("/home/dmercado/ISIMIP_Lake_Sector_inputs/HydroLAKES_polys_v10.shp")
#subset the previous selected id in grand
HL_dam <- subset(HL, HL$Grand_id %in% grand_id)
#no concuerdan en numero, hay menos en HL de los seleccionados
length(HL_dam[[1]])==length(grand_id)
#years with dam construction beoynd 1850
years_dam <- sort(unique(grand$YEAR[grand$GRAND_ID %in% grand_id]))

for (y in years_dam){  
  grand_id_y <- grand$GRAND_ID[grand$YEAR==y]
  HL_dam_y <- subset(HL_dam, HL_dam$Grand_id %in% grand_id_y)
  
  raster_matrix <- matrix(data = 1, nrow=360, ncol = 720)
  rasterHL <- raster(raster_matrix)
  extent(rasterHL) <- extent(c(-180,180,-90,90))
  crs(rasterHL) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  
  #convert grid raster to polygons
  ras_pol <- rasterToPolygons(rasterHL)
  #convert to sf format
  ras_pol_sf <- st_as_sf(ras_pol) 
  
  #some polygons need repairing before intersection
  HL_dam_y_valid <- st_make_valid(HL_dam_y)
  st_agr(HL_dam_y_valid) = "constant"
  st_agr(ras_pol_sf) = "constant"
  
  #intersection of the two. Lake polygons get cutter at pixel borders
  union_st <- st_intersection(ras_pol_sf, HL_dam_y_valid)
  #we store the area of the intersected lake polygons as a new field
  union_st$area_cut <- st_area(union_st)
  
  #a handful (64 out of 1,515,528) of features are points or strings with no area, associated to more that one grid cell. Must be removed.
  #delete_positions <- which(as.numeric(union_st$area_cut)==0)
  #union_st <- union_st[-delete_positions,]
  
  #this tells you which grid cell covers each intersected polygon
  aaa <- st_covered_by(union_st, ras_pol_sf)
  #area of each grid cell
  ras_pol_sf$area_cell <- st_area(ras_pol_sf)
  #initializing the collector for the area of all lakes ina pixel
  ras_pol_sf$areas_lakes <- NA
  
  #loop for assigning intersectoed polygons areas to pixels
  for (i in 1:dim(ras_pol_sf)[1]){
    locations<-which(as.numeric(aaa)==i)
    ras_pol_sf$areas_lakes[i] <- sum(union_st$area_cut[locations])
  }
  #calculation of the fraction occupied by lakes
  ras_pol_sf$frac_water <- ras_pol_sf$areas_lakes / as.numeric(ras_pol_sf$area_cell)
  
  #preparing and saving the output raster with fraction water coverage
  raster_salida <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs", resolution=c(0.5,0.5), vals=NULL)
  raster_areas <- rasterize(ras_pol_sf, raster_salida,'frac_water')
  writeRaster(raster_areas,paste0("frac_areas_",y,".tif"), overwrite=T)
  print(Sys.time())
  
  ####################added a posteriori to remove zeros and small deviations >1 
  raster_frac <- raster(paste0("frac_areas_",y,".tif"))
  m <- c(0, 0, NA,1,2,1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  raster_frac_NA <- reclassify(raster_frac, rclmat, include.lowest=TRUE)
  writeRaster(raster_frac_NA,paste0("frac_areas_NA_",y,".tif"), overwrite=T)
  #################################################  
}



