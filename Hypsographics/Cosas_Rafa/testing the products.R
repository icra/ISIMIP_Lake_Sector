library(ncdf4)
library(raster)
library(viridis)

nc_data <- nc_open("/home/rmarce/ISIMIP_Lake_Sector/inputs_ISIMIP3_netcdf/Reference/biglakes_mask.nc")

{
  sink('gimms3g_ndvi_1982-2012_metadata.txt')
  print(nc_data)
  sink()
}

attri_nc <- ncvar_get(nc_data, "biglakes_mask")


tmp_raster <- brick(nc_data, varname="biglakes_mask")


raster_bl <- raster(attri_nc)
extent(raster_bl) <- extent(c(-180,180,-90,90))
crs(raster_bl) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(raster_bl,col=viridis(2))
