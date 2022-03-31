
"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : February 2022

functions to convert .tif files into netCDF

"""

# import modules 
import xarray as xr
import numpy as np
from datetime import date
import gdal
import os 

# ----------------------------------------------------------------
# Function to write 2d netCDF files

def write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global):

    # delete if file exists
    #if os.path.isfile(filename_netcdf):
       # os.system('rm '+filename_netcdf)


    # isimip resolution hardcoded
    resolution = 0.5

    # read raster 
    raster = gdal.Open(filename_raster)
    values = np.array(raster.GetRasterBand(1).ReadAsArray())
    values = np.flipud(values)

# replace missing values (-3.4e+38 in raster) with NaNs
    #values[values<-1000] = np.nan
    values[values<-1000] = 1.e+20


    lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
    lats= np.arange(-90+resolution/2,90+resolution/2,resolution)

    lon_da = xr.DataArray(lons, 
                            coords = {'lon':lons}, 
                            dims='lon', 
                            attrs={'units':'degrees_east', 'axis':"X"})

    lat_da = xr.DataArray(lats,
                            coords = {'lat':lats}, 
                            dims='lat', 
                            attrs={'units':'degrees_north', 'axis':"Y"})


    values_da = xr.DataArray(values, 
                            coords = {'lon':lons,'lat':lats},
                            dims=('lat','lon'),
                            attrs = attrs_variable)


    ds = xr.Dataset(data_vars={ 'lon' : lon_da,   
                                'lat' : lat_da,
                                variable_name : values_da},
                                attrs=attrs_global)
    
    
    
    ds.to_netcdf(filename_netcdf, format='NETCDF4_CLASSIC',mode='w')


# ----------------------------------------------------------------
# Function to write 3d netCDF files

def write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global):

    # delete if file exists
    if os.path.isfile(filename_netcdf):
        os.system('rm '+filename_netcdf)


    # isimip resolution hardcoded
    resolution = 0.5

    # read rasters - if it is a list, read different rasters and concatenate them
    values_dict = {}

    for n,filename_raster in enumerate(filename_rasters):
        raster = gdal.Open(filename_raster)
        values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
        values_dict[n] = np.flipud(values_individual) 
    values = np.stack(values_dict.values())


    # read raster level - if it is a list, read different rasters and concatenate them

    level_dict = {}
    for n,filename_raster_individual in enumerate(filename_rasters_level):
        raster = gdal.Open(filename_raster_individual)
        values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
        level_dict[n] = np.flipud(values_individual) 
    levels = np.stack(level_dict.values())


# replace missing values (-3.4e+38 in raster) with NaNs
    values[values<-1000] = 1.e+20
    levels[levels<-1000] = 1.e+20
 
    lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
    lats= np.arange(-90+resolution/2,90+resolution/2,resolution)
    
    levlaks = np.arange(1,np.shape(values)[0]+1)

    lon_da = xr.DataArray(lons, 
                            coords = {'lon':lons}, 
                            dims='lon', 
                            attrs={'units':'degrees_east', 'axis':"X"})

    lat_da = xr.DataArray(lats,
                            coords = {'lat':lats}, 
                            dims='lat', 
                            attrs={'units':'degrees_north', 'axis':"Y"})


    values_da = xr.DataArray(values, 
                            coords = {'levlak':levlaks,'lon':lons,'lat':lats},
                            dims=('levlak','lat','lon'),
                            attrs = attrs_variable)

    levels_da = xr.DataArray(levels,
                            coords = {'levlak':levlaks,'lon':lons,'lat':lats},
                            dims=('levlak','lat','lon'),
                            attrs = attrs_levels)

    ds = xr.Dataset(data_vars={ 'lon' : lon_da,   
                                'lat' : lat_da,
                                variable_name : values_da, 
                                'levels': levels_da},
                                attrs=attrs_global)
                                
    ds.to_netcdf(filename_netcdf, format='NETCDF4_CLASSIC',mode='w')
