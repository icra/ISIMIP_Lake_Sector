"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : February 2022

Preprocess lake input data for ISIMIP3: convert .tif files into netCDF

"""

#%%
# import modules 
import xarray as xr
import numpy as np
from datetime import date
import gdal
import os 
from functions import *

# add today (for saving to netCDF later)
today = date.today()
date = today.strftime("%c")

# define directories 
parent_directory= os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/rmarce/ISIMIP_Lake_Sector/'
directory_netcdf =  parent_directory+'inputs_ISIMIP3_netcdf/' 
directory_raster =  parent_directory+'inputs_ISIMIP3/'


#%%
# ----------------------------------------------------------------
# Mean lakedepth
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Dmean_raster.tif'
filename_netcdf = directory_netcdf + 'Basic_inputs/' + 'mean_lakedepth.nc'

variable_name = 'mean_lakedepth'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Mean lake depth'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Mean depth for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Max lakedepth
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Dmax_raster.tif'
filename_netcdf = directory_netcdf + 'Basic_inputs/' + 'max_lakedepth.nc'

variable_name = 'max_lakedepth'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Maximum lake depth'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Maximum depth for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)


# %%
# ----------------------------------------------------------------
# Id hydrolakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Hylak_id.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'hydrolakes_id.nc'

variable_name = 'id'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The ID is shared with the GLOBathy database, so they are interoperable. A representative lake for each pixel was selected as explained in the github page of the url',
                        'title': 'ID of the representative lake at each pixel, from HydroLAKES and GLOBathy',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)


#%%
# ----------------------------------------------------------------
# Id hydrolakes big lakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'HLid_biglakes.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'hydrolakes_id_biglakes.nc'

variable_name = 'id_biglakes'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID for big lakes'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The ID is shared with the GLOBathy database, so they are interoperable. ',
                        'title': 'ID of very big lakes in HydroLAKES, for plotting purposes',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Area raster
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'A_raster.tif'
filename_netcdf = directory_netcdf + 'Basic_inputs/' + 'surface_area.nc'

variable_name = 'surface_area'
# variable attributes
attrs_variable = {'units': 'km^2', 'long_name' : 'Lake Surface Area'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Surface Area for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Volume raster
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'V_raster.tif'
filename_netcdf = directory_netcdf + 'Basic_inputs/' + 'volume.nc'

variable_name = 'volume'
# variable attributes
attrs_variable = {'units': 'km^3', 'long_name' : 'Lake Volume'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Volume for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)


# RAFA: better we do not include this, because it is not going to be used during simulations 
# #%%
# # ----------------------------------------------------------------
# # Volume development
# # ----------------------------------------------------------------

# # define filename, variable name and attributes
# filename_raster = directory_raster + 'Vd_raster.tif'
# filename_netcdf = directory_netcdf +'volume_development.nc'

# variable_name = 'volume_development'
# # variable attributes
# attrs_variable = {'units': '-', 'long_name' : 'volume development parameter (Vd)'}

# # global attributes
# attrs_global = {'creation_date': date,
#                         'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
#                         'title': 'Max lake and reservoir depth calculated from HydroLAKES',
#                         'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
#                         'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
#                         'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

# write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Lake type
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Lake_type.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'lake_type.nc'

variable_name = 'lake_type'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES Lake Type'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019',
                        'title': 'Lake Type from HydroLAKES. 1: Lake; 2: Reservoir; 3: Lake control (i.e. natural lake with regulation structure)',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

# %%
# ----------------------------------------------------------------
# Big lakes mask
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'biglakes_mask.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'biglakes_mask.nc'

variable_name = 'biglakes_mask'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'Biglakes mask'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019',
                        'title': 'Raster mask for very big lakes from HydroLAKES, for plotting purposes',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)



#%%
# ----------------------------------------------------------------
# area per level
# ----------------------------------------------------------------

# define filename, variable name and attributes
nlevels = 11

filename_rasters = [directory_raster + 'rasters_hypsographic/area_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_rasters_level = [directory_raster + 'rasters_hypsographic/level_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf + 'Hypsographic_curves/' + 'hypso_area.nc'

variable_name = 'area'
# variable attributes
attrs_variable = {'units': 'km^2', 'long_name' : 'Area per level'}
attrs_levels = {'units': 'm', 'long_name' : 'lake level from bottom'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Level-Area hypsographic information for each representative lake.',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global)

# %%
# ----------------------------------------------------------------
# volume per level
# ----------------------------------------------------------------

# define filename, variable name and attributes
nlevels = 11

filename_rasters = [directory_raster + 'rasters_hypsographic/volume_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_rasters_level = [directory_raster + 'rasters_hypsographic/level_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf + 'Hypsographic_curves/' + 'hypso_volume.nc'

variable_name = 'volume'
# variable attributes
attrs_variable = {'units': 'km^3', 'long_name' : 'Volume per level'}
attrs_levels = {'units': 'm', 'long_name' : 'lake level from bottom'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Level-Volume hypsographic information for each representative lake.',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global)

#%%
# ----------------------------------------------------------------
# power fit for hypsographic curve for volume 
# two parameters of the power fit, the third is the R2 of the fit.  
# ----------------------------------------------------------------
nlevels = 3

filename_rasters = [directory_raster + 'rasters_hypsographic/fitvolume_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf + 'Power_fit_to_hypsographic_curves/' + 'volume_fit.nc'

# variable attributes
attrs_parameter1 = {'units': '-', 'long_name' : 'parameter a of power fit for hypsographic of volume'}
attrs_parameter2 = {'units': '-', 'long_name' : 'parameter b of power fit for hypsographic of volume'}
attrs_R2 = {'units': '-', 'long_name' : 'R2 of power fit for hypsographic of volume'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Parameters (a and b) and goodness of fit (R2) for the fit V=ah^b, where h is lake lavel from the bottom',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

# delete if file exists
if os.path.isfile(filename_netcdf):
    os.system('rm '+filename_netcdf)


# isimip resolution hardcoded
resolution = 0.5

# read rasters 

raster = gdal.Open(filename_rasters[0])
values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
values_parameter1 = np.flipud(values_individual) 

raster = gdal.Open(filename_rasters[1])
values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
values_parameter2 = np.flipud(values_individual) 

raster = gdal.Open(filename_rasters[2])
values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
values_R2 = np.flipud(values_individual) 

# replace missing values (-3.4e+38 in raster) with NaNs
values_parameter2[values_parameter2<-1000] = np.nan
values_parameter1[values_parameter1<-1000] = np.nan
values_R2[values_R2<-1000] = np.nan

lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
lats= np.arange(-90+resolution/2,90+resolution/2,resolution)

lon_da = xr.DataArray(lons, 
                        coords = {'lon':lons}, 
                        dims='lon', 
                        attrs={'units':'degrees_east', 'axis':"X"})

lat_da = xr.DataArray(lats,
                        coords = {'lat':lats}, 
                        dims='lat', 
                        attrs={'units':'degrees_north', 'axis':"X"})

values_da_parameter1 = xr.DataArray(values_parameter1, 
                        coords = {'lon':lons,'lat':lats},
                        dims=('lat','lon'),
                        attrs = attrs_parameter1)

values_da_parameter2 = xr.DataArray(values_parameter2, 
                        coords = {'lon':lons,'lat':lats},
                        dims=('lat','lon'),
                        attrs = attrs_parameter2)

values_da_R2 = xr.DataArray(values_R2, 
                        coords = {'lon':lons,'lat':lats},
                        dims=('lat','lon'),
                        attrs = attrs_R2)

ds = xr.Dataset(data_vars={ 'lon' : lon_da,   
                            'lat' : lat_da,
                            'fit_parameter_a' : values_da_parameter1, 
                            'fit_parameter_b' : values_da_parameter2, 
                            'fit_R2' : values_da_R2, 
                            },
                            attrs=attrs_global)
                            
ds.to_netcdf(filename_netcdf, format='NETCDF4_CLASSIC',mode='w')




#%%
# ----------------------------------------------------------------
# power fit for hypsographic area
# two parameters of the power fit, the third is the R2 of the fit.  
# ----------------------------------------------------------------
nlevels = 3

filename_rasters = [directory_raster + 'rasters_hypsographic/fitarea_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf + 'Power_fit_to_hypsographic_curves/' + 'area_fit.nc'

# variable attributes
attrs_parameter1 = {'units': '-', 'long_name' : 'parameter a of power fit for hypsographic of volume'}
attrs_parameter2 = {'units': '-', 'long_name' : 'parameter b of power fit for hypsographic of volume'}
attrs_R2 = {'units': '-', 'long_name' : 'R2 of power fit for hypsographic of area'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1.',
                        'title': 'Parameters (a and b) and goodness of fit (R2) for the fit A=ah^b, where h is lake lavel from the bottom',
                       'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

# delete if file exists
if os.path.isfile(filename_netcdf):
    os.system('rm '+filename_netcdf)


# isimip resolution hardcoded
resolution = 0.5

# read rasters 

raster = gdal.Open(filename_rasters[0])
values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
values_parameter1 = np.flipud(values_individual) 

raster = gdal.Open(filename_rasters[1])
values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
values_parameter2 = np.flipud(values_individual) 

raster = gdal.Open(filename_rasters[2])
values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
values_R2 = np.flipud(values_individual) 

# replace missing values (-3.4e+38 in raster) with NaNs
values_parameter2[values_parameter2<-1000] = np.nan
values_parameter1[values_parameter1<-1000] = np.nan
values_R2[values_R2<-1000] = np.nan

lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
lats= np.arange(-90+resolution/2,90+resolution/2,resolution)

lon_da = xr.DataArray(lons, 
                        coords = {'lon':lons}, 
                        dims='lon', 
                        attrs={'units':'degrees_east', 'axis':"X"})

lat_da = xr.DataArray(lats,
                        coords = {'lat':lats}, 
                        dims='lat', 
                        attrs={'units':'degrees_north', 'axis':"X"})

values_da_parameter1 = xr.DataArray(values_parameter1, 
                        coords = {'lon':lons,'lat':lats},
                        dims=('lat','lon'),
                        attrs = attrs_parameter1)

values_da_parameter2 = xr.DataArray(values_parameter2, 
                        coords = {'lon':lons,'lat':lats},
                        dims=('lat','lon'),
                        attrs = attrs_parameter2)

values_da_R2 = xr.DataArray(values_R2, 
                        coords = {'lon':lons,'lat':lats},
                        dims=('lat','lon'),
                        attrs = attrs_R2)

ds = xr.Dataset(data_vars={ 'lon' : lon_da,   
                            'lat' : lat_da,
                            'fit_parameter_a' : values_da_parameter1, 
                            'fit_parameter_b' : values_da_parameter2, 
                            'fit_R2' : values_da_R2, 
                            },
                            attrs=attrs_global)
                            
ds.to_netcdf(filename_netcdf, format='NETCDF4_CLASSIC',mode='w')
# %%
