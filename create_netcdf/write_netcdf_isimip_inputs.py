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
from osgeo import gdal
import os 
from functions import *

# add today (for saving to netCDF later)
today = date.today()
date = today.strftime("%c")

# define directories 
parent_directory= os.path.abspath(os.path.join(os.getcwd(), os.pardir)) 
directory_netcdf =  parent_directory+'/create_netcdf/netcdfs/' 
directory_raster =  parent_directory+'/inputs_ISIMIP3/'


#%%
# ----------------------------------------------------------------
# Mean lakedepth
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Dmean_raster.tif'
filename_netcdf = directory_netcdf +'mean_lakedepth.nc'

variable_name = 'mean_lakedepth'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Mean lake depth'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Mean lake and reservoir depth calculated from HydroLAKES',
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
filename_netcdf = directory_netcdf +'max_lakedepth.nc'

variable_name = 'max_lakedepth'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Maximum lake depth'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)


# %%
# ----------------------------------------------------------------
# Id hydrolakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Hylak_id.tif'
filename_netcdf = directory_netcdf +'hydrolakes_id.nc'

variable_name = 'id'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)


#%%
# ----------------------------------------------------------------
# Id hydrolakes big lakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'HLid_biglakes.tif'
filename_netcdf = directory_netcdf +'hydrolakes_id_biglakes.nc'

variable_name = 'id_biglakes'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID for big lakes'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Id hydrolakes big lakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'HLid_biglakes.tif'
filename_netcdf = directory_netcdf +'hydrolakes_id_biglakes.nc'

variable_name = 'id_biglakes'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID for big lakes'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
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
filename_netcdf = directory_netcdf +'surface_area.nc'

variable_name = 'surface_area'
# variable attributes
attrs_variable = {'units': 'm^2', 'long_name' : 'HydroLAKES ID for big lakes'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Area raster
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'V_raster.tif'
filename_netcdf = directory_netcdf +'lake_total_volume.nc'

variable_name = 'volume'
# variable attributes
attrs_variable = {'units': 'm^3', 'long_name' : 'HydroLAKES ID for big lakes'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Volume development
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Vd_raster.tif'
filename_netcdf = directory_netcdf +'volume_development.nc'

variable_name = 'volume_development'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'volume development parameter (Vd)'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Lake type
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Lake_type.tif'
filename_netcdf = directory_netcdf +'lake_type.nc'

variable_name = 'lake_type'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'Lake Type'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

# %%
# ----------------------------------------------------------------
# Big lakes mask
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'biglakes_mask.tif'
filename_netcdf = directory_netcdf +'biglakes_mask.nc'

variable_name = 'biglakes_mask'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'Biglakes mask'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
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

filename_rasters = [directory_raster + 'rasters/area_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_rasters_level = [directory_raster + 'rasters/level_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf +'lake_area.nc'

variable_name = 'area'
# variable attributes
attrs_variable = {'units': 'm^2', 'long_name' : 'area'}
attrs_levels = {'units': 'm', 'long_name' : 'lake depth level'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
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

filename_rasters = [directory_raster + 'rasters/volume_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_rasters_level = [directory_raster + 'rasters/level_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf +'lake_volume.nc'

variable_name = 'volume'
# variable attributes
attrs_variable = {'units': 'm^3', 'long_name' : 'lake volume per level '}
attrs_levels = {'units': 'm', 'long_name' : 'lake depth level'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global)

#%%
# ----------------------------------------------------------------
# fit volume 
# ----------------------------------------------------------------

# define filename, variable name and attributes
nlevels = 3

filename_rasters = [directory_raster + 'rasters/fitvolume_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_rasters_level = [directory_raster + 'rasters/level_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
filename_netcdf = directory_netcdf +'lake_volume.nc'

variable_name = 'volume'
# variable attributes
attrs_variable = {'units': 'm^3', 'long_name' : 'lake volume per level '}
attrs_levels = {'units': 'm', 'long_name' : 'lake depth level'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019, The lakedepth includes depths from reservoirs (included in HydroLAKES) and is rasterised for shallow lakes, while big lakes are assigned their unique value to all grid cells they cover.',
                        'title': 'Max lake and reservoir depth calculated from HydroLAKES',
                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global)


# %%
# ----------------------------------------------------------------
# fit area
# ----------------------------------------------------------------

