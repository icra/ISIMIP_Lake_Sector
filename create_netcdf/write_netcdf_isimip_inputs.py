"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be), modified by Rafael Marcé (rmarce@icra.cat) and Daniel Mercado-Bettín (dmercado@icra.cat)
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
directory_netcdf =  parent_directory+'/inputs_ISIMIP3_netcdf/' 
directory_raster =  parent_directory+'/inputs_ISIMIP3/'
directory_waterarea_per_pixel =  parent_directory+'/Water_area_per_pixel/yearly/final/'


#%%
# ----------------------------------------------------------------
# Mean lakedepth
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Dmean_raster.tif'
filename_netcdf = directory_netcdf + 'Alternative_inputs/' + 'mean_lakedepth.nc'

variable_name = 'mean_lakedepth'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Mean lake depth', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
                        'title': 'Mean depth for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'use': 'We reccommend using the area o volume hypsographic curves provided in this repository as inputs for your lake model. Use this file only if your lake model does not accept a full hypsographic curve as an input.',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)


#%%
# ----------------------------------------------------------------
# Max lakedepth
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Dmax_raster.tif'
filename_netcdf = directory_netcdf + 'Alternative_inputs/' + 'max_lakedepth.nc'

variable_name = 'max_lakedepth'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Maximum lake depth', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
                        'title': 'Maximum depth for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'use': 'We reccommend using the area o volume hypsographic curves provided in this repository as inputs for your lake model. Use this file only if your lake model does not accept a full hypsographic curve as an input.',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)

# %%
# ----------------------------------------------------------------
# Id hydrolakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'Hylak_id.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'hydrolakes_id.nc'

variable_name = 'id'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019',
                        'title': 'ID of the representative lake at each pixel, from HydroLAKES and GLOBathy',
                        'use': 'File provided for reference, to relate HydroLAKES and GLOBathy database fields to the ISIMIP3 representative lakes',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)

#%%
# ----------------------------------------------------------------
# Id hydrolakes big lakes
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'HLid_biglakes.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'hydrolakes_id_biglakes.nc'

variable_name = 'id_biglakes'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES ID for big lakes', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019',
                        'title': 'ID of very big lakes in HydroLAKES, for plotting purposes',
                        'use': 'File provided for reference, to produce global plots with conspicuous large lakes during figure production. To be used together with the file storing the big lakes mask',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)
#%%
# ----------------------------------------------------------------
# Area raster
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'A_raster.tif'
filename_netcdf = directory_netcdf + 'Alternative_inputs/' + 'surface_area.nc'

variable_name = 'surface_area'
# variable attributes
attrs_variable = {'units': 'km^2', 'long_name' : 'Lake Surface Area'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
                        'title': 'Surface Area for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'use': 'We reccommend using the area o volume hypsographic curves provided in this repository as inputs for your lake model. Use this file only if your lake model does not accept a full hypsographic curve as an input.',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)
#%%
# ----------------------------------------------------------------
# Volume raster
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'V_raster.tif'
filename_netcdf = directory_netcdf + 'Alternative_inputs/' + 'volume.nc'

variable_name = 'volume'
# variable attributes
attrs_variable = {'units': 'km^3', 'long_name' : 'Lake Volume', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
                        'title': 'Volume for ISIMIP3 representative lakes calculated from GLOBathy and HydroLAKES',
                        'use': 'We reccommend using the area o volume hypsographic curves provided in this repository as inputs for your lake model. Use this file only if your lake model does not accept a full hypsographic curve as an input.',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20  ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)

#%%
# ----------------------------------------------------------------
# Heigth or m.a.s.l (added by Daniel)
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'H_raster.tif'
filename_netcdf = directory_netcdf + 'Alternative_inputs/' + 'elevation.nc'

variable_name = 'elevation'
# variable attributes
attrs_variable = {'units': 'm', 'long_name' : 'Lake Heigth', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019', 
                        'title': 'Heigth (meters above sea level) for ISIMIP3 representative lakes calculated from HydroLAKES',
                        'use': 'We recommend using elevation provided in this repository as inputs for your lake model. Use this file only if your lake model needs elevation values.',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)




# RAFA: better we do not include this, because it is not going to be used during simulations 
# #%%
# # ----------------------------------------------------------------
# # Volume development
## RAFA: we finally decided not to include this datset in the final set to avoid confusions 
# # ----------------------------------------------------------------

# # define filename, variable name and attributes
# filename_raster = directory_raster + 'Vd_raster.tif'
# filename_netcdf = directory_netcdf +'volume_development.nc'

# variable_name = 'volume_development'
# # variable attributes
# attrs_variable = {'units': '-', 'long_name' : 'volume development parameter (Vd)'}

# # global attributes
# attrs_global = {'creation_date': date,
#                         'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
#                         'title': 'Max lake and reservoir depth calculated from HydroLAKES',
#                         'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be)',
#                         'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
#                         'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

# write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)

#%%
# ----------------------------------------------------------------
# Lake type
## RAFA: we finally decided not to include this dataset in the final set to avoid confusions 
# ----------------------------------------------------------------
#
## define filename, variable name and attributes
#filename_raster = directory_raster + 'Lake_type.tif'
#filename_netcdf = directory_netcdf + 'Reference/' + 'lake_type.nc'
#
#variable_name = 'lake_type'
## variable attributes
#attrs_variable = {'units': '-', 'long_name' : 'HydroLAKES Lake Type'}
#
## global attributes
#attrs_global = {'creation_date': date,
#                        'source': 'HydroLAKES polygons dataset v1.0 June 2019',
#                        'title': 'Lake Type from HydroLAKES. 1: Lake; 2: Reservoir; 3: Lake control (i.e. natural lake with regulation structure)',
#                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
#                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
#                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector'}
#
#write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
#os.system('cdo -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
#os.system('rm '+ filename_netcdf)
#os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)
# %%
# ----------------------------------------------------------------
# Big lakes mask
# ----------------------------------------------------------------

# define filename, variable name and attributes
filename_raster = directory_raster + 'biglakes_mask.tif'
filename_netcdf = directory_netcdf + 'Reference/' + 'biglakes_mask.nc'

variable_name = 'biglakes_mask'
# variable attributes
attrs_variable = {'units': '-', 'long_name' : 'Biglakes mask', '_FillValue' : 1e20}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES polygons dataset v1.0 June 2019',
                        'title': 'Raster mask for very big lakes from HydroLAKES, for plotting purposes',
                        'use': 'File provided for reference, to produce global plots with conspicuous large lakes during figure production. To be used together with the file storing the big lakes IDs',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_2d(filename_raster,filename_netcdf,attrs_variable,variable_name,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)


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
attrs_variable = {'units': 'km^2', 'long_name' : 'Area per level', '_FillValue' : 1e20}
attrs_levels = {'units': 'm', 'long_name' : 'lake level from bottom'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
                        'title': 'Level-Area hypsographic information for each representative lake.',
                        'use': 'Level-Area Hypsographic curves to be provided to your model (alternatively, you can also use the Level-Volume Hypsographic curves in this repository). Each hypsographic consist in 11 data pairs. Level refers to the level of the lake taking the lake bottom as the reference (in meters), Area is the area at the corresponding level (in km2).',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)
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
attrs_variable = {'units': 'km^3', 'long_name' : 'Volume per level', '_FillValue' : 1e20}
attrs_levels = {'units': 'm', 'long_name' : 'lake level from bottom'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
                        'title': 'Level-Volume hypsographic information for each representative lake.',
                        'use': 'Level-Volume Hypsographic curves to be provided to your model (alternatively, you can also use the Level-Area Hypsographic curve in this repository). Each hypsographic consist in 11 data pairs. Level refers to the level of the lake taking the lake bottom as the reference (in meters), Volume is the volume at the corresponding level (in km3).',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }

write_netcdf_3d(filename_rasters,filename_rasters_level,filename_netcdf,attrs_variable,variable_name,attrs_levels,attrs_global)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)


##%%
## ----------------------------------------------------------------
## power fit for hypsographic curve for volume 
## two parameters of the power fit, the third is the R2 of the fit. 
## RAFA: we finally decided not to include this option in the final set to avoid confusions 
## ----------------------------------------------------------------
#nlevels = 3
#
#filename_rasters = [directory_raster + 'rasters_hypsographic/fitvolume_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
#filename_netcdf = directory_netcdf + 'Power_fit_to_hypsographic_curves/' + 'volume_fit.nc'
#
## variable attributes
#attrs_parameter1 = {'units': '-', 'long_name' : 'parameter a of power fit for hypsographic of volume'}
#attrs_parameter2 = {'units': '-', 'long_name' : 'parameter b of power fit for hypsographic of volume'}
#attrs_R2 = {'units': '-', 'long_name' : 'R2 of power fit for hypsographic of volume'}
#
## global attributes
#attrs_global = {'creation_date': date,
#                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
#                        'title': 'Parameters (a and b) and goodness of fit (R2) for the fit V=ah^b, where h is lake lavel from the bottom',
#                        'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
#                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
#                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }
#
## delete if file exists
#if os.path.isfile(filename_netcdf):
#    os.system('rm '+filename_netcdf)
#
#
## isimip resolution hardcoded
#resolution = 0.5
#
## read rasters 
#
#raster = gdal.Open(filename_rasters[0])
#values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
#values_parameter1 = np.flipud(values_individual) 
#
#raster = gdal.Open(filename_rasters[1])
#values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
#values_parameter2 = np.flipud(values_individual) 
#
#raster = gdal.Open(filename_rasters[2])
#values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
#values_R2 = np.flipud(values_individual) 
#
## replace missing values (-3.4e+38 in raster) with NaNs
#values_parameter2[values_parameter2<-1000] = np.nan
#values_parameter1[values_parameter1<-1000] = np.nan
#values_R2[values_R2<-1000] = np.nan
#
#lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
#lats= np.arange(-90+resolution/2,90+resolution/2,resolution)
#
#lon_da = xr.DataArray(lons, 
#                        coords = {'lon':lons}, 
#                        dims='lon', 
#                        attrs={'units':'degrees_east', 'axis':"X"})
#
#lat_da = xr.DataArray(lats,
#                        coords = {'lat':lats}, 
#                        dims='lat', 
#                        attrs={'units':'degrees_north', 'axis':"X"})
#
#values_da_parameter1 = xr.DataArray(values_parameter1, 
#                        coords = {'lon':lons,'lat':lats},
#                        dims=('lat','lon'),
#                        attrs = attrs_parameter1)
#
#values_da_parameter2 = xr.DataArray(values_parameter2, 
#                        coords = {'lon':lons,'lat':lats},
#                        dims=('lat','lon'),
#                        attrs = attrs_parameter2)
#
#values_da_R2 = xr.DataArray(values_R2, 
#                        coords = {'lon':lons,'lat':lats},
#                        dims=('lat','lon'),
#                        attrs = attrs_R2)
#
#ds = xr.Dataset(data_vars={ 'lon' : lon_da,   
#                            'lat' : lat_da,
#                            'fit_parameter_a' : values_da_parameter1, 
#                            'fit_parameter_b' : values_da_parameter2, 
#                            'fit_R2' : values_da_R2, 
#                            },
#                            attrs=attrs_global)
#                            
#ds.to_netcdf(filename_netcdf, format='NETCDF4_CLASSIC',mode='w')
#os.system('cdo -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
#os.system('rm '+ filename_netcdf)
#os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)
#
#
#
##%%
## ----------------------------------------------------------------
## power fit for hypsographic area
## two parameters of the power fit, the third is the R2 of the fit.  
## RAFA: we finally decided not to include this option in the final set to avoid confusions 
## ----------------------------------------------------------------
#nlevels = 3
#
#filename_rasters = [directory_raster + 'rasters_hypsographic/fitarea_level'+str(n)+'.tif' for n in range(1,nlevels+1)]
#filename_netcdf = directory_netcdf + 'Power_fit_to_hypsographic_curves/' + 'area_fit.nc'
#
## variable attributes
#attrs_parameter1 = {'units': '-', 'long_name' : 'parameter a of power fit for hypsographic of volume'}
#attrs_parameter2 = {'units': '-', 'long_name' : 'parameter b of power fit for hypsographic of volume'}
#attrs_R2 = {'units': '-', 'long_name' : 'R2 of power fit for hypsographic of area'}
#
## global attributes
#attrs_global = {'creation_date': date,
#                        'source': 'GLOBathy, the Global Lakes Bathymetry Dataset, January 2022, https://doi.org/10.6084/m9.figshare.c.5243309.v1., and HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603', 
#                        'title': 'Parameters (a and b) and goodness of fit (R2) for the fit A=ah^b, where h is lake lavel from the bottom',
#                       'contact' : 'Daniel Mercado - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
#                        'references':'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022). GLOBathy, the global lakes bathymetry dataset. Scientific Data, 9(1), 36. https://doi.org/10.1038/s41597-022-01132-9',
#                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }
#
## delete if file exists
#if os.path.isfile(filename_netcdf):
#    os.system('rm '+filename_netcdf)
#
#
## isimip resolution hardcoded
#resolution = 0.5
#
## read rasters 
#
#raster = gdal.Open(filename_rasters[0])
#values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
#values_parameter1 = np.flipud(values_individual) 
#
#raster = gdal.Open(filename_rasters[1])
#values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
#values_parameter2 = np.flipud(values_individual) 
#
#raster = gdal.Open(filename_rasters[2])
#values_individual = np.array(raster.GetRasterBand(1).ReadAsArray())
#values_R2 = np.flipud(values_individual) 
#
## replace missing values (-3.4e+38 in raster) with NaNs
#values_parameter2[values_parameter2<-1000] = np.nan
#values_parameter1[values_parameter1<-1000] = np.nan
#values_R2[values_R2<-1000] = np.nan
#
#lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
#lats= np.arange(-90+resolution/2,90+resolution/2,resolution)
#
#lon_da = xr.DataArray(lons, 
#                        coords = {'lon':lons}, 
#                        dims='lon', 
#                        attrs={'units':'degrees_east', 'axis':"X"})
#
#lat_da = xr.DataArray(lats,
#                        coords = {'lat':lats}, 
#                        dims='lat', 
#                        attrs={'units':'degrees_north', 'axis':"X"})
#
#values_da_parameter1 = xr.DataArray(values_parameter1, 
#                        coords = {'lon':lons,'lat':lats},
#                        dims=('lat','lon'),
#                        attrs = attrs_parameter1)
#
#values_da_parameter2 = xr.DataArray(values_parameter2, 
#                        coords = {'lon':lons,'lat':lats},
#                        dims=('lat','lon'),
#                        attrs = attrs_parameter2)
#
#values_da_R2 = xr.DataArray(values_R2, 
#                        coords = {'lon':lons,'lat':lats},
#                        dims=('lat','lon'),
#                        attrs = attrs_R2)
#
#ds = xr.Dataset(data_vars={ 'lon' : lon_da,   
#                            'lat' : lat_da,
#                            'fit_parameter_a' : values_da_parameter1, 
#                            'fit_parameter_b' : values_da_parameter2, 
#                            'fit_R2' : values_da_R2, 
#                            },
#                            attrs=attrs_global)
#                            
#ds.to_netcdf(filename_netcdf, format='NETCDF4_CLASSIC',mode='w')
#os.system('cdo -O --history -setmissval,1e+20 ' + filename_netcdf + ' ' + filename_netcdf + '4')
#os.system('rm '+ filename_netcdf)
#os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)
## %%

# %%
# ----------------------------------------------------------------
# transient water area per year 
# ----------------------------------------------------------------

# invariant settings for writing files

variable_name = 'pctlake'

# variable attributes
attrs_variable = {'units': '%', 'long_name' : 'Grid cell lake and area fraction', '_FillValue' : 1e20}
attrs_years = {'units': 'years', 'long_name' : 'time'}

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603, GRanD dataset v1.3 2017,  https://doi.org/10.1890/100125', 
                        'title': 'Lake and reservoir percent, transient over time',
                        'use': 'Lake and reservoir area fraction calculated from the HydroLAKES and GRanD database mapped to 0.5 degrees resolution, time steps represent increasing reservoirs in their construction year, as given by GRanD. Reservoirs from 2017 onwards are kept constant',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Lehner, B., C. Reidy Liermann, C. Revenga, C. Vörösmarty, B. Fekete, P. Crouzet, P. Döll, M. Endejan, K. Frenken, J. Magome, C. Nilsson, J.C. Robertson, R. Rodel, N. Sindorf, and D. Wisser. 2011. High-resolution mapping of the world’s reservoirs and dams for sustainable river-flow management. Frontiers in Ecology and the Environment 9 (9): 494-502, https://doi.org/10.1890/100125. Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022).',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }


# -----------------------------------------------------------
# ISIMIP3a histsoc 1850-1900

years_start = 1850
years_end   = 1900

filename_rasters = [directory_waterarea_per_pixel + 'frac_areas_NA_'+str(year)+'.tif' for year in range(years_start,years_end+1)]

filename_netcdf = directory_netcdf + 'Water_area_per_pixel/' + 'pctlake_'+str(years_start)+'_'+str(years_end)+'.nc'


write_netcdf_3d_lakepct(filename_rasters, filename_netcdf, attrs_variable,variable_name, attrs_global, years_start, years_end)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -setcalendar,360_day -setreftime,1850-01-01,00:00:00,1years -settaxis,1850-01-01,00:00:00,1years ' + filename_netcdf + ' ' + filename_netcdf + '4')

os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)

# -----------------------------------------------------------
# ISIMIP3a histsoc 1901-2020


years_start = 1901
years_end   = 2020
years_files_end = 2017

filename_rasters = [directory_waterarea_per_pixel + 'frac_areas_NA_'+str(year)+'.tif' for year in range(years_start,years_files_end+1)]


# extend filenames from 2017 to 2020 by repeating last 3 years
for year in range(years_files_end,years_end): 
    filename_rasters.append(directory_waterarea_per_pixel + 'frac_areas_NA_'+str(years_files_end)+'.tif')

filename_netcdf = directory_netcdf + 'Water_area_per_pixel/' + 'pctlake_'+str(years_start)+'_'+str(years_end)+'.nc'


write_netcdf_3d_lakepct(filename_rasters, filename_netcdf, attrs_variable,variable_name, attrs_global, years_start, years_end)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat -setcalendar,360_day -setreftime,1901-01-01,00:00:00,1years -settaxis,1901-01-01,00:00:00,1years ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)


# -----------------------------------------------------------
# ISIMIP3b histsoc 1850-2014

# define filename, variable name and attributes


years_start = 1850
years_end   = 2014

filename_rasters = [directory_waterarea_per_pixel + 'frac_areas_NA_'+str(year)+'.tif' for year in range(years_start,years_end+1)]


filename_netcdf = directory_netcdf + 'Water_area_per_pixel/' + 'pctlake_'+str(years_start)+'_'+str(years_end)+'.nc'


write_netcdf_3d_lakepct(filename_rasters, filename_netcdf, attrs_variable, variable_name, attrs_global, years_start, years_end)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat -setcalendar,360_day -setreftime,1850-01-01,00:00:00,1years -settaxis,1850-01-01,00:00:00,1years ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)


# %%


# -----------------------------------------------------------
# by Daniel to upload file requested by Matthias Buechner
# ISIMIP3b histsoc 1850-2020

# define filename, variable name and attributes


years_start = 1850
years_end   = 2020
years_files_end = 2017

filename_rasters = [directory_waterarea_per_pixel + 'frac_areas_NA_'+str(year)+'.tif' for year in range(years_start,years_files_end+1)]


# extend filenames from 2017 to 2020 by repeating last 3 years
for year in range(years_files_end,years_end): 
    filename_rasters.append(directory_waterarea_per_pixel + 'frac_areas_NA_'+str(years_files_end)+'.tif')

filename_netcdf = directory_netcdf + 'Water_area_per_pixel/' + 'pctlake_'+str(years_start)+'_'+str(years_end)+'.nc'


write_netcdf_3d_lakepct(filename_rasters, filename_netcdf, attrs_variable,variable_name, attrs_global, years_start, years_end)
os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat -setcalendar,360_day -setreftime,1850-01-01,00:00:00,1years -settaxis,1850-01-01,00:00:00,1years ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)


# %%




# -----------------------------------------------------------
# ISIMIP3b 1850soc


year = 1850


variable_name = 'pctlake'
# variable attributes

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603, GRanD dataset v1.3 2017,  https://doi.org/10.1890/100125', 
                        'title': 'Lake and reservoir percent for 1850 (natural lakes excluding reservoirs)',
                        'use': 'Lake and reservoir area fraction calculated from the HydroLAKES and GRanD database mapped to 0.5 degrees resolution for 1850, mostly including natural lakes',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Lehner, B., C. Reidy Liermann, C. Revenga, C. Vörösmarty, B. Fekete, P. Crouzet, P. Döll, M. Endejan, K. Frenken, J. Magome, C. Nilsson, J.C. Robertson, R. Rodel, N. Sindorf, and D. Wisser. 2011. High-resolution mapping of the world’s reservoirs and dams for sustainable river-flow management. Frontiers in Ecology and the Environment 9 (9): 494-502, https://doi.org/10.1890/100125. Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022).',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }


filename_raster = directory_waterarea_per_pixel + 'frac_areas_NA_'+str(year)+'.tif'


filename_netcdf = directory_netcdf + 'Water_area_per_pixel/' + 'pctlake_'+str(year)+'soc.nc'

write_netcdf_2d_pctlake(filename_raster,  filename_netcdf,attrs_variable,variable_name,attrs_global)

os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat -setcalendar,360_day -setreftime,1850-01-01,00:00:00,1years -settaxis,1850-01-01,00:00:00,1years ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)



# -----------------------------------------------------------
# ISIMIP3b 2015soc


year = 2015


variable_name = 'pctlake'
# variable attributes

# global attributes
attrs_global = {'creation_date': date,
                        'source': 'HydroLAKES dataset v1.0 June 2019, https://doi.org/10.1038/ncomms13603, GRanD dataset v1.3 2017,  https://doi.org/10.1890/100125', 
                        'title': 'Lake and reservoir percent for 2015 ',
                        'use': 'Lake and reservoir area fraction calculated from the HydroLAKES and GRanD database mapped to 0.5 degrees resolution for 1850, mostly including natural lakes',
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Inne Vanderkelen - VUB (inne.vanderkelen@vub.be); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'references':'Lehner, B., C. Reidy Liermann, C. Revenga, C. Vörösmarty, B. Fekete, P. Crouzet, P. Döll, M. Endejan, K. Frenken, J. Magome, C. Nilsson, J.C. Robertson, R. Rodel, N. Sindorf, and D. Wisser. 2011. High-resolution mapping of the world’s reservoirs and dams for sustainable river-flow management. Frontiers in Ecology and the Environment 9 (9): 494-502, https://doi.org/10.1890/100125. Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603, Khazaei, B., Read, L. K., Casali, M., Sampson, K. M., & Yates, D. N. (2022).',
                        'url' : 'https://github.com/icra/ISIMIP_Lake_Sector' }


filename_raster = directory_waterarea_per_pixel + 'frac_areas_NA_'+str(year)+'.tif' 


filename_netcdf = directory_netcdf + 'Water_area_per_pixel/' + 'pctlake_'+str(year)+'soc.nc'

write_netcdf_2d_pctlake(filename_raster,  filename_netcdf,attrs_variable,variable_name,attrs_global)

os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat -setcalendar,360_day -setreftime,2015-01-01,00:00:00,1years -settaxis,2015-01-01,00:00:00,1years ' + filename_netcdf + ' ' + filename_netcdf + '4')
os.system('rm '+ filename_netcdf)
os.system('mv '+ filename_netcdf+ '4' + ' ' + filename_netcdf)



# %%
