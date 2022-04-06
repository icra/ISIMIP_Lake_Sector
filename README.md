# Inputs for lake morphology in ISIMIP3 runs of the global lake sector <br />
## Method implemented to obtain the lake collection and morphology inputs for ISIMIP3 simulations in global lakes <br />

### This is a brief guide to follow the steps implemented, if you want more details please check the three folders containing the full calculations:<br />

**[lake_identification]:** selection of the representative lake for each pixel. We took the 1.4 million lakes in [HydroLAKES](https://www.hydrosheds.org/pages/hydrolakes) and calculated the depth weighted median (weighted by area of the lakes) for all the lakes contained in each pixel with a 0.5º resolution.  <br />
**[Hypsographics]:** we extract the Volume, Area, mean and maximum Depth, and hypsographic curves from [GLOBathy](https://www.nature.com/articles/s41597-022-01132-9) (find here [repository](https://springernature.figshare.com/collections/GLOBathy_the_Global_Lakes_Bathymetry_Dataset/5243309)) for each representative lake selected in the previous step, using the lake ID to cross reference the two databases (GLOBathy and HydroLAKES use the same id attribute). <br />
**[create_netcdf]:** Phyton scripts to produce the final netcdf files available in [inputs_ISIMIP_netcdf]

**[Water_area_per_pixel]:**: Under development. Will contain the fraction of ecah pixel covered by lakes present in the HydroLAKES database. <br />

### Folders containing the lake morphology inputs for ISIMIP3 sumulations in the global lakes sector <br />

**[inputs_ISIMIP3_netcdf]:** Final collection of inputs for the ISIMIP3 runs of the global lake sector in netcdf format. THESE ARE THE RECCOMENDED FILES TO USE. See the Readme.md file in the folder for further details.

**[inputs_ISIMIP3]:** Final input datasets in raster and shapefile format. In this folder you will find information in raster and shapefile format on the representative lakes per pixel, basic morphological features for each lake, and hypsographic curves. PLease, consider there are no supporting metadata for these files, so we reccommend using the netcdf collection at [inputs_ISIMIP3_netcdf] <br />

## Summary of the approach

The representative lakes is a collection of lakes to be used for global lakes ISIMIP3 runs. There is one lake per 0.5º pixel, and the final selection was the result of a procedure that tended to pick conspicuous lakes in the landscape (in terms of both area and depth). You will find an assessment of the implications of our choice in the folder [Hypsographics/GLOBathy]. The final set considers 41449 lakes, one per pixel. Very big lakes occupying more than one pixel are associated to one single pixel as well, but the area and volume information stored in that pixel refers to the whole lake. We have included two files (biglakes_mask.nc, hydrolakes_id_biglakes.nc) in case you want to plot the results of very big lakes on their real extension, not just one pixel. We have included those files to help producing plots that would look realistic (i.e., with very big lakes conspicuous in a global or regional map), but we urge modellers to use the original approach (1 lake, 1 pixel) for model runs and statistical analyses of the results.

You can identify each representative lake with the hydrolakes_id.nc. Since this variable stores the Lake ID used in HydroLAKES, all attributes in this databse can be easily recovered. Please, be mindful that representative lake definition in ISIMIP3 selected real lakes included in HydroLAKES, therefore the representative lakes in ISIMIP3 ARE NOT a statistical construct, but an actual lake in the landscape.    

For the morphology of each representative lake, we included information to allow calculations considering different sources of information. Modellers can pick the method more suited to the model to be applied:

1. Hysographic curve: this must be the default method if your model accepts an hypsographic curve as an input. We have included a level-area and a level-volume hypsographic curve for each of the 41449 lakes. Each hypsographic curve consist of 11 data pairs relating level (bottom of the lake as the reference) with area or volume. This is stored in hypso_area.nc and hypso_volume.nc. 
 
2. If an hypsographic curve cannot be input to your model, we also offer information on lake total volume, maximum area, and mean and maximum depth. You will need an assumption of lake shape (cyclyndrical, conical, etc.) to use this approach. You will find the corresponding netcdf files in [ISIMIP_Lake_Sector/inputs_ISIMIP3_netcdf/Basic_inputs] 

