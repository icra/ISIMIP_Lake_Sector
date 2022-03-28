# Inputs for lake morphology in ISIMIP3 runs of the global lake sector <br />
## Method implemented to obtain the lake collection and morphology inputs for ISIMIP3 simulations in global lakes <br />

### This is a brief guide to follow the steps implemented to obtain the results, if you want more details please check the three folders containing the full calculations:<br />

**[lake_identification]:** selection of the representative lake for each pixel. We took the 1.4 million lakes in [HydroLAKES](https://www.hydrosheds.org/pages/hydrolakes) and calculated the depth weighted median (weighted by area of the lakes) for all the lakes contained in each pixel with a 0.5º resolution.  <br />
**[Hypsographics]:** we extract the Volume, Area, mean and maximum Depth, and hypsographic curves from [GLOBathy](https://www.nature.com/articles/s41597-022-01132-9) (find here [repository](https://springernature.figshare.com/collections/GLOBathy_the_Global_Lakes_Bathymetry_Dataset/5243309)) for each representative lake selected in the previous step, using the lake ID to cross reference the two databases (GLOBathy and HydroLAKES use the same id attribute). <br />
**[Water_area_per_pixel]:**: Under development. Will contain the fraction of ecah pixel covered by lakes present in the HydroLAKES database. <br />

### Final results containing the lake morphology inputs for ISIMIP3 sumulations in the global lakes sector <br />
**[inputs_ISIMIP3]:** final datasets ready to use for ISIMIP3 simulations <br />

In this folder you will find information in raster and shapefile format on the representative lakes per pixel, and basic morphological features for each lake.

The representative lakes is the collection of lakes to be used for global lakes ISIMIP3 runs. There is one lake per 0.5º pixel, and the final selection was the result of a procedure that tended to pick conspicuous lakes in the landscape (in terms of both area and depth). You will find an assessment of the implications of our choice in the folder [lake_identification]. The final set considers 41449 lakes. 

For the morphology of eah representative lake, we included information to allow calculations considering different sources of information. Modellers can pick the method more suited to the model to be applied:

1. Hysographic curve: this must be the default method if your model accepts an hypsographic curve as an input. We have included a level-area and a level-volume hypsographic curve for each of the 41449 lakes. Each hypsographic curve consist of 11 data pairs relating level (bottom of the lake as the reference) with area or volume. 
 
2. We also included a power law fit of the previous hypsographic curves, defining the relationship between level (bottom of the lake as the reference) and area or volume. Those fits are generally good, but deviations from the reference hypsographic can be large in some cases. We reccommend using the hypsographic curve whenever possible.

3. If an hypsographic curve cannot be input to your model, we also offer information on lake total volume, maximum area, and mean and maximum depth. You will need an assumption of lake shape (cyclyndrical, conical, etc.) to use this approach.

