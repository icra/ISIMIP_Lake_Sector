# Inputs for ISIMIP3 global lake sector <br />
## Explanation of the method implemented to obtain the inputs for ISIMIP3 simulations <br />

### This is a brief guide to follow the steps implemented to obtain the results, if you want more details about please check the three folders containing the full calculations and processes:<br />

**[lake_identification]:** the first calculation to select the representative lakes for each pixel. We took the 1.4 million lakes in HydroLAKES and calculated the depth weighted median (weighted by area of the lakes) for all the lakes contained in each pixel with a 0.5º resolution convering -180 to 180 longitude and -90 to 90 latitude.  <br />
**[Hypsographics]:** with the selected representative lakes for each pixel in the previous step, we extract the Volume, Area, mean and maximum Depth, and bathymetry from GLOBathy (https://www.nature.com/articles/s41597-022-01132-9) for each representative lake selected in the previous step using id (GLOBathy and HydroLAKES use the same id attribute). 
**[Water_area_per_pixel]:**: it's stil under development. This will be usefull to implement future analysis with the final data from ISIMIP3, so the simulation round doesn't depend on it. <br />

### The final results containing the inputs for ISIMIP3 simulation <br />
**[inputs_ISIMIP3]:** final results from processing, ready to use for ISIMIP3 simulations <br />
