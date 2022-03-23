## Input data to use ISIMIP3 simulations <br />

**[lake_identification]** files that are outcomes from the first step (lake_identification folder in main) calculated from HydroLAKES <br />
**[Hypsographics]** files that are outcomes from the second step (Hypsographics folder in main) extracted from GLOBathy<br />

**Hylak_id.tif [lake_identification]:** unique lake identifier<br />
**Lake_type.tif [lake_identification]:** indicator for lake type; 1:Lake; 2:Reservoir; 3:Lake control (i.e. natural lake with regulation structure)<br />
**HL_selected [lake_identification]:** shapefile of the final selected lakes for each pixel. Please find the HL_selected.shp file in lake_identification/output<br />
**HLid_biglakes.tif [lake_identification]:** ID of big lakes, i.e., with area greater than 0.5 degrees <br />
**biglakes_mask [lake_identification]:**  mask for postprocessing and plotting big lakes (> 0.5º) <br />

**A_raster.tif [Hypsographics]:** lake surface area, in square kilometers <br />
**Dmax_raster [Hypsographics]:** maximum lake depth, in meters <br />
**Dmean_raster [Hypsographics]:** mean lake depth, in meters <br />
**Vd_raster.tif [Hypsographics]:** Total lake volume, in million cubic meters (1 mcm = 0.001 km3)<br />
