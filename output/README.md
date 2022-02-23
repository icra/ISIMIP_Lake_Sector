
## Variables from HydroLAKES <br />
(more detail in manual of the same database)<br />
**Hylak_id [ISIMIP]:** Unique lake identifier<br />
**Depth_avg.tif [ISIMIP]:**  Average lake depth, in meters<br />
**Elevation.tif:** Elevation of lake surface, in meters above sea level<br />
**Lake_area.tif [ISIMIP]:** Lake surface area (i.e. polygon area), in square kilometers<br />
**Res_time.tif:** Average residence time of the lake water, in days<br />
**Shore_len.tif:** Length of shoreline (i.e. polygon outline), in kilometers<br />
**Vol_res.tif:** Reported reservoir volume, or storage volume of added lake regulation, in million cubic meters (1 mcm = 0.001 km 3 ); 0: no reservoir volume<br />
**Vol_total.tif [ISIMIP]:** Total lake or reservoir volume, in million cubic meters (1 mcm = 0.001 km3)<br />
**Dis_avg.tif:** Average long-term discharge flowing through the lake, in cubic meters per second<br />
**Grand_id.tif:** ID of the corresponding reservoir in the GRanD database, or value 0 for no corresponding GRanD record<br />
**Lake_type.tif:** Indicator for lake type; 1:Lake; 2:Reservoir; 3:Lake control (i.e. natural lake with regulation structure)<br />
**Shore_dev.tif:** Shoreline development, measured as the ratio between shoreline length and the circumference of a circle with the same area<br />
**Slope_100.tif:** Average slope within a 100 meter buffer around the lake polygon, in degrees<br />
**Vol_src.tif:** 1: ‘Vol_total’ is the reported total lake volume from literature; 2: ‘Vol_total’ is the reported total reservoir volume from GRanD or literature; 3: ‘Vol_total’ is the estimated total lake volume using the geostatistical modeling approach by Messager et al. (2016)<br />
**Wshd_area.tif:** Area of the watershed associated with the lake, in square kilometers<br />

## Other variables<br />
**HLid_biglakes.tif:** ID of big lakes, i.e., with area greater than 0.5 degrees <br />
**HL_selected.shp:** shapefile of the final selected lakes for each pixel<br />
**HL_big.shp:** shapefile containing the big lakes (with area great than 0.5 degrees), this can be used as mask for postprocessing and plotting<br />
**HL_big_raster.tif:** outcome of rasterising HL_big.shp />
**HL_big_raster_total.tif [ISIMIP]:** outcome of merging HL_big_raster.tif and HLid_biglakes.tif. This merge was implemented because when doing the rasterisation of HL_big.shp some of the lakes were not represented in the final raster />
**HL_cent.shp** file is not in this folder due to its size, but it could be requested to the authors if needed />

## Plots<br />
**area_density.pdf:** density plot of area values of: HydroLAKES (black, 1.4 million of lakes), representative lakes after weighted median (blue, 41000 pixels) and representative lakes after median (red, 41000 pixels)<br />
**depth_density.pdf:** density plot of depth values of: HydroLAKES (black, 1.4 million of lakes), representative lakes after weighted median (blue, 41000 pixels) and representative lakes after median (red, 41000 pixels)<br />

**[ISIMIP]** files with this are part of ISIMIP3 input data for global lake sector
