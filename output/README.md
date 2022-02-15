
## Variables from HydroLAKES (more detail in manual of the same database<br />
**Hylak_id:** Unique lake identifier<br />
**Depth_avg.tif:**  Average lake depth, in meters<br />
**Elevation.tif:** Elevation of lake surface, in meters above sea level<br />
**Lake_area.tif:** Lake surface area (i.e. polygon area), in square kilometers
**Res_time.tif:** Average residence time of the lake water, in days
**Shore_len.tif:** Length of shoreline (i.e. polygon outline), in kilometers
**Vol_res.tif:** Reported reservoir volume, or storage volume of added lake regulation, in million cubic meters (1 mcm = 0.001 km 3 ); 0: no reservoir volume
**Vol_total.tif:** Total lake or reservoir volume, in million cubic meters (1 mcm = 0.001 km3)
**Dis_avg.tif:** Average long-term discharge flowing through the lake, in cubic meters per second
**Grand_id.tif:** ID of the corresponding reservoir in the GRanD database, or value 0 for no corresponding GRanD record.
**Lake_type.tif:** Indicator for lake type; 1:Lake; 2:Reservoir; 3:Lake control (i.e. natural lake with regulation structure)
**Shore_dev.tif:** Shoreline development, measured as the ratio between shoreline length and the circumference of a circle with the same area.
**Slope_100.tif:** Average slope within a 100 meter buffer around the lake polygon, in degrees
**Vol_src.tif:** 1: ‘Vol_total’ is the reported total lake volume from literature; 2: ‘Vol_total’ is the reported total reservoir volume from GRanD or literature; 3: ‘Vol_total’ is the estimated total lake volume using the geostatistical modeling approach by Messager et al. (2016)
**Wshd_area.tif:**

##Other variables
**HLid_biglakes.tif:** ID of big lakes, i.e., with area great than 0.5 degrees
**HL_selected.shp:** shapefile of the final selected lakes for each pixel
**HL_big.shp:** shapefile containing the big lakes (with area great than 0.5 degrees), this can be used as mask for postprocessing and plotting
**HL_cent.shp** file is not in this folder due to its size, but it could be requested to the authors if needed

##Plots
**area_density.pdf:** density plot of area values of: HydroLAKES (black, 1.4 million of lakes), representative lakes after weighted median (blue, 41000 pixels) and representative lakes after median (red, 41000 pixels)
**depth_density.pdf:** density plot of depth values of: HydroLAKES (black, 1.4 million of lakes), representative lakes after weighted median (blue, 41000 pixels) and representative lakes after median (red, 41000 pixels)
