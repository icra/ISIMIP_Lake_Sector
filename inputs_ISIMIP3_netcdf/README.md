# Final collection of netcdf files containing the representative lakes information for ISIMIP3 global lake runs  <br />

**[Basic_inputs]:** Area (km2), volume (km3), and maximum and mean depth (m) of the representative lakes. <br />

**[Hypsographic_curves]:** Hypsographic curves for area and volume for the representative lakes. Each hypsographic curve consist of 11 data pairs relating level (bottom of the lake as the reference, in meters) with area (km2) or volume (km3).  <br />

**[Power_fit_to_hypsographic_curves]:** Power law fit of the hypsographic curves in [Hypsographic_curves], defining the relationship between level (bottom of the lake as the reference, in m) and area (km2) or volume (km3). The power fits are of the form Volume=a*h^b or Area=a*h^b, where h is level and a and b are fitted parameters. The netcdf files store the values of the two parameters for each representative lake, and also a goodness of fit indication (R2). Those fits are generally good, but deviations from the reference hypsographic can be large in some cases. Therefore, we reccommend using the hypsographic curve in [Hypsographic_curves] whenever possible.  <br />

**[Reference]:** netcdf files containing reference information for the representative lakes.


