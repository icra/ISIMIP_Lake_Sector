Analysis of Hydrolakes and Globathy databases to generate lake
morphometry for ISIMIP3
================
Rafael Marcé and Daniel Mercado,
March 2022

## Antecedents

During the ISIMIP Lake Sector in February 2022 the community decides to
change the strategy to define representative lakes per pixel for ISIMIP3
runs. In a separate prevuious analysis, the representative lakes per
pixel have been identified calculating the weighted median lake depth,
using lake area as weights. In this way, individual lakes have been
selected at each pixel, for a total of 41449 lakes. The selected lakes
correspond to real lakes in the HydroLakes database.

In this part of the analysis, we will define the morphometry of each
lake (hypsographic curve). In ISIMIP2 all lakes had a cylyndrical shape,
and the majority of modellers were in favour of trying to define a more
realistic shape for each representative lake.

Taking advantage of the recent publication of the GLOBATHY database
(<https://www.nature.com/articles/s41597-022-01132-9>), which estimates
maximum depth and bathymetric information for all lakes in Hydrolakes,
we analyze which is the best alternative to represent lake morphometry
in the most efficient way in ISIMIP3 runs.

This task was led by the authors above, with contributions by Inne
Vanderkelen, Maddalena Tigli, Iestyn Woolway, Annette Janssen, Benjamin
Kraemer, Mahtab Yaghouti, Sebastiano Piccolroaz, Wim Thiery, and Don
Pierson

## Opening databases and checking coherence

Loading the attribute table of the HydroLakes polygons (~1.4M lakes).
You will likely have to insert your local path to find the file - too
big to be stored at
Github

``` r
HydroLakes <- read.dbf("/home/rmarce/Cloud/a. WATExR/ISIMIP/Area_depth February 2022/Hydrolakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.dbf")
```

Loading basic morphometric parameters from Khazaei et al.
(2022)

``` r
Globathy_basic <-read.csv("/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv",header=T)
```

Checking order to make sure we do not mess with different lakes, this
sum must be zero

``` r
sum(HydroLakes$Hylak_id-Globathy_basic$Hylak_id)
```

    ## [1] 0

Adding maximum depth (Dmax) from Globathy to HydroLakes. Using the Dmax
result reccommended in Khazaei et al. (2022) from several alternatives

``` r
HydroLakes$Dmax_Khazaei <- Globathy_basic$Dmax_use_m
```

Comparing Dmax in Globathy vs Dmean in HydroLakes. Some impossible
situations (Dmax\<mean depth)

Using ratios:

``` r
ratio_mm <-HydroLakes$Dmax_Khazaei/HydroLakes$Depth_avg 
summary(ratio_mm)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##   0.02705   2.22827   3.24813   3.94063   4.75387 227.74006

``` r
hist(log10(ratio_mm))
```

![](Morphometry_ISIMIP3_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Using differences:

``` r
dif_mm <- HydroLakes$Dmax_Khazaei-HydroLakes$Depth_avg
summary(dif_mm)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -556.330    4.822    7.611    7.570   10.056  903.300

``` r
hist(dif_mm)
```

![](Morphometry_ISIMIP3_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
boxplot(dif_mm)
```

![](Morphometry_ISIMIP3_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

Percent of impossible values using Dmax from Globathy and Dmean from
HydroLakes:

``` r
length(which(dif_mm<0))/length(dif_mm)*100
```

    ## [1] 3.123091

Cheking the volume development parameter (Vd) using Dmax from Globathy
and Dmean from HydroLakes (should between 0 and 3) Note the nonsense
values (\>3) and the median \<1, which is not very realistic considering
classical literature (for instance:
<https://doi.org/10.1016/B978-012370626-3.00024-7>,
<https://link.springer.com/article/10.1007/s10666-006-9069-z>,
<http://www.jstor.org/stable/520579> )

``` r
HydroLakes$Vd  <- 3*HydroLakes$Depth_avg/HydroLakes$Dmax_Khazaei
summary(HydroLakes$Vd)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##   0.01317   0.63107   0.92361   1.11600   1.34633 110.89223

``` r
hist(log10(HydroLakes$Vd)) 
```

![](Morphometry_ISIMIP3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
hist((HydroLakes$Vd),xlim=c(0,3),breaks=400)
```

![](Morphometry_ISIMIP3_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

Loading extended morphological information from Globathy (NetCDF
file)

``` r
nc_data <- nc_open('/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_hAV_relationships.nc')
# Save the print(nc) dump to a text file
{
  sink('gimms3g_ndvi_1982-2012_metadata.txt')
  print(nc_data)
  sink()
}
```

Extracting the lake basic attributes: max depth, mean depth, surface
area, total volume (m, m, km2, km3)

``` r
attri_nc <- ncvar_get(nc_data, "lake_attributes")
attri_nc[1:4,1:3]#just to look at aspect
```

    ##             [,1]       [,2]       [,3]
    ## [1,]   1025.0000   446.0000   614.0000
    ## [2,]    317.5679   108.5791   133.3336
    ## [3,] 376904.5992 30574.0173 26664.3815
    ## [4,] 119692.8818  3319.7045  3555.2596

Collecting maximum volume and area, and maximum and mean depth

``` r
A_Khazaei <- attri_nc[3,] 
V_Khazaei <- attri_nc[4,]
Dmax_Khazaei <- attri_nc[1,] 
Dmean_Khazaei <-attri_nc[2,] 
```

\#comparing Hydrolakes and Khazaei volumes. Volume in hydrolakes is in
0.001 km3
plot(data\(Vol_total*0.001,V_Khazaei,log="xy") abline(0,1) #comparing Hydrolakes and Khazaei areas plot(data\)Lake\_area,A\_Khazaei,log=“xy”)
abline(0,1) \#comparing Hydrolakes and Khazaei mean depth
plot(data$Depth\_avg,Dmean\_Khazaei,log=“xy”) abline(0,1)

\#difference in total volum and area, in percent
(sum(V\_Khazaei)-sum(data\(Vol_total)*0.001)/(sum(data\)Vol\_total)*0.001)*100
(sum(A\_Khazaei)-sum(data\(Lake_area))/(sum(data\)Lake\_area))\*100

\#distribution of individual differences, in percent dif\_area \<-
(data\(Lake_area-A_Khazaei)/data\)Lake\_area\*100 hist(dif\_area)

dif\_volume \<-
(data\(Vol_total*0.001-V_Khazaei)/(data\)Vol\_total*0.001)*100
hist(dif\_volume)

dif\_Dmean \<- (data\(Depth_avg-Dmean_Khazaei)/data\)Depth\_avg\*100
hist(dif\_Dmean)

# volume development using data from Khazaei exclusively

Vd\_Khazaei \<- 3\*Dmean\_Khazaei/Dmax\_Khazaei hist((Vd\_Khazaei))
summary(Vd\_Khazaei)

\#ANalyzing the results for the selected representative lakes
\#recovering lakes ID from Khazaei databse nc\_lake\_ID \<-
ncvar\_get(nc\_data, “lake\_id”) dim(nc\_lake\_ID)

\#loading the selected representative lakes for ISIMIP3
data\_selected\<-read.dbf(“/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP
Lake Sector Feb 2022 - ICRA PC/HL\_selected.dbf”)
head(data\_selected$Hylak\_id)

\#attributes of the selected representative lakes  
attri\_nc\_selected \<- attri\_nc\[,data\_selected$Hylak\_id\]
dim(attri\_nc\_selected) attri\_nc\_selected\[,1:10\]

\#Atttributes to variables A\_Khazaei\_selected \<-
attri\_nc\_selected\[3,\] V\_Khazaei\_selected \<-
attri\_nc\_selected\[4,\] Dmax\_Khazaei\_selected \<-
attri\_nc\_selected\[1,\] Dmean\_Khazaei\_selected
\<-attri\_nc\_selected\[2,\]

\#Volume development for selected lakes Vd\_Khazaei\_selected \<-
3\*Dmean\_Khazaei\_selected/Dmax\_Khazaei\_selected
hist((Vd\_Khazaei\_selected))

\#comparing with whole-database results summary(Vd\_Khazaei\_selected)
summary(Vd\_Khazaei)

plot(density(Vd\_Khazaei), col=“red”)
lines(density(Vd\_Khazaei\_selected))
plot(quantile(Vd\_Khazaei,p=seq(0,1,0.01)),type=“l”)
lines(quantile(Vd\_Khazaei\_selected,p=seq(0,1,0.01)), col=“red”)

\#comparing volumes and areas with Hydrolakes for selected
representative lakes
plot(data\(Vol_total[data_selected\)Hylak\_id\]\*0.001,V\_Khazaei\_selected,log=“xy”)
abline(0,1)
plot(data\(Lake_area[data_selected\)Hylak\_id\],A\_Khazaei\_selected,log=“xy”)
abline(0,1)

\#save.image(“.RData”)

plot(Dmax\_Khazaei\_selected,Vd\_Khazaei\_selected,log=“x”)
plot(A\_Khazaei\_selected,Vd\_Khazaei\_selected,log=“x”)

\#saving rasters for ISIMIP3

\#loading lake\_ID of representative lakes for reference to assign
values to rasters raster\_id \<-
raster(“/home/rmarce/ISIMIP\_Lake\_Sector/output/Hylak\_id.tif”)
raster\_id\_asmatrix \<- as.matrix(raster\_id)\#to speed up

\#initializing matrices that will eventually become the rasters
matrix\_vd \<- matrix(data = NA, nrow=360, ncol = 720) \#empty matrix
matrix\_A \<- matrix(data = NA, nrow=360, ncol = 720) \#empty matrix
matrix\_V \<- matrix(data = NA, nrow=360, ncol = 720) \#empty matrix
matrix\_Dmax \<- matrix(data = NA, nrow=360, ncol = 720) \#empty matrix
matrix\_Dmean \<- matrix(data = NA, nrow=360, ncol = 720) \#empty matrix

\#loop for assigning values in the matrices for (i in 1:360){ for (j in
1:720){ \#print(i);print(j) if (is.na(raster\_id\_asmatrix\[i,j\]))
{\#to avoid the loop crashing }else{ \#identifying which lake goes to
the current position position \<- which(data\_selected$Hylak\_id ==
raster\_id\_asmatrix\[i,j\]) matrix\_vd\[i,j\] \<-
Vd\_Khazaei\_selected\[position\] matrix\_A\[i,j\] \<-
A\_Khazaei\_selected\[position\] matrix\_V\[i,j\] \<-
V\_Khazaei\_selected\[position\] matrix\_Dmax\[i,j\] \<-
Dmax\_Khazaei\_selected\[position\] matrix\_Dmean\[i,j\] \<-
Dmean\_Khazaei\_selected\[position\] } } }

## converting to raster and writing data

Vd\_raster \<- raster(matrix\_vd) A\_raster \<- raster(matrix\_A)
V\_raster \<- raster(matrix\_V) Dmax\_raster \<- raster(matrix\_Dmax)
Dmean\_raster \<- raster(matrix\_Dmean)

extent(Vd\_raster) \<- extent(c(-180,180,-90,90)) extent(A\_raster) \<-
extent(c(-180,180,-90,90)) extent(V\_raster) \<-
extent(c(-180,180,-90,90)) extent(Dmax\_raster) \<-
extent(c(-180,180,-90,90)) extent(Dmean\_raster) \<-
extent(c(-180,180,-90,90))

crs(Vd\_raster) \<- “+proj=longlat +datum=WGS84 +no\_defs +ellps=WGS84
+towgs84=0,0,0” crs(A\_raster) \<- “+proj=longlat +datum=WGS84 +no\_defs
+ellps=WGS84 +towgs84=0,0,0” crs(V\_raster) \<- “+proj=longlat
+datum=WGS84 +no\_defs +ellps=WGS84 +towgs84=0,0,0” crs(Dmax\_raster)
\<- “+proj=longlat +datum=WGS84 +no\_defs +ellps=WGS84 +towgs84=0,0,0”
crs(Dmean\_raster) \<- “+proj=longlat +datum=WGS84 +no\_defs
+ellps=WGS84 +towgs84=0,0,0”

writeRaster(Vd\_raster,“./Results/Vd\_raster.tif”, overwrite=T)
writeRaster(A\_raster,“./Results/A\_raster.tif”, overwrite=T)
writeRaster(V\_raster,“./Results/V\_raster.tif”, overwrite=T)
writeRaster(Dmax\_raster,“./Results/Dmax\_raster.tif”, overwrite=T)
writeRaster(Dmean\_raster,“./Results/Dmean\_raster.tif”, overwrite=T)

plot(Vd\_raster,col=viridis(2)) plot(log10(Dmax\_raster),col=viridis(2))

aa\<-(V\_raster/(A\_raster*Dmean\_raster*0.001)) summary(aa) hist(aa)
\#looks like rounding errors

bb\<-(Vd\_raster/(3\*Dmean\_raster/Dmax\_raster)) summary(bb)

\#Collating hypsografics

Extraxting the depth index of the hypsografics (11 depths per lake) id
\<- ncvar\_get(nc\_data, “h”) (id\[1:11,1:11\])

V \<- ncvar\_get(nc\_data, “V”) V\[1:11,1:3\] A \<- ncvar\_get(nc\_data,
“A”) A\[1:11,1:3\] You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](Morphometry_ISIMIP3_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
