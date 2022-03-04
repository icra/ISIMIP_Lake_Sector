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
data <- read.dbf("/home/rmarce/Cloud/a. WATExR/ISIMIP/Area_depth February 2022/Hydrolakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.dbf")
```

Loading basic morphometric parameters from Khazaei et al.
(2022)

``` r
Dmax <-read.csv("/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv",header=T)
```

Checking order to make sure we do not mess with different lakes, this
sum must be zero

``` r
sum(data$Hylak_id-Dmax$Hylak_id)
```

    ## [1] 0

\#adding Dmax to Hydrolakes data\(Dmax_Khazaei <- Dmax\)Dmax\_use\_m
\#this is the one that paper reccomends
data\(Dmax_Khazaei_cone <- Dmax\)Dmax\_cone\_m \#assuming a cone,
provided in that same paper

\#comparing Dmax\_Khazaei vs mean depth in Hydrolakes. Some impossible
situations. \#with ratios ratio\_mm
\<-data\(Dmax_Khazaei/data\)Depth\_avg summary(ratio\_mm)
hist(log10(ratio\_mm))

\#with differences dif\_mm \<- data\(Dmax_Khazaei-data\)Depth\_avg
summary(dif\_mm) hist(dif\_mm) boxplot(dif\_mm)

\#now with the cone assumption, better behaviour dif\_mm\_cone \<-
data\(Dmax_Khazaei_cone-data\)Depth\_avg summary(dif\_mm\_cone)
hist(dif\_mm\_cone)

\#percent of impossible values with Dmax\_Khazaei
length(which(dif\_mm\<0))/length(dif\_mm)\*100

\#cheking the volume development using Dmax\_Khazaei (should between 0
and 3)
data\(Vd <- 3*data\)Depth\_avg/data\(Dmax_Khazaei summary(data\)Vd)
hist(log10(data\(Vd)) hist((data\)Vd),xlim=c(0,3),breaks=400) \#nonsense
values (\>3) and median \<1, not very realistic distribution \# \#
\#volume of a cone \# pi*R<sup>2*h/3 \# \#volume of a truncated cone \#
pi*h*(r^2+r*R+R</sup>2)/3 \# \# \#volume \# Area*mean\_depth \#
pi*R<sup>2*mean\_depth \# \# pi*R</sup>2*mean\_depth \# pi*R^2*h/3 \# \#
mean\_depth \# h/3 \# \# 3\*mean\_depth/h \#

# Testing hypsografics in Khazaei

nc\_data \<- nc\_open(‘/home/rmarce/Cloud/a. WATExR/ISIMIP/ISIMIP Lake
Sector Feb 2022 - ICRA PC/GLOBathy/GLOBathy\_hAV\_relationships.nc’) \#
Save the print(nc) dump to a text file {
sink(‘gimms3g\_ndvi\_1982-2012\_metadata.txt’) print(nc\_data) sink()
} \#load the depth index of the hypsografics (11 depths per) id \<-
ncvar\_get(nc\_data, “h”) (id\[1:11,1:11\])

\#load the lake attributes:max depth, mean depth, surface area, total
volume (m, m, km2, km3) attri\_nc \<- ncvar\_get(nc\_data,
“lake\_attributes”) attri\_nc\[1:4,1:3\]

\#Collecting maximum volume and area, max and mean depth A\_Khazaei \<-
attri\_nc\[3,\] V\_Khazaei \<- attri\_nc\[4,\] Dmax\_Khazaei \<-
attri\_nc\[1,\] Dmean\_Khazaei \<-attri\_nc\[2,\]

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

![](xxx_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
