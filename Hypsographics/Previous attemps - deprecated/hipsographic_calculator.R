

 volume <- 168555237 #m3
 surface_area <- 5999766.8 #m2
 max_depth <- 59.54 #m
 mean_depth <- volume/surface_area

 sau <- read.csv("sau_hipso.csv",header=F)
 
 hipsographicR_Hd <-  function(volume,surface_area,max_depth,mean_depth){ 
 
   if (max_depth > floor(max_depth)){
     hipsographic <- matrix(NA, floor(max_depth)+2,3) 
   } else {
     hipsographic <- matrix(NA, floor(max_depth)+1,3) 
   }
   
   
   z_m <- mean_depth
   Vd <- 3*z_m/max_depth
   b <- log10(Vd)
   #if per les diferents possibilitats
   Hd <- 10^(-34*b^5-18*b^4-6.3*b^3-1.9*b^2-4*b)
   
   Hd <- 10^(-19261*b^5+27179*b^4-15506*b^3+4432.6*b^2-639.51*b+36.442)
   
   serie=sau[,1]
   
   if (floor(max_depth)>=1){
     #area across slices of the truncated cone
     for (i in 0:floor(max_depth)){   
       
       hipsographic[i+1,1]=i
       hipsographic[i+1,2]=surface_area*(  ( Hd^(-i/max_depth) - Hd^(-max_depth/max_depth) ) / (Hd^(-0/max_depth) - Hd^(-max_depth/max_depth) ) )^2
       hipsographic[i+1,3]= (surface_area/((1-(Hd^-1))^2)) * (  (2*max_depth*(Hd^(-i/max_depth-1)))/log(Hd) - (max_depth*(Hd^(-2*i/max_depth)))/(2*log(Hd)) + (Hd^-2)*i )     
     }
     if (max_depth>floor(max_depth)){
       hipsographic[i+2,1]= max_depth
       hipsographic[i+2,2]= surface_area*(  ( Hd^(-max_depth/max_depth) - Hd^(-max_depth/max_depth) ) / (Hd^(-0/max_depth) - Hd^(-max_depth/max_depth) ) )^2
       hipsographic[i+2,3]= (surface_area/((1-(Hd^-1))^2)) * (  (2*max_depth*(Hd^(-max_depth/max_depth-1)))/log(Hd) - (max_depth*(Hd^(-2*max_depth/max_depth)))/(2*log(Hd)) + (Hd^-2)* max_depth) 
     }
     
   } else{
     
     if (max_depth>floor(max_depth)){
       hipsographic[2,1]= max_depth
       hipsographic[2,2]= surface_area*(  ( Hd^(-max_depth/max_depth) - Hd^(-max_depth/max_depth) ) / (Hd^(-0/max_depth) - Hd^(-max_depth/max_depth) ) )^2
       hipsographic[2,3]= 0
     }
     
     
   }
   
 }
 
 
plot(sau[,2],sau[,2],type="l") 
points(sau[,2],hipso_TC[,2],col="green") 
points(sau[,2],hipsographic[,2],col="red")
summary(lm(sau[,2]~hipso_TC[,2]))
summary(lm(sau[,2]~hipsographic[,2]))

hipso_TC <-  hipsographicR_TruncatedCone(volume,surface_area,max_depth)
lines(hipso_TC[,2],hipso_TC[,1],type="l",ylim = rev(range(hipso_TC[,1]))) 

plot(hipsographic[,2],hipsographic[,1], type="l",col="red",ylim = rev(range(hipsographic[,1])),log="x")
points(sau[,2],sau[,1], col="green",ylim = rev(range(sau[,1])))
 
library(foreign)

#importing the dbf with the lake information
data<-read.dbf("/home/rmarce/Cloud/a. WATExR/ISIMIP/Area_depth February 2022/Hydrolakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.dbf")

#Assumption to go from mean depth to maximum depth
data$max_depth <- data$Depth_avg/0.47

#preparing a list vector to store the results (hypsografic curves)
hipso_list <- vector(mode = "list", length = dim(data)[1])
#hipso_list <- matrix(NA, dim(data)[1],1)

#calculating the hypsografic for every lake in the shapefile
#volume and area must be converted to m^3 and m^2 respectively
for (i in 1:dim(data)[1]){
  hipso_list[i]<-hipsographicR_TruncatedCone(data$Vol_total[i]*1000000,data$Lake_area[i]*1000000,data$max_depth[i])
  #print(hipso_list[i])
  #print(i)
}

# Function for the hypsographic curve assuming a truncated cone from volume (m^3), surface area (m^2), and maximum depth (m).
# Calculated every meter from surface and at maximum depth. 
hipsographicR_TruncatedCone <-  function(volume,surface_area,max_depth){
  
  #volume=total volume of lake (m^3)
  #surface_area= surface area (m^2)
  #max_depth=maximum depth (m)
  
  #initializing variable  
  ##fields for depth,area,accumulated volume
  
  #conditional for including the maximum depth in case it is not an integer depth
  if (max_depth > floor(max_depth)){
    hipsographic <- matrix(NA, floor(max_depth)+2,3) 
  } else {
    hipsographic <- matrix(NA, floor(max_depth)+1,3) 
  }
  
  #solving a truncated cone
  Radius <- sqrt(surface_area/pi) #radius of the large end of the truncated cone (surface of the lake) 
  term <- volume/max_depth/pi*3-Radius^2 #solving the truncated cone equation [=0])
  radius_imaginary <- polyroot(c(-term,Radius,1)) #finding the roots of the equation
  real<-Re(radius_imaginary) # taking the real part of the result
  #the positive real result is the result we want:
  radius <-real[which(real>0)] #radius of the small end of the truncated cone (bottom of the lake)
  
  #checks if the solution makes physical sense and print an error otherwise:
  
  #the root must have a real part
  if (length(radius)==0){
    a=2
    print("Error 2")
    return(a)
  }
  stopifnot(length(radius)>0)
  
  #surface radius must be larger or equal than bottom radius
  if (Radius<radius){
    a=1
    print("Error 1")
    return(a)
  }
  stopifnot(Radius>=radius)
  
  #Surface area
  area_layer <- pi*Radius^2 #must be equal to surface_area
  
  #Result at the surface
  hipsographic[1,1]=0
  hipsographic[1,2]=area_layer 
  hipsographic[1,3]=volume
  
  slope <- (Radius-radius)/max_depth #for calculating the radius at each depth
  
  if (floor(max_depth)>=1){
    #area across slices of the truncated cone
    for (i in 1:floor(max_depth)){   
      
      new_radius <- Radius-i*slope
      area_layer <- pi*new_radius^2  
      hipsographic[i+1,1]=i
      hipsographic[i+1,2]=area_layer
      hipsographic[i+1,3]=pi/3*(max_depth-i)*(new_radius^2+new_radius*radius+radius^2)
    }
    
    if (max_depth>floor(max_depth)){
      hipsographic[i+2,1]= max_depth
      hipsographic[i+2,2]= pi*radius^2
      hipsographic[i+2,3]= 0
    }
    
  } else{
    
    if (max_depth>floor(max_depth)){
      hipsographic[2,1]= max_depth
      hipsographic[2,2]= pi*radius^2
      hipsographic[2,3]= 0
    }
    
    
  }
  
  #a=0
  #return(a)
  return(hipsographic)
}
