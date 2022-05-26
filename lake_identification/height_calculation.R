library(raster); library(foreign)

HL_id <- raster("/home/ry4902/ISIMIP_Lake_Sector/inputs_ISIMIP3/Hylak_id.tif")
HL_id_NoNA <- HL_id[][!is.na(HL_id[])]
HL_sel <- read.dbf("/home/ry4902/ISIMIP_Lake_Sector/inputs_ISIMIP3/HL_selected.dbf")
HL_ele <- HL_id
for (l in 1:length(HL_id_NoNA)){
  pos_HL_sel <- which(HL_sel$Hylak_id == HL_id_NoNA[l])
  pos_HL_ele <- which(HL_id[] == HL_id_NoNA[l])
  HL_ele[pos_HL_ele] <- HL_sel$Elevation[pos_HL_sel]
  print(l)
}

writeRaster(HL_ele, "/home/ry4902/ISIMIP_Lake_Sector/lake_identification/output/Height.tif")
writeRaster(HL_ele, "/home/ry4902/ISIMIP_Lake_Sector/inputs_ISIMIP3/H_raster.tif")