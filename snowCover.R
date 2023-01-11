suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(gdalUtils))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(sf))

print("Starting Workflow")

##### Input files:

# Sentinel-2 L2A (mandatory)
# FMASK output (mandatory)
# DEM (mandatory)
# AOI (optional)

###### Set working directory

repo <- getwd()
print(repo)
repo_data <- paste0(repo,"/data")
setwd(repo_data)

###### Read AOI (if present) and DEM

aoizip <- file.exists(pattern=glob2rx("*.zip"))
aoifile <- FALSE

if (aoizip == TRUE){
  aoifile <- TRUE
  unzip("aoi.zip", files = NULL, list = FALSE, overwrite = TRUE,
        junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)
  aoi <- shapefile(pattern=glob2rx(paste0("*.shp")))
}

print("aoi ok")
dem <- setMinMax(raster("DEM.tif"))
print("dem ok")

####### Read Sentinel-2 L2A images

elenco_file_zip <- list.files(pattern=glob2rx('S2*.zip'))
num_zip <- length(elenco_file_zip)

for (i in elenco_file_zip[1:num_zip]){
  s2_zip = i
  print(s2_zip)
  unzip(s2_zip, files = NULL, list = FALSE, overwrite = TRUE,
        junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)
}

elenco_file_SAFE <- list.files(pattern=glob2rx('S2*_MSIL2A*.SAFE'))

num_SAFE=length(elenco_file_SAFE)
print(num_SAFE)


x <- elenco_file_SAFE[1]

sensore <- substr(x,1,3)
tile <- substr(x,39,44)
ripresa <- substr(x,34,37)
date <- substr(x,12,19)
anno <- substr(x,12,15)
ora_acq <- substr(x,21,26)
ora_proc <- substr(x,55,60)
inter_name_file <- substr(x,12,27)
n_code <- substr(x,28,32)
end_name_file<-substr(x,12,60)

# Find Fmask
Fmask_20m <- raster(list.files(pattern=glob2rx('*Fmask*.tif')))
print("raster fmask ok")
if (aoifile == TRUE) {
  Fmask_20m <- crop(Fmask_20m, aoi)
}

path_bands_20m <- "IMG_DATA/R20m"
image_l2_name <- paste0(sensore,"_MSIL2A_",end_name_file,".SAFE",collapse = NULL, recycle0 = FALSE)
path_bands_1 <- dir(path=paste0(repo_data,"/",image_l2_name,"/GRANULE"),full.names = FALSE,recursive = FALSE,pattern = "*")
path_bands_def_20m <- paste0(repo_data,"/",image_l2_name,"/GRANULE/",path_bands_1,"/",path_bands_20m)
print(path_bands_def_20m)
print(getwd())
setwd(path_bands_def_20m)
print(getwd())
print(list.files(pattern=glob2rx("*")))

### Read INPUT BANDS
B3_20m <- raster(list.files(pattern=glob2rx('*B03*.jp2')))
print("raster bands ok b3")
B4_20m <- raster(list.files(pattern=glob2rx('*B04*.jp2')))
print("raster bands ok b4")
B11_20m <- raster(list.files(pattern=glob2rx("*B11*.jp2")))
print("raster bands ok b11")
print(list.files(pattern=glob2rx("*")))
SCL_20m <- raster(list.files(pattern=glob2rx("*SCL*.jp2")))
print("raster scl20 ok")

### Crop on AOI if presents

if (aoifile == TRUE){
  B3_20m <- crop(B3_20m, aoi)
  B4_20m <- crop(B4_20m, aoi)
  B11_20m <- crop(B11_20m, aoi)
  SCL_20m <- crop(SCL_20m, aoi)
  
  print("data cropped")
}

##################### Pre-processing DEM

# check CRS
if (compareCRS(dem, B3_20m) == FALSE){
  dem <- projectRaster(dem, crs = crs(B3_20m), res=20)
}
# check resolution
rres <- xres(dem) / xres(B3_20m)
if (rres != 1){
  dem <- resample(dem, B3_20m, method="bilinear")
}
# check extent
if (extent(B3_20m) != extent(dem)){
  dem <- crop(dem, B3_20m)
}
if (extent(dem) != extent(B3_20m)){
  dem <- setExtent(dem, extent(B3_20m), keepres=FALSE)
}

print("dem ready")

###### No data masking
partial <- FALSE

total_NA_cells <- sum(SCL_20m [!is.na(SCL_20m)] == 0)
total_raster_cells <- ncell(SCL_20m)
total_NA_fraction <- round(total_NA_cells*100/total_raster_cells,0)

if (total_NA_fraction >= 2){
  
  partial <- TRUE
  
  nd_20m <- calc(SCL_20m, fun=function(x){x == 0})
  nodata_20m <- clump(nd_20m,directions=4, gaps=TRUE)
  mask_crop_nodata_20m <- (nodata_20m != 1)
  
  B3_20m <- mask(B3_20m, mask_crop_nodata_20m, inverse=TRUE, maskvalue=1, updatevalue=NA)
  B4_20m <- mask(B4_20m, mask_crop_nodata_20m, inverse=TRUE, maskvalue=1, updatevalue=NA)
  B11_20m <- mask(B11_20m, mask_crop_nodata_20m, inverse=TRUE, maskvalue=1, updatevalue=NA)
  
  dem_masking <- mask(dem,mask_crop_nodata_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
  print("no data masked")
  
}

if (partial == FALSE) {
  dem_masking <- dem
}

###### Fmask No data masking

if (partial == TRUE){
  maskFmask_crop_compl_20m <- (Fmask_20m != 255)
  Fmask_20m <- mask(Fmask_20m,maskFmask_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
}

############################################### Set output wd

setwd(repo_data)
output_folder <- as.character("output")
dir.create(output_folder,showWarnings = FALSE)
setwd(output_folder)

###############################################

B3_20m_FLOAT = B3_20m/10000
B4_20m_FLOAT = B4_20m/10000
B11_20m_FLOAT = B11_20m/10000

print("Snow cover detection starting")

##################################################### Snow Cover Processor - START

# PARAMETERS

rB4_darkcloud <- 0.300
n1 <- 4000
n2 <- 1500
r1 <- 0.200
r2 <- 0.040
dz <- 100
fsnow_lim <- 0.1
fclear_lim <- 0.100
fsnow_total_lim <- 0.001
rB4_backtocloud <- 0.100
s1 <- 0.100
s2 <- 0.250

### NDSI at 20m

NDSI_20m <- ((B3_20m_FLOAT - B11_20m_FLOAT)/(B3_20m_FLOAT + B11_20m_FLOAT))*10000
NDSI_20m[NDSI_20m < -10000] = -10000
NDSI_20m[NDSI_20m > 10000] = 10000

### CLOUD PASS 1
cloud_pass0 <- Fmask_20m

cloud_pass0[cloud_pass0 == 1] <- 0
cloud_pass0[cloud_pass0 == 4] <- 1
cloud_pass0[cloud_pass0 == 2] <- 1
cloud_pass0[cloud_pass0 != 1] <- 0

if (partial == TRUE){
 cloud_pass0 <- mask(cloud_pass0,mask_crop_nodata_20m,inverse=TRUE,maskvalue=1,updatevalue=NA) 
}


# B4 noise reduction
B4_20m_FLOAT_filtered <- focal(B4_20m_FLOAT, w=matrix(1/9, nc=3, nr=3), na.rm=TRUE, pad=TRUE, padValue=NA)

# Dark-Cloud pixel
dark_clouds <- overlay(Fmask_20m,
                       B4_20m_FLOAT_filtered,
                       fun = function(x,y) {return(x == 4 & y < rB4_darkcloud)})

cloud_pass1 <- cloud_pass0 - dark_clouds
cloud_pass1[is.na(cloud_pass1)] <- 0

if (partial == TRUE){
 cloud_pass1 <- mask(cloud_pass1,mask_crop_nodata_20m,inverse=TRUE,maskvalue=1,updatevalue=NA) 
}

#### SNOW PASS 1

snow_pass1 <- overlay(NDSI_20m,
                      B4_20m_FLOAT,
                      B11_20m_FLOAT,
                      fun=function(a,b,c,d) {return(a>n1 & b>r1 & c<s1)})

### TEST "ENOUGH SNOW?"

# calculates the fraction of total snow cover in the image
total_snow_cells <- sum(snow_pass1 [!is.na(snow_pass1)]  == 1)
total_cells <- ncell(snow_pass1[!is.na(snow_pass1)])
total_snow_fraction <- total_snow_cells*100/total_cells


if (total_snow_fraction > fsnow_total_lim){
  
  ### SNOWLINE ELEVATION

bh <- 100
zmax <- round(maxValue(dem_masking),-2)
zmin <- zmax - bh

#calcola estensione banda altitudinale libera da nuvole
dem_band <- calc(dem_masking, fun = function(x) {return(x > zmin & x < zmax)})
dem_band_pixel <- sum(dem_band[!is.na(dem_band)])

cloud_band <- overlay(cloud_pass0,
                      dem_band,
                      fun=function(x,y) {return(x == 1 & y == 1)})

cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel

while (dem_band_free_pixel < fclear_lim & zmin>0) {
  
  zmax <- zmax - bh
  zmin <- zmin - bh
  
  dem_band <- calc(dem_masking, fun = function(x) {return(x > zmin & x < zmax)})
  dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
  
  cloud_band <- overlay(cloud_pass0,
                        dem_band,
                        fun=function(x,y) {return(x == 1 & y ==1)})
  
  cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
  dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
  
}

#calcola estensione neve x banda altitudinale
snow_band <- overlay(snow_pass1,
                     dem_band,
                     fun=function(x,y) {return(x == 1 & y == 1)})

snow_band_pixel <- sum(snow_band[!is.na(snow_band)]==1)

#calcola la frazione neve x banda altitudinale
band_snow_fraction <- snow_band_pixel/dem_band_pixel*100

#trova la banda altitudinale più bassa per la quale la snow cover fraction è > 0.1
while (band_snow_fraction > fsnow_lim & snow_band_pixel > 0){
  
  zmax <- zmax - bh
  zmin <- zmin - bh
  
  dem_band <- calc(dem, fun = function(x) {return(x > zmin & x < zmax)})
  
  cloud_band <- overlay(cloud_pass0,
                        dem_band,
                        fun=function(x,y) {return(x == 1 & y == 1)})
  
  dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
  cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
  dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
  
  while (dem_band_free_pixel < fclear_lim & zmin>0) {
    
    zmax <- zmax - bh
    zmin <- zmin - bh
    
    dem_band <- calc(dem_masking, fun = function(x) {return(x > zmin & x < zmax)})
    dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
    
    cloud_band <- overlay(cloud_pass0,
                          dem_masking,
                          fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
    
    cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
    dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
    
  }
  
  snow_band <- overlay(snow_pass1,
                       dem_band,
                       fun=function(x,y) {return(x == 1 & y == 1)})
  
  snow_band_pixel <- sum(snow_band[!is.na(snow_band)]==1)
  band_snow_fraction <- snow_band_pixel/dem_band_pixel*100
  
}

#definisci limite inferiore neve come due bande inferiori a quella precedentemente trovata
zs <- zmin - (2*bh)

if (zs < round(minValue(dem_masking),-2)){
  zs <- minValue(dem_masking)
}
  
  ### SNOW PASS 2
  snow_pass2 <- overlay(NDSI_20m,
                        B4_20m_FLOAT,
                        dem_masking,
                        B11_20m_FLOAT,
                        fun=function(a,b,c,d) {return(a>n2 & b>r2 & c>zs & d<s2)})
  
  ### CLOUD PASS2
  cloud_pass2 <- overlay(snow_pass2,
                         dark_clouds,
                         B4_20m_FLOAT,
                         fun=function(x,y,z) {return(x != 1 & y == 1 & z > rB4_backtocloud)})
  
  cloud_pass2[is.na(cloud_pass2)] <- 0
  
  if (partial == TRUE){
    cloud_pass2 <- mask(cloud_pass2,mask_crop_nodata_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  }
  
  
  ###CLOUD FINAL
  cloud_final <- overlay(cloud_pass1,
                         cloud_pass2,
                         fun = function(x, y) {x == 1 | y == 1})
  
  
  ############################# FINAL OUTPUT CREATION
  
  # CREATE RASTER FROM SNOW PASS2 AND CLOUD FINAL -> WITH 4 values: 100 SNOW; 0 NO-SNOW; 205 CLOUD; 254 NO DATA
  snow_pass2[snow_pass2 == 1] <- 100
  cloud_final[cloud_final == 1] <- 205
  
  output_SnowCloud_20m <- sum(snow_pass2,cloud_final)
  output_SnowCloud_20m[output_SnowCloud_20m == 305] <- 100
  output_SnowCloud_20m[is.na(output_SnowCloud_20m)] <- 254
  
  ### Remove small patches
  land <- calc(output_SnowCloud_20m, fun=function(x){x == 0})
  patches <- clump(land,directions=4, gaps=TRUE)
  # get frequency table    
  f<-freq(patches)
  # save frequency table as data frame
  f<-as.data.frame(f)
  # which rows of the data.frame are only represented by clumps under 4 pixels?
  str(which(f$count <= 5))
  # which values do these correspond to?
  str(f$value[which(f$count <= 5)])
  # put these into a vector of clump ID's to be removed
  excludeID <- f$value[which(f$count <= 5)]
  # assign NA to all clumps whose IDs are found in excludeID
  patches[patches %in% excludeID] <- 0
  # reclassify: 0 are pixels to remove and all the rest is 1
  patches[patches != 0] <- 1
  patches[is.na(patches)] <- 1
  # produce mask
  patches_mask <- (patches != 0)
  output_SnowCloud_20m_gap <- mask(output_SnowCloud_20m,patches_mask,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
  # replace only NA values
  output_SnowCloud_20m_gap_filled <- focal(output_SnowCloud_20m_gap,
                                           w=matrix(1,nrow=7,ncol=7),
                                           fun=max,
                                           na.rm=TRUE,
                                           NAonly=TRUE,
                                           pad=TRUE)
  
  output_SnowCloud_20m_name <- paste0("SnowCloud_mask_20m.tif")
  writeRaster(output_SnowCloud_20m_gap_filled, output_SnowCloud_20m_name, datatype='INT2S', format = "GTiff", overwrite=TRUE)
  
} else {
  
  ############################# FINAL OUTPUT CREATION
  
  cloud_pass3 <- overlay(snow_pass1,
                         dark_clouds,
                         B4_20m_FLOAT,
                         fun = function(x,y,z) {return(x != 1 & y == 1 & z > rB4_backtocloud)})
  
  cloud_pass3[is.na(cloud_pass3)] <- 0
  
  if (partial == TRUE){
    cloud_pass3 <- mask(cloud_pass3,mask_crop_nodata_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  }
  
  
  cloud_final <- overlay(cloud_pass1,
                         cloud_pass3,
                         fun = function(x, y) {return(x == 1 | y == 1)})
  
  
  # CREATE RASTER FROM SNOW PASS2 AND CLOUD FINAL -> WITH 4 values: 100 SNOW; 0 NO-SNOW; 205 CLOUD; 254 NO DATA
  snow_pass1[snow_pass1 == 1] <- 100
  cloud_final[cloud_final == 1] <- 205
  
  output_SnowCloud_20m <- sum(snow_pass1,cloud_final)
  output_SnowCloud_20m[output_SnowCloud_20m == 305] <- 100
  output_SnowCloud_20m[is.na(output_SnowCloud_20m)] <- 254
  
  ### Remove small patches
  land <- calc(output_SnowCloud_20m, fun=function(x){x == 0})
  patches <- clump(land,directions=4, gaps=TRUE)
  # get frequency table    
  f<-freq(patches)
  # save frequency table as data frame
  f<-as.data.frame(f)
  # which rows of the data.frame are only represented by clumps under 4 pixels?
  str(which(f$count <= 5))
  # which values do these correspond to?
  str(f$value[which(f$count <= 5)])
  # put these into a vector of clump ID's to be removed
  excludeID <- f$value[which(f$count <= 5)]
  # assign NA to all clumps whose IDs are found in excludeID
  patches[patches %in% excludeID] <- 0
  # reclassify: 0 are pixels to remove and all the rest is 1
  patches[patches != 0] <- 1
  patches[is.na(patches)] <- 1
  # produce mask
  patches_mask <- (patches != 0)
  output_SnowCloud_20m_gap <- mask(output_SnowCloud_20m,patches_mask,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
  # replace only NA values
  output_SnowCloud_20m_gap_filled <- focal(output_SnowCloud_20m_gap,
                                           w=matrix(1,nrow=7,ncol=7),
                                           fun=max,
                                           na.rm=TRUE,
                                           NAonly=TRUE,
                                           pad=TRUE)
  
  output_SnowCloud_20m_name <- paste0("SnowCloud_mask_20m.tif")
  writeRaster(output_SnowCloud_20m_gap_filled, output_SnowCloud_20m_name, datatype='INT2S', format = "GTiff", overwrite=TRUE)
  
}

############################################################## Snow Cover Processor - END

print("Processing finished")

############ remove temp files

dir_base<-getwd()

setwd(dirname(rasterTmpFile()))
file.remove(list.files(pattern=glob2rx('*.*')))

setwd(dir_base)

############

setwd(repo_data)

print("Workflow finished")

