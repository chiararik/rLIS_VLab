suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(gdalUtils))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(sf))

print("Starting Workflow")

##### Input files:

# Sentinel-2 L2A with FMASK output (mandatory)
# DEM (mandatory)
# AOI (optional)

###### Set working directory

repo <- getwd()
print(repo)
repo_data <- paste0(repo,"/data")
setwd(repo_data)

###### Read AOI (if present) and DEM

aoizip <- file.exists("aoi.zip")
aoifile <- FALSE

if (aoizip == TRUE){
  aoifile <- TRUE
  unzip("aoi.zip", files = NULL, list = FALSE, overwrite = TRUE,
        junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)
  aoi <- shapefile("aoi.shp")
}

dem <- raster("DEM.tif")

####### Read Sentinel-2 L2A images

elenco_file_zip <- list.files(pattern=glob2rx('S2*_MSIL2A*.zip'))
num_zip <- length(elenco_file_zip)

for (i in elenco_file_zip){
  s2_zip = i
  unzip(s2_zip, files = NULL, list = FALSE, overwrite = TRUE,
        junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)
}

elenco_file_SAFE <- list.files(pattern=glob2rx('S2*_MSIL2A*.SAFE'))

num_SAFE=length(elenco_file_SAFE)

n <- 1

for (i in elenco_file_SAFE[1:num_SAFE]){
  
  x=i
  
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
  image_l2_name <- paste0(sensore,"_MSIL2A_",end_name_file,".SAFE",collapse = NULL, recycle0 = FALSE)
  path_bands_1 <- dir(path=paste0(repo_data,"/",image_l2_name,"/GRANULE"),full.names = FALSE,recursive = FALSE,pattern = "L2A_*")
  path_mask_cloud <- paste0(repo_data,"/",image_l2_name,"/GRANULE/",path_bands_1,"/FMASK_DATA")
  print(path_mask_cloud)
  setwd(path_mask_cloud)
  
  Fmask_20m <- raster(list.files(pattern=glob2rx('*_Fmask.tif')))
  if (aoifile == TRUE) {
    Fmask_20m <- crop(Fmask_20m, aoi)
  }
  
  setwd(repo_data)
  
  path_bands_20m <- "IMG_DATA/R20m"
  path_bands_def_20m <- paste0(repo_data,"/",image_l2_name,"/GRANULE/",path_bands_1,"/",path_bands_20m)
  setwd(path_bands_def_20m)
  print(path_bands_def_20m)
  print(getwd())
  
  ### Read INPUT BANDS
  B3_20m <- raster(list.files(pattern=glob2rx('*B03_20m.jp2')))
  B4_20m <- raster(list.files(pattern=glob2rx('*B04_20m.jp2')))
  B11_20m <- raster(list.files(pattern=glob2rx("*B11_20m.jp2")))
  #gdal_translate("*SCL_20m.jp2","L2A_SCL_20m.tif")
  gdal_translate(list.files(pattern=glob2rx("*SCL_20m.jp2"))[1],"L2A_SCL_20m.tif")
  SCL_20m <- raster(list.files(pattern=glob2rx('L2A_SCL_20m.tif')))
  
  ### Crop on AOI if presents
  
  if (aoifile == TRUE){
    B3_20m <- crop(B3_20m, aoi)
    B4_20m <- crop(B4_20m, aoi)
    B11_20m <- crop(B11_20m, aoi)
    SCL_20m <- crop(SCL_20m, aoi)
  }
  
  ###### Mask no data from SCL 
  
  maskSCL_crop_compl_20m <- (SCL_20m != 0)
  
  B3_crop_20m_masking_SCL <- mask(B3_20m, maskSCL_crop_compl_20m, inverse=TRUE, maskvalue=1, updatevalue=NA)
  B4_crop_20m_masking_SCL <- mask(B4_20m, maskSCL_crop_compl_20m, inverse=TRUE, maskvalue=1, updatevalue=NA)
  B11_crop_20m_masking_SCL <- mask(B11_20m, maskSCL_crop_compl_20m, inverse=TRUE, maskvalue=1, updatevalue=NA)
  
  ###### Fmask no data masking
  
  maskFmask_crop_compl_20m <- (Fmask_20m != 255)
  Fmask_crop_20m_masking_SCL <- mask(Fmask_20m,maskFmask_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
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
  
  ###### No data masking from DEM
  
  dem_masking_SCL <- mask(dem,maskSCL_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
  ############################################### Set output wd 

  setwd(repo_data)
  output_folder <- as.character("output")
  dir.create(output_folder,showWarnings = FALSE)
  setwd(output_folder)
  
  ###############################################
  
  B3_crop_20m_masking_SCL_FLOAT = B3_crop_20m_masking_SCL/10000
  B4_crop_20m_masking_SCL_FLOAT = B4_crop_20m_masking_SCL/10000
  B11_crop_20m_masking_SCL_FLOAT = B11_crop_20m_masking_SCL/10000
  
  message3 <- paste0("Input data of image ",n,"/",num_SAFE," ready, snow cover detection starting")
  print(message3)
  
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
  
  NDSI_crop_20m <- ((B3_crop_20m_masking_SCL_FLOAT - B11_crop_20m_masking_SCL_FLOAT)/(B3_crop_20m_masking_SCL_FLOAT + B11_crop_20m_masking_SCL_FLOAT))*10000
  NDSI_crop_20m[NDSI_crop_20m < -10000] = -10000
  NDSI_crop_20m[NDSI_crop_20m > 10000] = 10000
  
  ### CLOUD PASS 1
  cloud_pass0 <- Fmask_crop_20m_masking_SCL
  
  cloud_pass0[cloud_pass0 == 1] <- 0
  cloud_pass0[cloud_pass0 == 4] <- 1
  cloud_pass0[cloud_pass0 == 2] <- 1
  cloud_pass0[cloud_pass0 != 1] <- 0
  
  cloud_pass0_mask <- mask(cloud_pass0,maskSCL_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
  # B4 noise reduction
  B4_crop_20m_masking_SCL_FLOAT_filtered <- focal(B4_crop_20m_masking_SCL_FLOAT, w=matrix(1/9, nc=3, nr=3), fun=median, na.rm=TRUE, pad=TRUE, padValue=NA)
  
  # Dark-Cloud pixel
  dark_clouds <- overlay(Fmask_crop_20m_masking_SCL,
                         B4_crop_20m_masking_SCL_FLOAT_filtered,
                         fun = function(x,y) {return(x == 4 & y < rB4_darkcloud)})
  
  cloud_pass1 <- cloud_pass0_mask - dark_clouds
  cloud_pass1[is.na(cloud_pass1)] <- 0
  cloud_pass1_mask <- mask(cloud_pass1,maskSCL_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
  
  #### SNOW PASS 1
  
  snow_pass1 <- overlay(NDSI_crop_20m,
                        B4_crop_20m_masking_SCL_FLOAT,
                        B11_crop_20m_masking_SCL_FLOAT,
                        cloud_pass1_mask,
                        fun=function(a,b,c,d) {return(a>n1 & b>r1 & c<s1 & d != 1)})
  
  ### TEST "ENOUGH SNOW?"
  
  # calculates the fraction of total snow cover in the image
  total_snow_cells <- sum(snow_pass1 [!is.na(snow_pass1)]  == 1)
  total_cells <- ncell(snow_pass1[!is.na(snow_pass1)])
  total_snow_fraction <- total_snow_cells*100/total_cells
  
  
  if (total_snow_fraction > fsnow_total_lim){
    
    ### SNOWLINE ELEVATION
    
    bh <- 100
    zmax <- round(maxValue(dem_masking_SCL),-2)
    zmin <- zmax - bh
    
    # calculates cloud-free altitude band extension
    dem_band <- calc(dem_masking_SCL, fun = function(x) {return(x > zmin & x < zmax)})
    dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
    
    cloud_band <- overlay(cloud_pass0_mask,
                          dem_masking_SCL,
                          fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
    
    cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
    dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
    
    while (cloudy_band_pixel == dem_band_pixel & dem_band_pixel > 0) {
      
      zmax <- zmax - bh
      zmin <- zmin - bh
      
      dem_band <- calc(dem_masking_SCL, fun = function(x) {return(x > zmin & x < zmax)})
      dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
      
      cloud_band <- overlay(cloud_pass0_mask,
                            dem_masking_SCL,
                            fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
      
      cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
      dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
      
    }
    
    # calculate snow extension x altitude band
    snow_band <- overlay(snow_pass1,
                         dem_masking_SCL,
                         fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
    
    snow_band_pixel <- sum(snow_band[!is.na(snow_band)]==1)
    
    # calculates the snow fraction x altitude band
    band_snow_fraction <- snow_band_pixel/dem_band_free_pixel*100
    
    # find the lowest altitudinal band for which the snow cover fraction is > 0.1
    while (band_snow_fraction > fsnow_lim & snow_band_pixel > 0){
      
      zmax <- zmax - bh
      zmin <- zmin - bh
      
      snow_band <- overlay(snow_pass1,
                           dem_masking_SCL,
                           fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
      
      dem_band <- calc(dem, fun = function(x) {return(x > zmin & x < zmax)})
      
      cloud_band <- overlay(cloud_pass0_mask,
                            dem_masking_SCL,
                            fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
      
      dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
      cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
      dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
      
      while (cloudy_band_pixel == dem_band_pixel & dem_band_pixel > 0) {
        
        zmax <- zmax - bh
        zmin <- zmin - bh
        
        dem_band <- calc(dem_masking_SCL, fun = function(x) {return(x > zmin & x < zmax)})
        dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
        
        cloud_band <- overlay(cloud_pass0_mask,
                              dem_masking_SCL,
                              fun=function(x,y) {return(x == 1 & (y >= zmin & y < zmax))})
        
        cloudy_band_pixel <- sum(cloud_band[!is.na(cloud_band)]==1)
        dem_band_free_pixel <- dem_band_pixel - cloudy_band_pixel
        
      }
      
      snow_band_pixel <- sum(snow_band[!is.na(snow_band)]==1)
      band_snow_fraction <- snow_band_pixel/dem_band_free_pixel*100
      
    }
    
    # define lower snow limit as two bands lower than the one previously found
    zs <- zmin - (2*bh)
    
    ### SNOW PASS 2    
    snow_pass2 <- overlay(NDSI_crop_20m,
                          B4_crop_20m_masking_SCL_FLOAT,
                          dem_masking_SCL,
                          B11_crop_20m_masking_SCL_FLOAT,
                          fun=function(a,b,c,d) {return(a>n2 & b>r2 & c>zs & d<s2)})
    
    ### CLOUD PASS2
    cloud_pass2 <- overlay(snow_pass2,
                           dark_clouds,
                           B4_crop_20m_masking_SCL_FLOAT,
                           fun=function(x,y,z) {return(x != 1 & y == 1 & z > rB4_backtocloud)})
    
    cloud_pass2[is.na(cloud_pass2)] <- 0
    cloud_pass2_mask <- mask(cloud_pass2,maskSCL_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
    
    ###CLOUD FINAL
    cloud_final <- overlay(cloud_pass1_mask,
                           cloud_pass2_mask, 
                           fun = function(x, y) {x == 1 | y == 1})
    
    
    ############################# FINAL OUTPUT CREATION
    
    # CREATE RASTER FROM SNOW PASS2 AND CLOUD FINAL -> WITH 4 values: 100 SNOW; 0 NO-SNOW; 205 CLOUD; 255 NO DATA
    snow_pass2[snow_pass2 == 1] <- 100
    cloud_final[cloud_final == 1] <- 205
    
    output_SnowCloud_20m <- sum(snow_pass2,cloud_final)
    output_SnowCloud_20m[output_SnowCloud_20m == 305] <- 100
    output_SnowCloud_20m[is.na(output_SnowCloud_20m)] <- 254
    
    output_SnowCloud_20m_name <- paste0("Snow&Cloud_mask_20m_",date,"_",sensore,"_L2A_",tile,"_",ripresa,"_TA",ora_acq,"_TP",ora_proc,".tif")
    writeRaster(output_SnowCloud_20m, output_SnowCloud_20m_name, datatype='INT2S', format = "GTiff", overwrite=TRUE)
    
  } else { 
    
    ############################# FINAL OUTPUT CREATION
    
    cloud_pass3 <- overlay(snow_pass1,
                           dark_clouds,
                           B4_crop_20m_masking_SCL_FLOAT,
                           fun = function(x,y,z) {return(x != 1 & y == 1 & z > rB4_backtocloud)})
    
    cloud_pass3[is.na(cloud_pass3)] <- 0
    cloud_pass3_mask <- mask(cloud_pass3,maskSCL_crop_compl_20m,inverse=TRUE,maskvalue=1,updatevalue=NA)
    
    cloud_final <- overlay(cloud_pass1_mask,
                           cloud_pass3_mask, 
                           fun = function(x, y) {return(x == 1 | y == 1)})
    
    
    # CREATE RASTER FROM SNOW PASS2 AND CLOUD FINAL -> WITH 4 values: 100 SNOW; 0 NO-SNOW; 205 CLOUD; 255 NO DATA
    snow_pass1[snow_pass1 == 1] <- 100
    cloud_final[cloud_final == 1] <- 205
    
    output_SnowCloud_20m <- sum(snow_pass1,cloud_final)
    output_SnowCloud_20m[output_SnowCloud_20m == 305] <- 100
    output_SnowCloud_20m[is.na(output_SnowCloud_20m)] <- 254
    
    output_SnowCloud_20m_name <- paste0("Snow&Cloud_mask_20m_",date,"_",sensore,"_L2A_",tile,"_",ripresa,"_TA",ora_acq,"_TP",ora_proc,".tif")
    writeRaster(output_SnowCloud_20m, output_SnowCloud_20m_name, datatype='INT2S', format = "GTiff", overwrite=TRUE)
    
  }
  
  ############################################################## Snow Cover Processor - END
  
  message4 <- paste0("Processing of image ",n,"/",num_SAFE," finished")
  print(message4)
  
  n <- n+1
  
  ############ remove temp files

  dir_base<-getwd()
  
  setwd(dirname(rasterTmpFile()))
  file.remove(list.files(pattern=glob2rx('*.*')))
  
  setwd(dir_base)
  
  ############
  
  setwd(repo_data)
  
}

print("Workflow Finished")