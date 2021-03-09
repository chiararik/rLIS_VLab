# Snow Cover Processor
 Snow Cover Processor (SCP) repository
 
 ## Input file  
 As input use a zip file named after "input.zip" with:
 - Sentinel-2 L2A Bottom-of-atmosphere with the original name unzipped,  
 - Sentinel-2 L1C pre-processed with FMASK with the original name unzipped,  
 - a DEM, named after "DEM.tif",  
 - a shapefile of the area of interest, named after "aoi.shp", with the same CRS of the Sentinel-2 (optional)  
  
The DEM will be pre-processed in the workflow to meet the following conditions: the same CSR of the S2 image, the same extent and resolution as the SWIR band (20 m), therefore, to decrease the processing time it is useful to provide a DEM already projected and with the required resolution, in addition to providing an area of interest.
