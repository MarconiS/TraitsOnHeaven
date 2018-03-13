def h5refl2array(refl_filename):
    #Read in reflectance hdf5 file (include full or relative path if data is located in a different directory)
    hdf5_file = h5py.File(refl_filename,'r')

    #Get the site name
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[1]
    
    #Extract the reflectance & wavelength datasets
    refl = hdf5_file[sitename]['Reflectance']
    reflArray = refl['Reflectance_Data']
    refl_shape = reflArray.shape
    wavelengths = refl['Metadata']['Spectral_Data']['Wavelength']
    
    #Create dictionary containing relevant metadata information
    metadata = {}
    metadata['shape'] = reflArray.shape
    metadata['mapInfo'] = refl['Metadata']['Coordinate_System']['Map_Info'].value

    #Extract no data value & set no data value to NaN
    metadata['noDataVal'] = float(reflArray.attrs['Data_Ignore_Value'])
    metadata['scaleFactor'] = float(reflArray.attrs['Scale_Factor'])
    
    #Extract bad band windows
    metadata['bad_band_window1'] = (refl.attrs['Band_Window_1_Nanometers'])
    metadata['bad_band_window2'] = (refl.attrs['Band_Window_2_Nanometers'])
    
    #Extract projection information
    metadata['projection'] = refl['Metadata']['Coordinate_System']['Proj4'].value
    metadata['epsg'] = int(refl['Metadata']['Coordinate_System']['EPSG Code'].value)
    
    #Extract map information: spatial extent & resolution (pixel size)
    mapInfo = refl['Metadata']['Coordinate_System']['Map_Info'].value
    mapInfo_string = str(mapInfo); 
    mapInfo_split = mapInfo_string.split(",")
    
    #Extract the resolution & convert to floating decimal number
    metadata['res'] = {}
    metadata['res']['pixelWidth'] = float(mapInfo_split[5])
    metadata['res']['pixelHeight'] = float(mapInfo_split[6])
    
    #Extract the upper left-hand corner coordinates from mapInfo
    xMin = float(mapInfo_split[3]) #convert from string to floating point number
    yMax = float(mapInfo_split[4])
    #Calculate the xMax and yMin values from the dimensions
    xMax = xMin + (refl_shape[1]*metadata['res']['pixelWidth']) #xMax = left edge + (# of columns * resolution)",
    yMin = yMax - (refl_shape[0]*metadata['res']['pixelHeight']) #yMin = top edge - (# of rows * resolution)",
    metadata['extent'] = (xMin,xMax,yMin,yMax) #useful format for plotting
    metadata['ext_dict'] = {}
    metadata['ext_dict']['xMin'] = xMin
    metadata['ext_dict']['xMax'] = xMax
    metadata['ext_dict']['yMin'] = yMin
    metadata['ext_dict']['yMax'] = yMax
    hdf5_file.close        
    
    return reflArray, metadata, wavelengths

def calc_clip_index(clipExtent, h5Extent, xscale=1, yscale=1):
    
    h5rows = h5Extent['yMax'] - h5Extent['yMin']
    h5cols = h5Extent['xMax'] - h5Extent['xMin']    
    
    ind_ext = {}
    ind_ext['xMin'] = round((clipExtent['xMin']-h5Extent['xMin'])/xscale)
    ind_ext['xMax'] = round((clipExtent['xMax']-h5Extent['xMin'])/xscale)
    ind_ext['yMax'] = round(h5rows - (clipExtent['yMin']-h5Extent['yMin'])/xscale)
    ind_ext['yMin'] = round(h5rows - (clipExtent['yMax']-h5Extent['yMin'])/yscale)
    
    return ind_ext



def stack_subset_bands(reflArray,reflArray_metadata,bands,clipIndex):
    
    subArray_rows = clipIndex['yMax'] - clipIndex['yMin']
    subArray_cols = clipIndex['xMax'] - clipIndex['xMin']
    
    stackedArray = np.zeros((subArray_rows,subArray_cols,len(bands)),dtype = np.int16)
    band_clean_dict = {}
    band_clean_names = []
    
    for i in range(len(bands)):
        band_clean_names.append("b"+str(bands[i])+"_refl_clean")
        band_clean_dict[band_clean_names[i]] = subset_clean_band(reflArray,reflArray_metadata,clipIndex,bands[i])
        stackedArray[...,i] = band_clean_dict[band_clean_names[i]]
                        
    return stackedArray


def subset_clean_band(reflArray,reflArray_metadata,clipIndex,bandIndex):
    
    bandCleaned = reflArray[clipIndex['yMin']:clipIndex['yMax'],clipIndex['xMin']:clipIndex['xMax'],bandIndex].astype(np.int16)
    
    return bandCleaned 




def array2raster(newRaster,reflBandArray,reflArray_metadata, extent, ras_dir): 
    
    NP2GDAL_CONVERSION = {
      "uint8": 1,
      "int8": 1,
      "uint16": 2,
      "int16": 3,
      "uint32": 4,
      "int32": 5,
      "float32": 6,
      "float64": 7,
      "complex64": 10,
      "complex128": 11,
    }
    
    pwd = os.getcwd()
    os.chdir(ras_dir) 
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = reflBandArray.shape[2]
    pixelWidth = float(reflArray_metadata['res']['pixelWidth'])
    pixelHeight = -float(reflArray_metadata['res']['pixelHeight'])
    originX = extent['xMin']
    originY = extent['yMax']
    
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    #outband = outRaster.GetRasterBand(1)
    #outband.WriteArray(reflBandArray[:,:,x])
    for band in range(bands):
        outRaster.GetRasterBand(band+1).WriteArray(reflBandArray[:,:,band])
    
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg']) 
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd) 

def stripe2Raser(f, pt):
    
    full_path = pt+f
    refl, refl_md, wavelengths = h5refl2array(full_path)
    refl_md['extent']
    rgb = np.r_[0:425]
    rgb = np.delete(rgb, np.r_[419:426])
    rgb = np.delete(rgb, np.r_[283:315])
    rgb = np.delete(rgb, np.r_[192:210])    
    xmin, xmax, ymin, ymax = refl_md['extent']
    clipExtent = {}
    clipExtent['xMin'] = int(xmin) 
    clipExtent['yMin'] = int(ymin) 
    clipExtent['yMax'] = int(ymax) 
    clipExtent['xMax'] = int(xmax) 

    subInd = calc_clip_index(clipExtent,refl_md['ext_dict']) 
    subInd['xMax'] = int(subInd['xMax'])
    subInd['xMin'] = int(subInd['xMin'])
    subInd['yMax'] = int(subInd['yMax'])
    subInd['yMin'] = int(subInd['yMin'])
    reflBandArray = stack_subset_bands(refl,refl_md,rgb,subInd)

    subArray_rows = subInd['yMax'] - subInd['yMin']
    subArray_cols = subInd['xMax'] - subInd['xMin']
    hcp = np.zeros((subArray_rows,subArray_cols,len(rgb)),dtype = np.int16)

    band_clean_dict = {}
    band_clean_names = []
    for i in range(len(rgb)):
        band_clean_names.append("b"+str(rgb[i])+"_refl_clean")
        band_clean_dict[band_clean_names[i]] = refl[:,:,rgb[i]].astype(np.int16)
        hcp[...,i] = band_clean_dict[band_clean_names[i]]

    reflArray_metadata = refl_md
    newRaster = f.replace(' ','')[:-3].upper()
    newRaster = str(clipExtent['xMin'])+str((clipExtent['yMax']))+'.tif'

    ras_dir = '/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/rasters/'
    array2raster(newRaster,hcp,reflArray_metadata, clipExtent, ras_dir)


import numpy as np
import h5py
import gdal, osr
import sys
import ogr, os

#pt = "/ufrc/ewhite/s.marconi/Marconi2018/spatialPhase/OSBS/Reflectance/"
#f = "NEON_D03_OSBS_DP3_394000_3280000_reflectance.h5"
#stripe2Raser(f,  pt)

stripe2Raser(sys.argv[1],  sys.argv[2])
print("at least raster is good")

