def h5refl2array(refl_filename):


    hdf5_file = h5py.File(refl_filename,'r')
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[1]
    reflArray = hdf5_file['Reflectance']
    reflArray = np.swapaxes(reflArray, 0,2)
    reflArray = np.swapaxes(reflArray, 0,1)
    reflArray.shape
    refl_shape = reflArray.shape
    wavelengths = hdf5_file['wavelength']
    #Create dictionary containing relevant metadata information
    metadata = {}
    metadata['shape'] = reflArray.shape
    metadata['mapInfo'] = hdf5_file['map info'].value
    #Extract no data value & set no data value to NaN
    metadata['noDataVal'] = float(max((reflArray[1,1,:])))
    metadata['scaleFactor'] = 10000.0
    #Extract bad band windows
    metadata['bad_band_window1'] = np.array([1340, 1445])
    metadata['bad_band_window2'] = np.array([1790, 1955])
    #Extract map information: spatial extent & resolution (pixel size)
    mapInfo = hdf5_file['map info'].value
    mapInfo_string = str(mapInfo); 
    mapInfo_split = mapInfo_string.split(",")
    #Extract projection information
    metadata['epsg'] = 4326
    #metadata['projection'] = 
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
    
    bandCleaned = reflArray[clipIndex['yMin']:clipIndex['yMax'],clipIndex['xMin']:clipIndex['xMax'],bandIndex-1].astype(np.int16)
    
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


  
    
import numpy as np
import h5py
import gdal, osr
import matplotlib.pyplot as plt
import sys
import ogr, os

f = sys.argv[1]
pt = sys.argv[2]
xx = sys.argv[3]
yy = sys.argv[4]

full_path = pt+f
refl, refl_md, wavelengths = h5refl2array(full_path)
refl_md['extent']
rgb = np.r_[0:425]
rgb = np.delete(rgb, np.r_[419:426])
rgb = np.delete(rgb, np.r_[283:315])
rgfull_path = pt+f
refl, refl_md, wavelengths = h5refl2array(full_path)
refl_md['extent']
rgb = np.r_[0:425]
rgb = np.delete(rgb, np.r_[419:426])
rgb = np.delete(rgb, np.r_[283:315])
rgb = np.delete(rgb, np.r_[192:210])    
xmin, xmax, ymin, ymax = refl_md['extent']
clipExtent = {}
clipExtent['xMin'] = max(int(xmin/1000)*1000 +1000*int(xx), int(xmin))
clipExtent['yMin'] = max(int(ymin/1000)*1000 +1000*int(yy), int(ymin))
clipExtent['yMax'] = min(int(ymin/1000)*1000 +1000*(int(yy)+1), int(ymax))
clipExtent['xMax'] = min(int(xmin/1000)*1000 +1000*(int(xx)+1), int(xmax))

subInd = calc_clip_index(clipExtent,refl_md['ext_dict']) 
subInd['xMax'] = int(subInd['xMax'])
subInd['xMin'] = int(subInd['xMin'])
subInd['yMax'] = int(subInd['yMax'])
subInd['yMin'] = int(subInd['yMin'])

refl = refl[(subInd['yMin']):subInd['yMax'],(subInd['xMin']):subInd['xMax'],:]

#hcp = stack_subset_bands(refl,refl_md,rgb,subInd)
subArray_rows = subInd['yMax'] - subInd['yMin']
subArray_cols = subInd['xMax'] - subInd['xMin']
hcp = np.zeros((subArray_rows,subArray_cols,len(rgb)),dtype = np.int16)

band_clean_dict = {}
band_clean_names = []
for i in range(len(rgb)):
	band_clean_names.append("b"+str(rgb[i])+"_refl_clean")
	band_clean_dict[band_clean_names[i]] = refl[:,:,rgb[i]].astype(np.int16)
	hcp[...,i] = band_clean_dict[band_clean_names[i]]

sub_meta = refl_md
ii = f.replace(' ','')[:-3].upper()
ii = str(clipExtent['xMin'])+str((clipExtent['yMax']))+'.tif'
ras_dir = "/ufrc/ewhite/s.marconi/Marconi2018/tmp"
array2raster(ii,hcp,sub_meta, clipExtent, ras_dir)
