//Region of interest selection
var roi = ee.FeatureCollection(table);
Map.addLayer(roi, {}, 'ROI')
Map.centerObject(roi, 8)

// Image selection and filtering
var R1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate('2018-08-01','2018-08-15')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))

Map.addLayer(R1.first(),{bands: 'VH',min: -20, max: 0}, 'SAR image')

// speckle filter to remove noise effect
var filter_Speckles = function(img) {
  var vh = img.select('VH') //select the VV polarization band
  var vh_smoothed = vh.focal_median(100,'circle','meters').rename('VH_Filtered') //Apply a focal median filter
  return img.addBands(vh_smoothed) // Add filtered VH band to original image
}

R1 = R1.map(filter_Speckles)

// Adding speckle filtered image 
Map.addLayer(R1.first(),{bands: 'VH_Filtered',min: -20, max: 0}, 'Filtered SAR image')

//Classification of Waterbody (flood)
var classify_Water = function(img) {
  var vh = img.select('VH_Filtered')
  var water = vh.lt(-14.5).rename('Water')  //Identify all pixels below threshold and keep them equal to 1. All other pixels set to 0
  water = water.updateMask(water) //Remove all pixels equal to 0
  return img.addBands(water)  //Return image with added classified water band
}

//Map classification of sentinel-1 
R1 = R1.map(classify_Water)
print(R1)


var classification = ee.Image(R1.first()).clip(roi).select('Water'); 
var SARimage = ee.Image(R1.first());

var R1_Layer = ui.Map.Layer(SARimage, {
  bands: ['VH'],
  max: 0,
  min: -20
});
Map.layers().reset([R1_Layer]);
 var visParams = {
  min: 0,
  max: 1,
  palette: ['#FFFFFF','#0000FF']
}
//Add water classification on top of SAR image
Map.addLayer(classification,visParams,'Water')
    
Export.image.toDrive({
  image: classification, 
  description: 'Flood_extent_raster',
  fileNamePrefix: 'flooded',
  scale: 10,
  region: roi, 
  maxPixels: 1e10
});

// Export flooded area as shapefile (for further analysis in e.g. QGIS or ArcGIS)
// Convert flood raster to polygons
var classification_vec = classification.reduceToVectors({
  scale: 10,
  geometryType:'polygon',
  geometry: roi,
  eightConnected: false,
  bestEffort:true,
  tileScale:2,
});

// Export flood polygons as shape-file
Export.table.toDrive({
  collection:classification_vec,
  description:'Flood_extent_vector',
  fileFormat:'SHP',
  fileNamePrefix:'flooded_vec'
});
