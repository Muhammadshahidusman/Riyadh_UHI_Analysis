//================================================//
//    *** Define the Study Area and Display***     //
//================================================//
// Define a study area of interest (Note: assets are public )
var roi = ee.FeatureCollection('projects/ee-shahidusmaninfo/assets/Riyadh_city');

Map.centerObject(roi, 8);

// Display study area 
Map.addLayer(roi, {}, "roi");

// ✅ Updated to year 2025
var startDate = '2025-06-01';
var endDate   = '2025-08-31';

// *****************************************************************************************

// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

// Cloud mask
function maskL8sr(col) {
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask      = (1 << 5);
  var qa   = col.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return col.updateMask(mask);
}

// Filter the collection
var image = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterDate(startDate, endDate)
  .filterBounds(roi)
  .map(applyScaleFactors)
  .map(maskL8sr)
  .median()
  .clip(roi);

var visualization = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
Map.addLayer(image, visualization, 'True Color (432)', false);

// NDVI
var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
Map.addLayer(ndvi, {min: -1, max: 1, palette: ['blue', 'white', 'green']}, 'ndvi', false);

// NDVI statistics
var ndvi_min = ee.Number(ndvi.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e9
}).values().get(0));

var ndvi_max = ee.Number(ndvi.reduceRegion({
  reducer: ee.Reducer.max(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e9
}).values().get(0));

// Fraction of vegetation
var fv = (ndvi.subtract(ndvi_min).divide(ndvi_max.subtract(ndvi_min)))
           .pow(ee.Number(2))
           .rename('FV');

// Emissivity
var em = fv.multiply(ee.Number(0.004)).add(ee.Number(0.986)).rename('EM');

var thermal = image.select('ST_B10').rename('thermal');

// LST
var lst = thermal.expression(
  '(tb / (1 + ((11.5 * (tb / 14380)) * log(em)))) - 273.15',
  {
    'tb': thermal.select('thermal'),
    'em': em
  }
).rename('LST');

var lst_vis = {
  min: 7,
  max: 50,
  palette: [
    '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
    '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
    '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
    'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
    'ff0000', 'de0101', 'c21301', 'a71001', '911003'
  ]
};
Map.addLayer(lst, lst_vis, 'LST AOI');
Map.centerObject(roi, 10);

// Urban Heat Island ***************************************************************

var lst_mean = ee.Number(lst.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e9
}).values().get(0));

var lst_std = ee.Number(lst.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e9
}).values().get(0));

print('Mean LST in AOI', lst_mean);
print('STD LST in AOI', lst_std);

var uhi = lst.subtract(lst_mean).divide(lst_std).rename('UHI');
var uhi_vis = {
  min: -4,
  max: 4,
  palette: ['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
};
Map.addLayer(uhi, uhi_vis, 'UHI AOI');

// UTFVI
var utfvi = lst.subtract(lst_mean).divide(lst).rename('UTFVI');
var utfvi_vis = {
  min: -1,
  max: 0.3,
  palette: ['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
};
Map.addLayer(utfvi, utfvi_vis, 'UTFVI AOI');

// *********************************************************************************

// ✅ Updated export descriptions to 2025
Export.image.toDrive({
  image: ndvi,
  description: 'ndvi_2025',
  folder: 'Riyadh',
  region: roi,
  scale: 30,
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13,
});

Export.image.toDrive({
  image: uhi,
  description: 'uhi_2025',
  folder: 'Riyadh',
  region: roi,
  scale: 30,
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13,
});

Export.image.toDrive({
  image: lst,
  description: 'lst_2025',
  folder: 'Riyadh',
  region: roi,
  scale: 30,
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13,
});