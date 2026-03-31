// ============================================================
// LULC Classification using Random Forest
// Landsat 8 SR Imagery — 4 Classes (2014)
// Classes: Water | Builtup | Vegetation | Barren Land
// Study Area: Riyadh City
// ============================================================


// ============================================================
// STEP 1: DEFINE STUDY AREA
// ============================================================

var aoi = ee.FeatureCollection('projects/ee-shahidusmaninfo/assets/Riyadh_city');

Map.centerObject(aoi, 8);
Map.addLayer(aoi, {color: 'red'}, 'Study Area Boundary');


// ============================================================
// STEP 2: LOAD AND FILTER LANDSAT 8 DATA
// ============================================================

var startDate = '2014-06-01';
var endDate   = '2014-08-31';

// Apply scaling factors (Landsat Collection 2 Level-2)
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

// Cloud & shadow masking using QA_PIXEL
function maskL8sr(image) {
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask      = (1 << 5);
  var qa   = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

var collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterDate(startDate, endDate)
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUD_COVER', 15))
  .map(applyScaleFactors)
  .map(maskL8sr);

print('Total Landsat 8 Scenes:', collection.size());

var composite = collection.median().clip(aoi);

// Visualize True Color (RGB = B4, B3, B2)
Map.addLayer(composite, {
  bands: ['SR_B5', 'SR_B4', 'SR_B3'],
  min: 0.140216,
  max: 0.575637,
  gamma: 1.4
}, 'True Color (432)');

// ============================================================
// STEP 3: CALCULATE SPECTRAL INDICES
// Landsat 8 Band Reference:
//   SR_B2 = Blue  | SR_B3 = Green | SR_B4 = Red
//   SR_B5 = NIR   | SR_B6 = SWIR1 | SR_B7 = SWIR2
// ============================================================

// NDVI — Vegetation index
var ndvi = composite.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

// NDWI — Water index (McFeeters)
var ndwi = composite.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');

// MNDWI — Modified water index (Xu)
var mndwi = composite.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI');

// NDBI — Built-up index
var ndbi = composite.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');

// BSI — Bare Soil Index
var bsi = composite.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))',
  {
    'BLUE':  composite.select('SR_B2'),
    'RED':   composite.select('SR_B4'),
    'NIR':   composite.select('SR_B5'),
    'SWIR1': composite.select('SR_B6')
  }
).rename('BSI');

// UI — Urban Index
var ui_index = composite.normalizedDifference(['SR_B6', 'SR_B4']).rename('UI');

// Stack all bands + indices into one multi-band image
var inputBands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];

var inputImage = composite.select(inputBands)
  .addBands(ndvi)
  .addBands(ndwi)
  .addBands(mndwi)
  .addBands(ndbi)
  .addBands(bsi)
  .addBands(ui_index);

print('Input image bands:', inputImage.bandNames());


// ============================================================
// STEP 4: DRAW TRAINING SAMPLES IN GEE
//
//   HOW TO DRAW:
//   1. In the GEE Code Editor, click the geometry tools (top-left of map)
//   2. Draw polygons for each class
//   3. In the "Geometry Imports" panel, rename each geometry:
//      → water        (for Class 0)
//      → builtup      (for Class 1)
//      → vegetation   (for Class 2)
//      → barrenLand   (for Class 3)
//   4. Set each geometry type to "Feature Collection" (not Geometry)
//
//   CLASS TABLE:
//   ┌──────┬──────────────┬─────────┐
//   │ Code │ Class Name   │ Color   │
//   ├──────┼──────────────┼─────────┤
//   │  0   │ Water        │ Blue    │
//   │  1   │ Builtup      │ Red     │
//   │  2   │ Vegetation   │ Green   │
//   │  3   │ Barren Land  │ Yellow  │
//   └──────┴──────────────┴─────────┘
//
//   TIP: Draw at least 10–20 polygons per class for better accuracy.
//   TIP: For Riyadh 2014, Water is rare — focus on seasonal streams/reservoirs.
// ============================================================

var waterFC      = ee.FeatureCollection([ee.Feature(water,      {landcover: 0, className: 'Water'})]);
var builtupFC    = ee.FeatureCollection([ee.Feature(builtup,    {landcover: 1, className: 'Builtup'})]);
var vegetationFC = ee.FeatureCollection([ee.Feature(vegetation, {landcover: 2, className: 'Vegetation'})]);
var barrenLandFC = ee.FeatureCollection([ee.Feature(barrenLand, {landcover: 3, className: 'BarrenLand'})]);


// ============================================================
// STEP 5: MERGE ALL TRAINING DATA
// ============================================================

var trainingData = waterFC
  .merge(builtupFC)
  .merge(vegetationFC)
  .merge(barrenLandFC);

print('=== TRAINING DATA SUMMARY ===');
print('Total training polygons:', trainingData.size());
print('Class distribution:', trainingData.aggregate_histogram('landcover'));

// Uncomment to visualize training polygons:
// Map.addLayer(waterFC,      {color: '0000FF'}, 'Train — Water');
// Map.addLayer(builtupFC,    {color: 'FF0000'}, 'Train — Builtup');
// Map.addLayer(vegetationFC, {color: '00AA00'}, 'Train — Vegetation');
// Map.addLayer(barrenLandFC, {color: 'FFFF00'}, 'Train — Barren Land');


// ============================================================
// STEP 6: SAMPLE THE INPUT IMAGE
// ============================================================

var trainingSamples = inputImage.sampleRegions({
  collection: trainingData,
  properties: ['landcover', 'className'],
  scale:      30,       // Landsat native resolution
  tileScale:  16
});

print('Total training samples extracted:', trainingSamples.size());


// ============================================================
// STEP 7: SPLIT INTO TRAIN (70%) AND VALIDATION (30%)
// ============================================================

var withRandom = trainingSamples.randomColumn('random', 42);

var trainSet = withRandom.filter(ee.Filter.lte('random', 0.7));
var testSet  = withRandom.filter(ee.Filter.gt('random',  0.7));

print('Train set size:', trainSet.size());
print('Test set size:',  testSet.size());


// ============================================================
// STEP 8: TRAIN RANDOM FOREST CLASSIFIER
// ============================================================

print('=== TRAINING RANDOM FOREST MODEL ===');

var classifier = ee.Classifier.smileRandomForest({
  numberOfTrees:     150,
  variablesPerSplit: null,
  minLeafPopulation: 1,
  bagFraction:       0.5,
  maxNodes:          null,
  seed:              42
})
.train({
  features:        trainSet,
  classProperty:   'landcover',
  inputProperties: inputImage.bandNames()
});

print('Model trained successfully ✓');


// ============================================================
// STEP 9: CLASSIFY THE IMAGE
// ============================================================

var classified = inputImage.classify(classifier);

var classVis = {
  min: 0,
  max: 3,
  palette: [
    '0000FF',  // 0 — Water       (Blue)
    'FF0000',  // 1 — Builtup     (Red)
    '00AA00',  // 2 — Vegetation  (Green)
    'FFFF00'   // 3 — Barren Land (Yellow)
  ]
};

Map.addLayer(classified.clip(aoi), classVis, 'LULC Classification — 2014');
print('Classification complete ✓');


// ============================================================
// STEP 10: ACCURACY ASSESSMENT
// ============================================================

print('=== ACCURACY ASSESSMENT ===');

var validated = testSet.classify(classifier);

var confusionMatrix = validated.errorMatrix('landcover', 'classification');

print('Confusion Matrix:',        confusionMatrix);
print('Overall Accuracy:',        confusionMatrix.accuracy());
print('Kappa Coefficient:',       confusionMatrix.kappa());
print("Producer's Accuracy:",     confusionMatrix.producersAccuracy());
print("Consumer's Accuracy:",     confusionMatrix.consumersAccuracy());


// ============================================================
// STEP 11: AREA STATISTICS (sq km per class) — Labeled Output
// ============================================================

print('=== AREA STATISTICS ===');

var classNames = {
  '0': 'Water',
  '1': 'Builtup',
  '2': 'Vegetation',
  '3': 'Barren Land'
};

var areaImage = ee.Image.pixelArea().divide(1e6)  // convert m² → sq km
  .addBands(classified.rename('landcover'));

var areaStats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName:  'landcover'
  }),
  geometry: aoi,
  scale:     30,
  maxPixels: 1e13
});

// Extract the list of groups
var groups = ee.List(areaStats.get('groups'));

// Map over each group and format as a readable dictionary
var labeledAreas = groups.map(function(item) {
  var dict      = ee.Dictionary(item);
  var classNum  = ee.Number(dict.get('landcover')).int();
  var area      = ee.Number(dict.get('sum')).round();

  var className = ee.Dictionary({
    '0': 'Water',
    '1': 'Builtup',
    '2': 'Vegetation',
    '3': 'Barren Land'
  }).get(classNum.format('%d'));

  return ee.Dictionary({
    'Class ID':     classNum,
    'Class Name':   className,
    'Area (sq km)': area
  });
});

print('─────────────────────────────────');
print('Class-wise Area Summary (sq km):');
print('─────────────────────────────────');
print(labeledAreas);
print('─────────────────────────────────');


// ============================================================
// STEP 12: MAP LEGEND
// ============================================================

var legend = ui.Panel({
  style: {
    position:        'bottom-left',
    padding:         '10px 15px',
    backgroundColor: 'white'
  }
});

legend.add(ui.Label({
  value: 'LULC — Riyadh 2014',
  style: { fontWeight: 'bold', fontSize: '15px', margin: '0 0 8px 0' }
}));

var classes = [
  [0, 'Water',       '0000FF'],
  [1, 'Builtup',     'FF0000'],
  [2, 'Vegetation',  '00AA00'],
  [3, 'Barren Land', 'FFFF00']
];

classes.forEach(function(cls) {
  var colorBox = ui.Label({
    style: {
      backgroundColor: '#' + cls[2],
      padding: '10px',
      margin:  '2px 6px 2px 0',
      border:  '1px solid #888'
    }
  });
  var label = ui.Label({
    value: cls[0] + ' — ' + cls[1],
    style: { fontSize: '13px', margin: '4px 0' }
  });
  legend.add(ui.Panel({
    widgets: [colorBox, label],
    layout:  ui.Panel.Layout.Flow('horizontal')
  }));
});

Map.add(legend);


// ============================================================
// STEP 13: EXPORT TO GOOGLE DRIVE
// ============================================================

// Export LULC Classification Map
Export.image.toDrive({
  image:          classified.clip(aoi).toByte(),
  description:    'LULC_Riyadh_2014_RF',
  folder:         'Riyadh',
  fileNamePrefix: 'LULC_Riyadh_2014_RF',
  region:         aoi,
  scale:          30,
  maxPixels:      1e13,
  crs:            'EPSG:4326'
});

print('=== SCRIPT COMPLETE ✓ ===');
print('Go to Tasks tab → click RUN on each export job');