// JavaScript code to be implemented in Google Earth Engine(c) developed by H.A. Orengo
// to accompany the paper:

// Orengo, H.A. and Petrie, C.A. 2018. 'Multi-Scale Relief Model (MSRM): a new algorithm
// for the analysis of subtle topographic change in digital elevation models'
// published in Earth Surface Processes and Landforms, 43(6): 1361-1369.

//                ------------- ooo -------------

// TO EXECUTE THE ALGORITHM PASTE THIS CODE INTO GOOGLE EARTH ENGINE CODE EDITOR AND PRESS 'Run'

//                ------------- ooo -------------

// For more information on how this algortihm works and how to apply them refer to the text of the 
// article, which is freely available at: https://onlinelibrary.wiley.com/doi/full/10.1002/esp.4317
// Suggestions and code improvementents are welcome! Please, contact Hector A. Orengo

// NOTES, READ CAREFULLY!
// 1. Visualisation parameters for each layer can be adjusted in the Layers dialogue in the map
//    area of Goggle Earth Engine. This might be convenient when visualising at different scales.
// 2. To apply these analyses in other areas simply change the 'geometry' variable and the
//    'Map.setCenter' coordinates. These two parameters can be deleted. Doing so will apply the
//    algorithms in whatever area is being displayed in the Google Earth Engine's Map area.
//    However, if the geometry variable is eliminated it will need to be also deleted from the
//    scripts incorporating it (for example when it is employed to clip the analysis area).
//    The user can also draw his/her own polygon using the tool provided in the top left of Google
//    Earth Engine's map window. This will create a polygon named 'geometry' by default. The user
//    can then delete the variables called 'geometry' and 'geometry_xprt' in (code lines 46-53 and
//    59-66) and the analysis will be performed in the newly defined area. Make sure that you (1)
//    create a second smaller area named 'geometry_xprt' to be used to clip the first polygon and
//    avoid edge effects and (2) the central coordinates defined by 'Map.setCenter' (line 42) are
//    also changed so the map window will zoom to these coordinates when executing the scripts. 
// 3. The raster resulting from this script can be transferred to your Google Drive running the
//    analysis from the 'task' tab on the top right of the screen.
// 4. Specific instructions of how the algorithm works have been included within the code, please,
//    read these carefully and consult the paper's text to understand how to best apply it.

//                ------------- ooo -------------

// Central coordinates of the Area of Interest (AOI) and zoom level for the map visualisation window,
// please change accordingly when applying MSRM to a different area.
Map.setCenter(74.7125, 29.3511, 8);

// Geometry of the area (defined by WGS84 coordinate pairs) in which to perform the analysis
// The user can define the coordinates of her/his AOI below.
var geometry = ee.Geometry.Polygon(
      [[[77.36572265625,30.230594564932193],
      [76.805419921875,30.590637026892942],
      [76.48681640625,31.024694128525134],
      [72.916259765625,31.00586290462423],
      [72.92724609375,27.98470011861268],
      [77.354736328125,27.994401411046173],
      [77.36572265625,30.230594564932193]]]);

// Geometry of the area (defined by WGS84 coordinate pairs) to cut from the original AOI
// to avoid an edge effect. The user should define the coordinates of her/his AOI below.
// Please, note that an edge effect will still be visible in the map window as GEE map
// visualisation works at screen resolution.
var geometry_xprt = ee.Geometry.Polygon(
      [[[77.27783203125,28.062285999812186],
      [77.288818359375,30.192618218499273],
      [76.739501953125,30.543338954230222],
      [76.44287109375,30.958768570779874],
      [73.004150390625,30.93992433102344],
      [72.9931640625,28.05259082333986],
      [77.27783203125,28.062285999812186]]]);

// Below are the sources imported from GEE repository. The user could add her/his own
// asssets or other sources available in GEE repository.
var SRTM30 = ee.Image('USGS/SRTMGL1_003');
var SRTM90 = ee.Image('CGIAR/SRTM90_V4');
var ALOS30 = ee.Image('JAXA/ALOS/AW3D30_V1_1').select('AVE');

// Please, select the source you want to use from the list of variables above (ALOS30
// selected by default).
var DSM = ALOS30 // Substitute for one of the previously defined DSMs or use your own
  .clip(geometry); // The DSM is clipped using the extent of the AOI defined above.

// The user will be asked about the size of the features s/he wants to detect and the
// Scaling factor to use (please, refer to the text of the article for more information).
var fmax = prompt('What is maximum size of the feature you want to detect in meters?');
var fminPrmpt = prompt('What is minimum size of the feature you want to detect in meters?');
var x = prompt('Which scaling factor do you want to use?');

// The following line of code will extract the raster resolution (rr).
var rr = Math.round((DSM.projection().nominalScale().getInfo())*1000) / 1000;


// The following lines of code will make sure that a minimum feature size is selected that
// can be computed by the algorithm.
if (fminPrmpt <= rr) {
  var fmin = rr;
} else {
  var fmin = fminPrmpt;
}

// Output of information in the Console window.
print('Maximun feature size',fmax,
  'Minimun feature size',fmin,
  'Raster resolution',rr,
  'Scaling factor',x);

// Calculation of 'i' and 'n' values and output of the results in the Console window.
var i = Math.floor(Math.pow((fmin-rr)/(2*rr),1/x));
print('value for i: index of the summation and the radius for the Low Pass ' +
  'filter to be used in the generation of the initial surface to be employed', i);
var n = Math.ceil(Math.pow((fmax-rr)/(2*rr),1/x));
print('Value for n: ', n, 'Number of relief surfaces to be employed', n-i);

// Creation of an array for the storage of the filtered surfaces.
var arrayLP = [];

// Variable for the output of information on the Low Pass filter radii.
var LPrs = '';

// Generation of filtered surfaces using a 'for' loop.
for (var ndx = i; ndx <= n; ndx++) {
  var boxcar = ee.Kernel.square({
  radius: Math.pow(ndx,x), units: 'pixels', normalize: true});
  var LP = DSM.convolve(boxcar);
  LPrs = LPrs + '' + Math.pow(ndx,x) + ',';
  arrayLP.push(LP);
}

// Output of information on the radii of the filtered surfaces generated.
print('Low Pass filter radii to be employed', LPrs);

// The array of filtered images is transformed in a multiband image.
var compLP = ee.Image([arrayLP]);

// Array for the storage of relief models resulting from the subtraction of consecutive
// filtered surfaces.
var arrayRM = [];

// Subtraction of consecutive filtered surfaces for the creation of Relief Models.
for (var ndx2 = 0; ndx2 < n-i; ndx2++) {
  var RM = compLP.expression('((b1 - b2))', {
      'b1': compLP.select(ndx2),
      'b2': compLP.select(ndx2+1)
  });
  arrayRM.push(RM); // Addition of the Relief Models in the previously created array
}

// Average of the relief models.
var MSRMraw = ee.Image([arrayRM])
.reduce(ee.Reducer.mean())
.clip(geometry_xprt); // The resulting raster is clipped using the clipping area defined above
                      // to avoid an edge effect.

// Reduction in the number of decimal places of the values of the resulting raster
// This will not reduce noticeably the quality of the data but it will reduce significantly
// the size of the resulting raster.
var MSRM = ee.Image(0).expression(
    'round(img * 1000) / 1000', {
      'img': MSRMraw
    });

// Statistical data on the range of values of the resulting raster, which will
// be used to strecth the MSRM.
// Please, note these reducers are computationally intensive and the scale might need to be
// increased for larger areas or rasters with higher resolution (he default provided equals
// 4 times the raster size, rr*4). Scale level can also be reduced (up to the raster
// resolution or 'rr') to increase the accuracy of the analysis but this will increase
// processing time and GEE might run out of memory.
var mean = MSRM.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometry,
      scale: rr*4,
      bestEffort: true,
      maxPixels: 1e9,
      tileScale: 10
  })
  .get('mean');
var sigma = MSRM.reduceRegion({
      reducer: ee.Reducer.stdDev(),
      geometry: geometry,
      scale: rr*4,
      bestEffort: true,
      maxPixels: 1e9,
      tileScale: 10
  })
  .get('mean');

// Stretch with using the resulting stats to normalize the image so that 2*sigma fits
// within [0:1].
// Please, note that the exported MSRM is not stretched (so the user can employ her/
// his preferred method).
var MSRMstretch = ee.Image(0).expression(
    '((img - mean) / (sigma * 2)) + 0.5', {
      'img': MSRM,
      'mean': ee.Image.constant(mean),
      'sigma': ee.Image.constant(sigma)
    });

// Generation of a shaded relief model from the stretched MSRM.
var shad =  ee.Terrain.hillshade(MSRMstretch.multiply(1000));

// Add the MSRM (with a 40% transparency) and the shaded relief model to the map view.
Map.addLayer(shad, {}, 'Shaded relief');
Map.addLayer(MSRMstretch,
  {palette: "000000,2207ff,00fff3,03ff00,fbff00,ffc800,ff0000,e300ff,ffffff",
  opacity: 0.4}, 'MSRM');

// Export the unstretched MSRM. Please, use the task window to run the analysis and get the
// resulting image in the user's Google Drive. Please note that the shaded relief model is
// not exported.
Export.image.toDrive({
  image: MSRM,
  description: 'MSRM_x' + Math.round(x) + '_fmn' + Math.round(fmin) + '_fmx' + Math.round(fmax) + '_rr' + Math.round(rr) + 'm',
  scale: rr,
  maxPixels: 9e12, 
  region: geometry_xprt,
});
