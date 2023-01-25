import ee 
from ee_plugin import Map

#*** Start of imports. If edited, may not auto-convert in the playground. ***#
wiosym_mask = ee.Image("users/PichayaMelody123/wiosymmarin"),
    wiosym_roi =

    # shown: False #
    # displayProperties: [
      {
        "type": "rectangle"
      }
    ] #
    ee.Geometry.Polygon(
        [[[7.481914062500019, 16.19710227789944],
          [7.481914062500019, -40.39396681696679],
          [78.32175781250001, -40.39396681696679],
          [78.32175781250001, 16.19710227789944]]], None, False)
#**** End of imports. If edited, may not auto-convert in the playground. ****#
roi= ee.FeatureCollection(wiosym_roi).geometry()

collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2') \
  .select('tropospheric_NO2_column_number_density') \
  .filterDate('2021-01-01', '2021-12-01')


wiosym = ee.Image ('users/PichayaMelody123/wiosymmarin')


print(collection)

band_viz = {
  'min': 0,
  'max': 0.0002,
  'palette': ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
}
Map.addLayer(wiosym)
Map.addLayer(collection.mean().clip(roi), band_viz, 'S5P N02')

Export.image.toDrive({
  'image' : collection.mean().clip(roi),
  'description': 'NO2_wio',
  'maxPixels': 3784216672400,
  'region' : roi})
Map