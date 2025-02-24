# -*- coding: utf-8 -*-
"""
Generated by ArcGIS ModelBuilder on : 2022-06-18 15:52:58
"""
import arcpy
from arcpy.ia import *
from arcpy.ia import *

def ship():  # ship

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("ImageAnalyst")

    shipdensity_global_grid1km_tif = arcpy.Raster("shipdensity_global_grid1km.tif")
    grid_1km_v01_1_tif = arcpy.Raster("grid_1km_v01.1.tif")
    ShipDensity_Commercial1_grid1km_tif = arcpy.Raster("ShipDensity_Commercial1_grid1km.tif")
    ShipDensity_Fishing1_grid1km_tif = arcpy.Raster("ShipDensity_Fishing1_grid1km.tif")
    ShipDensity_OilGas1_grid1km_tif = arcpy.Raster("ShipDensity_OilGas1_grid1km.tif")
    ShipDensity_Leisure1_grid1km_tif = arcpy.Raster("ShipDensity_Leisure1_grid1km.tif")
    ShipDensity_Passenger1_grid1km_tif = arcpy.Raster("ShipDensity_Passenger1_grid1km.tif")
    proc = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc"

    # Process: Extract by Mask (4) (Extract by Mask) (sa)
    ship_all_global_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_all_global_mask.tif"
    Extract_by_Mask_4_ = ship_all_global_mask_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_all_global_mask_tif = arcpy.sa.ExtractByMask(in_raster=shipdensity_global_grid1km_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_all_global_mask_tif.save(Extract_by_Mask_4_)


    # Process: Extract by Mask (2) (Extract by Mask) (sa)
    ship_commer_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_commer_mask.tif"
    Extract_by_Mask_2_ = ship_commer_mask_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_commer_mask_tif = arcpy.sa.ExtractByMask(in_raster=ShipDensity_Commercial1_grid1km_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_commer_mask_tif.save(Extract_by_Mask_2_)


    # Process: Extract by Mask (3) (Extract by Mask) (sa)
    ship_fishing_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_fishing_mask.tif"
    Extract_by_Mask_3_ = ship_fishing_mask_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_fishing_mask_tif = arcpy.sa.ExtractByMask(in_raster=ShipDensity_Fishing1_grid1km_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_fishing_mask_tif.save(Extract_by_Mask_3_)


    # Process: Extract by Mask (5) (Extract by Mask) (sa)
    ship_oilgas_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_oilgas_mask.tif"
    Extract_by_Mask_5_ = ship_oilgas_mask_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_oilgas_mask_tif = arcpy.sa.ExtractByMask(in_raster=ShipDensity_OilGas1_grid1km_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_oilgas_mask_tif.save(Extract_by_Mask_5_)


    # Process: Extract by Mask (6) (Extract by Mask) (sa)
    ship_leisure_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_leisure_mask.tif"
    Extract_by_Mask_6_ = ship_leisure_mask_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_leisure_mask_tif = arcpy.sa.ExtractByMask(in_raster=ShipDensity_Leisure1_grid1km_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_leisure_mask_tif.save(Extract_by_Mask_6_)


    # Process: Extract by Mask (Extract by Mask) (sa)
    ship_passenger_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_passenger_mask.tif"
    Extract_by_Mask = ship_passenger_mask_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_passenger_mask_tif = arcpy.sa.ExtractByMask(in_raster=ShipDensity_Passenger1_grid1km_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_passenger_mask_tif.save(Extract_by_Mask)


    # Process: Mosaic To New Raster (Mosaic To New Raster) (management)
    with arcpy.EnvManager(pyramid="PYRAMIDS -1 NEAREST LZ77 75 NO_SKIP NO_SIPS"):
        ship_combine_mosaic_tif = arcpy.management.MosaicToNewRaster(input_rasters=[ship_commer_mask_tif, ship_fishing_mask_tif, ship_oilgas_mask_tif, ship_leisure_mask_tif, ship_passenger_mask_tif], output_location=proc, raster_dataset_name_with_extension="ship_combine_mosaic.tif", coordinate_system_for_the_raster="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", pixel_type="2_BIT", cellsize=None, number_of_bands=1, mosaic_method="SUM", mosaic_colormap_mode="FIRST")[0]
        ship_combine_mosaic_tif = arcpy.Raster(ship_combine_mosaic_tif)

    # Process: Raster Calculator (Raster Calculator) (ia)
    ship_combine_mosaic_norm_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_combine_mosaic_norm.tif"
    Raster_Calculator = ship_combine_mosaic_norm_tif
    ship_combine_mosaic_norm_tif =  ((ship_combine_mosaic_tif -ship_combine_mosaic_tif.minimum)/(ship_combine_mosaic_tif.maximum - ship_combine_mosaic_tif.minimum))*(1)
    ship_combine_mosaic_norm_tif.save(Raster_Calculator)


    # Process: Focal Statistics (Focal Statistics) (ia)
    ship_mosaic_combine_focal_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_mosaic_combine_focal.tif"
    Focal_Statistics = ship_mosaic_combine_focal_tif
    ship_mosaic_combine_focal_tif = arcpy.ia.FocalStatistics(in_raster=ship_combine_mosaic_norm_tif, neighborhood="Annulus 25 50 CELL", statistics_type="SUM", ignore_nodata="DATA", percentile_value=80)
    ship_mosaic_combine_focal_tif.save(Focal_Statistics)


    # Process: Extract by Mask (7) (Extract by Mask) (sa)
    ship_pollution_ralative_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\proc\\ship_pollution_ralative_.tif"
    Extract_by_Mask_7_ = ship_pollution_ralative_tif
    with arcpy.EnvManager(mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
        ship_pollution_ralative_tif = arcpy.sa.ExtractByMask(in_raster=ship_mosaic_combine_focal_tif, in_mask_data=grid_1km_v01_1_tif)
        ship_pollution_ralative_tif.save(Extract_by_Mask_7_)


    # Process: Raster Calculator (2) (Raster Calculator) (ia)
    ship_pollution_norm_1km_v01_1_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\act\\ship\\v01\\ship_pollution_norm_1km_v01.1.tif"
    Raster_Calculator_2_ = ship_pollution_norm_1km_v01_1_tif
    with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        ship_pollution_norm_1km_v01_1_tif =  ((ship_pollution_ralative_tif -ship_pollution_ralative_tif.minimum)/( ship_pollution_ralative_tif.maximum- ship_pollution_ralative_tif.minimum))*100
        ship_pollution_norm_1km_v01_1_tif.save(Raster_Calculator_2_)


if __name__ == '__main__':
    # Global Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\act\ship\v01\proc\ship\Default.gdb", workspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\act\ship\v01\proc\ship\Default.gdb"):
        ship()
