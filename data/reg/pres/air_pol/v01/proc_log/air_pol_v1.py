# -*- coding: utf-8 -*-
"""
Generated by ArcGIS ModelBuilder on : 2022-06-18 13:41:09
"""
import arcpy
from arcpy.ia import *
from arcpy.ia import *
from arcpy.ia import *
from arcpy.ia import *

def air_pol_v1():  # air_pol_v1

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("ImageAnalyst")

    # Model Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb", workspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb"):
        SO2_wio21_tif = arcpy.Raster("SO2_wio21.tif")
        grid_1km_v01_1_tif = arcpy.Raster("grid_1km_v01.1.tif")
        CO_wio21_1km_tif_2_ = arcpy.Raster("CO_wio21_1km.tif")
        NO2_wio_tif = arcpy.Raster("NO2_wio.tif")

        # Process: Project Raster (6) (Project Raster) (management)
        SO2_project_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\SO2_project.tif"
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", nodata="MAP_UP", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              pyramid="PYRAMIDS -1 NEAREST LZ77 75 NO_SKIP NO_SIPS", scratchWorkspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb", snapRaster="grid_1km_v01.1.tif", 
                              workspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb"):
            arcpy.management.ProjectRaster(in_raster=SO2_wio21_tif, out_raster=SO2_project_tif, out_coor_system="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", resampling_type="NEAREST", cell_size="1000 1000", geographic_transform=[], Registration_Point="", in_coor_system="GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]", vertical="NO_VERTICAL")
            SO2_project_tif = arcpy.Raster(SO2_project_tif)

        # Process: Extract by Mask (2) (Extract by Mask) (sa)
        SO2_project_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\SO2_project_mask.tif"
        Extract_by_Mask_2_ = SO2_project_mask_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            SO2_project_mask_tif = arcpy.sa.ExtractByMask(in_raster=SO2_project_tif, in_mask_data=grid_1km_v01_1_tif)
            SO2_project_mask_tif.save(Extract_by_Mask_2_)


        # Process: Focal Statistics (2) (Focal Statistics) (ia)
        SO2_project_mask_focal_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\SO2_project_mask_focal.tif"
        Focal_Statistics_2_ = SO2_project_mask_focal_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\Air_pol_v01.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            SO2_project_mask_focal_tif = arcpy.ia.FocalStatistics(in_raster=SO2_project_mask_tif, neighborhood="Wedge 50 0 360 CELL", statistics_type="PERCENTILE", ignore_nodata="DATA", percentile_value=80)
            SO2_project_mask_focal_tif.save(Focal_Statistics_2_)


        # Process: Raster Calculator (Raster Calculator) (ia)
        SO2_project_mask_focal_norm_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\SO2_project_mask_focal_norm.tif"
        Raster_Calculator = SO2_project_mask_focal_norm_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask=r"ModelBuilder\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            SO2_project_mask_focal_norm_tif =  ((SO2_project_mask_focal_tif - SO2_project_mask_focal_tif.minimum)/(SO2_project_mask_focal_tif.maximum - SO2_project_mask_focal_tif.minimum))*100
            SO2_project_mask_focal_norm_tif.save(Raster_Calculator)


        # Process: Project Raster (2) (Project Raster) (management)
        CO_project_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\CO_project.tif"
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", nodata="MAP_UP", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              pyramid="PYRAMIDS -1 NEAREST LZ77 75 NO_SKIP NO_SIPS", scratchWorkspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb", snapRaster="grid_1km_v01.1.tif", 
                              workspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb"):
            arcpy.management.ProjectRaster(in_raster=CO_wio21_1km_tif_2_, out_raster=CO_project_tif, out_coor_system="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", resampling_type="NEAREST", cell_size="1000 1000", geographic_transform=[], Registration_Point="", in_coor_system="GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]", vertical="NO_VERTICAL")
            CO_project_tif = arcpy.Raster(CO_project_tif)

        # Process: Extract by Mask (3) (Extract by Mask) (sa)
        CO_project_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\CO_project_mask.tif"
        Extract_by_Mask_3_ = CO_project_mask_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            CO_project_mask_tif = arcpy.sa.ExtractByMask(in_raster=CO_project_tif, in_mask_data=grid_1km_v01_1_tif)
            CO_project_mask_tif.save(Extract_by_Mask_3_)


        # Process: Focal Statistics (3) (Focal Statistics) (ia)
        CO_project_mask_focal_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\CO_project_mask_focal.tif"
        Focal_Statistics_3_ = CO_project_mask_focal_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\Air_pol_v01.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            CO_project_mask_focal_tif = arcpy.ia.FocalStatistics(in_raster=CO_project_mask_tif, neighborhood="Wedge 50 0 360 CELL", statistics_type="PERCENTILE", ignore_nodata="DATA", percentile_value=80)
            CO_project_mask_focal_tif.save(Focal_Statistics_3_)


        # Process: Raster Calculator (2) (Raster Calculator) (ia)
        CO_project_mask_focal_norm_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\CO_project_mask_focal_norm.tif"
        Raster_Calculator_2_ = CO_project_mask_focal_norm_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask=r"ModelBuilder\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            CO_project_mask_focal_norm_tif =  ((CO_project_mask_focal_tif - CO_project_mask_focal_tif.minimum)/(CO_project_mask_focal_tif.maximum - CO_project_mask_focal_tif.minimum))*100
            CO_project_mask_focal_norm_tif.save(Raster_Calculator_2_)


        # Process: Project Raster (7) (Project Raster) (management)
        NO2_project_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\NO2_project.tif"
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", nodata="MAP_UP", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              pyramid="PYRAMIDS -1 NEAREST LZ77 75 NO_SKIP NO_SIPS", scratchWorkspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb", snapRaster="grid_1km_v01.1.tif", 
                              workspace=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\air_pol\air_pol.gdb"):
            arcpy.management.ProjectRaster(in_raster=NO2_wio_tif, out_raster=NO2_project_tif, out_coor_system="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", resampling_type="NEAREST", cell_size="1000 1000", geographic_transform=[], Registration_Point="", in_coor_system="GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]", vertical="NO_VERTICAL")
            NO2_project_tif = arcpy.Raster(NO2_project_tif)

        # Process: Extract by Mask (4) (Extract by Mask) (sa)
        NO2_project_mask_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\NO2_project_mask.tif"
        Extract_by_Mask_4_ = NO2_project_mask_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            NO2_project_mask_tif = arcpy.sa.ExtractByMask(in_raster=NO2_project_tif, in_mask_data=grid_1km_v01_1_tif)
            NO2_project_mask_tif.save(Extract_by_Mask_4_)


        # Process: Focal Statistics (4) (Focal Statistics) (ia)
        NO2_project_mask_focal_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\NO2_project_mask_focal.tif"
        Focal_Statistics_4_ = NO2_project_mask_focal_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\pres\air_pol\v01\proc\Air_pol_v01.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            NO2_project_mask_focal_tif = arcpy.ia.FocalStatistics(in_raster=NO2_project_mask_tif, neighborhood="Wedge 50 0 360 CELL", statistics_type="PERCENTILE", ignore_nodata="DATA", percentile_value=80)
            NO2_project_mask_focal_tif.save(Focal_Statistics_4_)


        # Process: Raster Calculator (3) (Raster Calculator) (ia)
        NO2_project_mask_focal_norm_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\NO2_project_mask_focal_norm.tif"
        Raster_Calculator_3_ = NO2_project_mask_focal_norm_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask=r"ModelBuilder\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            NO2_project_mask_focal_norm_tif =  ((NO2_project_mask_focal_tif - NO2_project_mask_focal_tif.minimum)/(NO2_project_mask_focal_tif.maximum - NO2_project_mask_focal_tif.minimum))*100
            NO2_project_mask_focal_norm_tif.save(Raster_Calculator_3_)


        # Process: Weighted Sum (Weighted Sum) (ia)
        air_pol_combine_weighted_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\proc\\air_pol\\air_pol_combine_weighted.tif"
        Weighted_Sum = air_pol_combine_weighted_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask=r"ModelBuilder\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            air_pol_combine_weighted_tif = arcpy.ia.WeightedSum(in_rasters=[[SO2_project_mask_focal_norm_tif, "VALUE", 3], [CO_project_mask_focal_norm_tif, "VALUE", 2], [NO2_project_mask_focal_norm_tif, "VALUE", 5]])
            air_pol_combine_weighted_tif.save(Weighted_Sum)


        # Process: Raster Calculator (4) (Raster Calculator) (ia)
        Air_pol_norm_v01_0_tif = "M:\\proj\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\air_pol\\v01\\Air_pol_norm_v01.0_.tif"
        Raster_Calculator_4_ = Air_pol_norm_v01_0_tif
        with arcpy.EnvManager(cellSize=r"M:\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask=r"ModelBuilder\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            Air_pol_norm_v01_0_tif =  ((air_pol_combine_weighted_tif-air_pol_combine_weighted_tif.minimum)/(air_pol_combine_weighted_tif.maximum- air_pol_combine_weighted_tif.minimum))*100
            Air_pol_norm_v01_0_tif.save(Raster_Calculator_4_)


if __name__ == '__main__':
    air_pol_v1()