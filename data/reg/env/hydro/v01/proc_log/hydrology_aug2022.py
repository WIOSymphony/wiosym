# -*- coding: utf-8 -*-
"""
Generated by ArcGIS ModelBuilder on : 2022-09-07 23:35:29
"""
import arcpy
from arcpy.ia import *

def hydrology_aug2022():  # hydrology_aug2022

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("3D")
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("ImageAnalyst")

    # Model Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"C:\Users\piczer\AppData\Local\Temp\ArcGISProTemp15916\80f4c96f-4192-46bf-9fa5-034acbd824e9\Default.gdb", workspace=r"C:\Users\piczer\AppData\Local\Temp\ArcGISProTemp15916\80f4c96f-4192-46bf-9fa5-034acbd824e9\Default.gdb"):
        gebco_land_sea_tif = arcpy.Raster("gebco_land_sea.tif")
        gebco_land_shallow_fill_n_tif = arcpy.Raster("gebco_land_shallow_fill_n.tif")
        gebco_bathymetry_tif = arcpy.Raster("gebco_bathymetry.tif")
        hydro_aug2022_update_2_ = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update"
        gebco_land_tif = arcpy.Raster("gebco_land.tif")
        eco_ch_shore_shoreline_proportion_norm01_1km_v01_0_tif_3_ = arcpy.Raster("eco_ch_shore_shoreline_proportion_norm01_1km_v01.0.tif")
        ch_mudflat_v01_0_tif = arcpy.Raster("ch_mudflat_v01.0.tif")
        hydro_aug2022_update = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update"
        grid_1km_v01_1_tif = arcpy.Raster("ModelBuilder\\grid_1km_v01.1.tif")
        flow_acc_grid_tif_flow_acc_grid_tif = arcpy.Raster("ModelBuilder\\flow_acc_grid.tif:flow_acc_grid.tif")
        flow_distance_tif_flow_distance_tif = arcpy.Raster("ModelBuilder\\flow_distance.tif:flow_distance.tif")
        DEM_land_shoreline_fill_flow_acc_tif_DEM_land_shoreline_fill_flow_acc_tif = arcpy.Raster("ModelBuilder\\DEM_land_shoreline_fill_flow_acc.tif:DEM_land_shoreline_fill_flow_acc.tif")

        # Process: Mosaic To New Raster (2) (Mosaic To New Raster) (management)
        gebco_land_bathymetry_tif = arcpy.management.MosaicToNewRaster(input_rasters=[gebco_land_shallow_fill_n_tif, gebco_bathymetry_tif], output_location=hydro_aug2022_update_2_, raster_dataset_name_with_extension="gebco_land_bathymetry.tif", coordinate_system_for_the_raster="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", pixel_type="32_BIT_FLOAT", cellsize=None, number_of_bands=1, mosaic_method="LAST", mosaic_colormap_mode="FIRST")[0]
        gebco_land_bathymetry_tif = arcpy.Raster(gebco_land_bathymetry_tif)

        # Process: Reclassify (Reclassify) (sa)
        shore_line_grid_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\shore_line_grid.tif"
        Reclassify = shore_line_grid_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            shore_line_grid_tif = arcpy.sa.Reclassify(in_raster=eco_ch_shore_shoreline_proportion_norm01_1km_v01_0_tif_3_, reclass_field="VALUE", remap="0 NODATA;1 1", missing_values="DATA")
            shore_line_grid_tif.save(Reclassify)


        # Process: Reclassify (2) (Reclassify) (sa)
        intertidal_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\intertidal.tif"
        Reclassify_2_ = intertidal_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            intertidal_tif = arcpy.sa.Reclassify(in_raster=ch_mudflat_v01_0_tif, reclass_field="VALUE", remap="0 100 1", missing_values="DATA")
            intertidal_tif.save(Reclassify_2_)


        # Process: Mosaic To New Raster (Mosaic To New Raster) (management)
        with arcpy.EnvManager(scratchWorkspace=r"C:\Users\piczer\AppData\Local\Temp\ArcGISProTemp15916\80f4c96f-4192-46bf-9fa5-034acbd824e9\Default.gdb", workspace=r"C:\Users\piczer\AppData\Local\Temp\ArcGISProTemp15916\80f4c96f-4192-46bf-9fa5-034acbd824e9\Default.gdb"):
            shore_line_land_tif = arcpy.management.MosaicToNewRaster(input_rasters=[gebco_land_tif, shore_line_grid_tif, intertidal_tif], output_location=hydro_aug2022_update, raster_dataset_name_with_extension="shore_line_land.tif", coordinate_system_for_the_raster="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", pixel_type="32_BIT_FLOAT", cellsize=None, number_of_bands=1, mosaic_method="SUM", mosaic_colormap_mode="FIRST")[0]
            shore_line_land_tif = arcpy.Raster(shore_line_land_tif)

        # Process: Reclassify (3) (Reclassify) (sa)
        shore_line_land_intertidal_mask1_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\shore_line_land_intertidal_mask1.tif"
        Reclassify_3_ = shore_line_land_intertidal_mask1_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            shore_line_land_intertidal_mask1_tif = arcpy.sa.Reclassify(in_raster=shore_line_land_tif, reclass_field="VALUE", remap="0 3 1", missing_values="DATA")
            shore_line_land_intertidal_mask1_tif.save(Reclassify_3_)


        # Process: Extract by Mask (Extract by Mask) (sa)
        DEM_land_shoreline_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\DEM_land_shoreline.tif"
        Extract_by_Mask = DEM_land_shoreline_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            DEM_land_shoreline_tif = arcpy.sa.ExtractByMask(in_raster=gebco_land_bathymetry_tif, in_mask_data=shore_line_land_intertidal_mask1_tif)
            DEM_land_shoreline_tif.save(Extract_by_Mask)


        # Process: Fill (Fill) (sa)
        DEM_land_shoreline_fill_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\DEM_land_shoreline_fill.tif"
        Fill = DEM_land_shoreline_fill_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            DEM_land_shoreline_fill_tif = arcpy.sa.Fill(in_surface_raster=DEM_land_shoreline_tif, z_limit=None)
            DEM_land_shoreline_fill_tif.save(Fill)


        # Process: Flow Direction (Flow Direction) (sa)
        DEM_land_shoreline_fill_flow_direction_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\DEM_land_shoreline_fill_flow_direction.tif"
        Flow_Direction = DEM_land_shoreline_fill_flow_direction_tif
        Output_drop_raster = ""
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            DEM_land_shoreline_fill_flow_direction_tif = arcpy.sa.FlowDirection(in_surface_raster=DEM_land_shoreline_fill_tif, force_flow="NORMAL", out_drop_raster=Output_drop_raster, flow_direction_type="D8")
            DEM_land_shoreline_fill_flow_direction_tif.save(Flow_Direction)


        # Process: Flow Accumulation (Flow Accumulation) (sa)
        DEM_land_shoreline_fill_flow_acc_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\DEM_land_shoreline_fill_flow_acc.tif"
        Flow_Accumulation = DEM_land_shoreline_fill_flow_acc_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            DEM_land_shoreline_fill_flow_acc_tif = arcpy.sa.FlowAccumulation(in_flow_direction_raster=DEM_land_shoreline_fill_flow_direction_tif, in_weight_raster="", data_type="FLOAT", flow_direction_type="D8")
            DEM_land_shoreline_fill_flow_acc_tif.save(Flow_Accumulation)


        # Process: Focal Statistics (Focal Statistics) (sa)
        DEM_land_shoreline_fill_flow_acc_focal_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\DEM_land_shoreline_fill_flow_acc_focal.tif"
        Focal_Statistics = DEM_land_shoreline_fill_flow_acc_focal_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            DEM_land_shoreline_fill_flow_acc_focal_tif = arcpy.sa.FocalStatistics(in_raster=DEM_land_shoreline_fill_flow_acc_tif, neighborhood="Circle 2 CELL", statistics_type="SUM", ignore_nodata="DATA", percentile_value=90)
            DEM_land_shoreline_fill_flow_acc_focal_tif.save(Focal_Statistics)


        # Process: Extract by Mask (2) (Extract by Mask) (sa)
        flow_acc_grid_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\flow_acc_grid.tif"
        Extract_by_Mask_2_ = flow_acc_grid_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            flow_acc_grid_tif = arcpy.sa.ExtractByMask(in_raster=DEM_land_shoreline_fill_flow_acc_focal_tif, in_mask_data=grid_1km_v01_1_tif)
            flow_acc_grid_tif.save(Extract_by_Mask_2_)


        # Process: Raster Calculator (Raster Calculator) (ia)
        v2_flow_acc_aug2022_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\v2_flow_acc_aug2022_.tif"
        Raster_Calculator = v2_flow_acc_aug2022_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask=r"ModelBuilder\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster=r"ModelBuilder\grid_1km_v01.1.tif"):
            v2_flow_acc_aug2022_tif = SetNull(flow_acc_grid_tif_flow_acc_grid_tif<300, flow_acc_grid_tif_flow_acc_grid_tif)
            v2_flow_acc_aug2022_tif.save(Raster_Calculator)


        # Process: Focal Flow (Focal Flow) (sa)
        focal_flow_flow_acc_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\env\\hydro\\v01\\proc\\hydro_aug2022_update\\focal_flow_flow_acc.tif"
        Focal_Flow = focal_flow_flow_acc_tif
        with arcpy.EnvManager(outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            focal_flow_flow_acc_tif = arcpy.sa.FocalFlow(in_surface_raster=DEM_land_shoreline_fill_flow_acc_tif, threshold_value=0)
            focal_flow_flow_acc_tif.save(Focal_Flow)


if __name__ == '__main__':
    hydrology_aug2022()