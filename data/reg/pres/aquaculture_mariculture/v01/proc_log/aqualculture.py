# -*- coding: utf-8 -*-
"""
Generated by ArcGIS ModelBuilder on : 2022-06-29 12:52:29
"""
import arcpy
from arcpy.ia import *
from arcpy.ia import *
from arcpy.ia import *

def aqualculture():  # aqualculture

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("ImageAnalyst")

    WIOSym_aquaculture_mariculture_selection = "WIOSym_aquaculture_mariculture selection"
    Chl_a_norm_int_tif = arcpy.Raster("Chl_a_norm_int.tif")
    v1_flow_acc_pres_norm_int_tif = arcpy.Raster("v1_flow_acc_pres_norm_int.tif")
    aq_buff_10km_shp_aq_buff_10km = "ModelBuilder\\aq_buff_10km.shp:aq_buff_10km"
    grid_1km_v01_1_tif_2_ = arcpy.Raster("grid_1km_v01.1.tif")
    proc = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc"

    # Process: Project (Project) (management)
    WIOSym_aquaculture_maricultu = "C:\\Users\\piczer\\AppData\\Local\\Temp\\ArcGISProTemp25096\\6a101358-0096-4787-9370-6b647554a71e\\Default.gdb\\WIOSym_aquaculture_maricultu"
    arcpy.management.Project(in_dataset=WIOSym_aquaculture_mariculture_selection, out_dataset=WIOSym_aquaculture_maricultu, out_coor_system="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", transform_method=[], in_coor_system="PROJCS[\"WGS_1984_EASE-Grid_2.0_Global\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",0.0],PARAMETER[\"Standard_Parallel_1\",30.0],UNIT[\"Meter\",1.0]]", preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")

    # Process: Buffer (Buffer) (analysis)
    aq_buff_10km_shp = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc\\aq_buff_10km.shp"
    with arcpy.EnvManager(outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]"):
        arcpy.analysis.Buffer(in_features=WIOSym_aquaculture_maricultu, out_feature_class=aq_buff_10km_shp, buffer_distance_or_field="10 Kilometers", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field=[], method="PLANAR")

    # Process: Weighted Sum (Weighted Sum) (ia)
    chla_alow_acc_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc\\chla_alow_acc.tif"
    Weighted_Sum = chla_alow_acc_tif
    chla_alow_acc_tif = arcpy.ia.WeightedSum(in_rasters=[[Chl_a_norm_int_tif, "Value", 6], [v1_flow_acc_pres_norm_int_tif, "Value", 4]])
    chla_alow_acc_tif.save(Weighted_Sum)


    # Process: Raster Calculator (Raster Calculator) (ia)
    chla_flow_acc_norm_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc\\chla_flow_acc_norm.tif"
    Raster_Calculator = chla_flow_acc_norm_tif
    with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        chla_flow_acc_norm_tif =  ((chla_alow_acc_tif-chla_alow_acc_tif.minimum)/ (chla_alow_acc_tif.maximum - chla_alow_acc_tif.minimum))*100
        chla_flow_acc_norm_tif.save(Raster_Calculator)


    # Process: Polygon to Raster (Polygon to Raster) (conversion)
    aquaculture_10km_buffer_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc\\aquaculture_10km_buffer.tif"
    with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
        arcpy.conversion.PolygonToRaster(in_features=aq_buff_10km_shp_aq_buff_10km, value_field="Id", out_rasterdataset=aquaculture_10km_buffer_tif, cell_assignment="MAXIMUM_COMBINED_AREA", priority_field="Id", cellsize="M:\\marin\\swoc\\work\\wiosym\\data\\reg\\grid\\grid\\v01\\grid_1km_v01.1.tif", build_rat="BUILD")

    # Process: Extract by Mask (Extract by Mask) (sa)
    aquaculture_eutro_relativ_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc\\aquaculture_eutro_relativ.tif"
    Extract_by_Mask = aquaculture_eutro_relativ_tif
    with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        aquaculture_eutro_relativ_tif = arcpy.sa.ExtractByMask(in_raster=chla_flow_acc_norm_tif, in_mask_data=aquaculture_10km_buffer_tif)
        aquaculture_eutro_relativ_tif.save(Extract_by_Mask)


    # Process: Raster Calculator (2) (Raster Calculator) (ia)
    aquaculture_eutro_relativ_1_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\proc\\aquaculture_eutro_relativ_+1.tif"
    Raster_Calculator_2_ = aquaculture_eutro_relativ_1_tif
    aquaculture_eutro_relativ_1_tif =  aquaculture_eutro_relativ_tif+1
    aquaculture_eutro_relativ_1_tif.save(Raster_Calculator_2_)


    # Process: Mosaic To New Raster (Mosaic To New Raster) (management)
    aquaculture_combine_grid_tif = arcpy.management.MosaicToNewRaster(input_rasters=[aquaculture_eutro_relativ_1_tif, grid_1km_v01_1_tif_2_], output_location=proc, raster_dataset_name_with_extension="aquaculture_combine_grid.tif", coordinate_system_for_the_raster="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", pixel_type="32_BIT_FLOAT", cellsize=None, number_of_bands=1, mosaic_method="SUM", mosaic_colormap_mode="FIRST")[0]
    aquaculture_combine_grid_tif = arcpy.Raster(aquaculture_combine_grid_tif)

    # Process: Raster Calculator (3) (Raster Calculator) (ia)
    aquaculture_norm_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\aquaculture\\v01\\aquaculture_norm.tif"
    Raster_Calculator_3_ = aquaculture_norm_tif
    with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                          snapRaster="grid_1km_v01.1.tif"):
        aquaculture_norm_tif =  ((aquaculture_combine_grid_tif- aquaculture_combine_grid_tif.minimum)/(aquaculture_combine_grid_tif.maximum-aquaculture_combine_grid_tif.minimum))*100
        aquaculture_norm_tif.save(Raster_Calculator_3_)


if __name__ == '__main__':
    # Global Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"M:\marin\swoc\work\wiosym\data\reg\pres\aquaculture\v01\proc\aquaculture\Default.gdb", workspace=r"M:\marin\swoc\work\wiosym\data\reg\pres\aquaculture\v01\proc\aquaculture\Default.gdb"):
        aqualculture()