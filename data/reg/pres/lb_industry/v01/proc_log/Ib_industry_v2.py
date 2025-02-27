# -*- coding: utf-8 -*-
"""
Generated by ArcGIS ModelBuilder on : 2022-09-06 13:48:47
"""
import arcpy
from arcpy.ia import *
from arcpy.ia import *
from arcpy.ia import *
from arcpy.ia import *
from arcpy.ia import *

def Ib_industry2():  # Ib_industry2

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("ImageAnalyst")

    # Model Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"M:\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"M:\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
        m3_tif = arcpy.Raster("Chl_a_b1mg/m3.tif")
        global_power_plant2021 = "global_power_plant2021"
        wiosym_grid_bounding_box_v01 = "wiosym_grid_bounding_box_v01"
        grid_1km_v01_1_tif = arcpy.Raster("grid_1km_v01.1.tif")
        TSM_mosaic_mg_l_tif = arcpy.Raster("TSM_mosaic_mg_l.tif")
        Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif_3_ = arcpy.Raster("Ib_industry_pol_v.0.1.0_1_km.tif:Ib_industry_pol_v.0.1.0_1_km.tif")
        Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif_4_ = arcpy.Raster("Ib_industry_pol_v.0.1.0_1_km.tif:Ib_industry_pol_v.0.1.0_1_km.tif")
        Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif = arcpy.Raster("Ib_industry_pol_v.0.1.0_1_km.tif:Ib_industry_pol_v.0.1.0_1_km.tif")
        Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif_2_ = arcpy.Raster("Ib_industry_pol_v.0.1.0_1_km.tif:Ib_industry_pol_v.0.1.0_1_km.tif")
        lb_industry = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry"
        industrial_pol = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industrial_pol"
        grid_1km_v01_1_tif_3_ = arcpy.Raster("grid_1km_v01.1.tif")
        industrial_pol_2_ = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industrial_pol"

        # Process: Project Raster (Project Raster) (management)
        chla_projected_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\chla_projected.tif"
        with arcpy.EnvManager(nodata="MAP_DOWN", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_na_v01.1.tif"):
            arcpy.management.ProjectRaster(in_raster=m3_tif, out_raster=chla_projected_tif, out_coor_system="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", resampling_type="NEAREST", cell_size="1000 1000", geographic_transform=[], Registration_Point="", in_coor_system="GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]", vertical="NO_VERTICAL")
            chla_projected_tif = arcpy.Raster(chla_projected_tif)

        # Process: Clip (Clip) (analysis)
        power_plant_clip_wio1_shp = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\power_plant_clip_wio1.shp"
        with arcpy.EnvManager(scratchWorkspace=r"C:\Users\piczer\AppData\Local\Temp\ArcGISProTemp11704\bb82b4f9-c166-400f-be37-6c5dc48b1b16\Default.gdb", workspace=r"C:\Users\piczer\AppData\Local\Temp\ArcGISProTemp11704\bb82b4f9-c166-400f-be37-6c5dc48b1b16\Default.gdb"):
            arcpy.analysis.Clip(in_features=global_power_plant2021, clip_features=wiosym_grid_bounding_box_v01, out_feature_class=power_plant_clip_wio1_shp, cluster_tolerance="")

        # Process: Buffer (Buffer) (analysis)
        power_plant_clip_wio_buffer50km_shp = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\power_plant_clip_wio_buffer50km.shp"
        with arcpy.EnvManager(outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            arcpy.analysis.Buffer(in_features=power_plant_clip_wio1_shp, out_feature_class=power_plant_clip_wio_buffer50km_shp, buffer_distance_or_field="50 Kilometers", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field=[], method="PLANAR")

        # Process: Polygon to Raster (Polygon to Raster) (conversion)
        _100km_buffer_raster_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\100km_buffer_raster.tif"
        arcpy.conversion.PolygonToRaster(in_features=power_plant_clip_wio_buffer50km_shp, value_field="Id", out_rasterdataset=_100km_buffer_raster_tif, cell_assignment="CELL_CENTER", priority_field="NONE", cellsize="100", build_rat="BUILD")

        # Process: Extract by Mask (Extract by Mask) (sa)
        chla_extracted_buffer_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\chla_extracted_buffer.tif"
        Extract_by_Mask = chla_extracted_buffer_tif
        chla_extracted_buffer_tif = arcpy.sa.ExtractByMask(in_raster=chla_projected_tif, in_mask_data=_100km_buffer_raster_tif)
        chla_extracted_buffer_tif.save(Extract_by_Mask)


        # Process: Resample (2) (Resample) (management)
        Chla_resample_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\Chla_resample.tif"
        with arcpy.EnvManager(nodata="MAP_DOWN", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", pyramid="PYRAMIDS -1 NEAREST LZ77 75 NO_SKIP NO_SIPS", 
                              snapRaster="grid_1km_na_v01.1.tif"):
            arcpy.management.Resample(in_raster=chla_extracted_buffer_tif, out_raster=Chla_resample_tif, cell_size="1000 1000", resampling_type="NEAREST")
            Chla_resample_tif = arcpy.Raster(Chla_resample_tif)

        # Process: Extract by Mask (3) (Extract by Mask) (sa)
        chla_grid_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\chla_grid.tif"
        Extract_by_Mask_3_ = chla_grid_tif
        with arcpy.EnvManager(scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            chla_grid_tif = arcpy.sa.ExtractByMask(in_raster=Chla_resample_tif, in_mask_data=grid_1km_v01_1_tif)
            chla_grid_tif.save(Extract_by_Mask_3_)


        # Process: Extract by Mask (2) (Extract by Mask) (sa)
        TSM_extract_buffer_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\TSM_extract_buffer.tif"
        Extract_by_Mask_2_ = TSM_extract_buffer_tif
        TSM_extract_buffer_tif = arcpy.sa.ExtractByMask(in_raster=TSM_mosaic_mg_l_tif, in_mask_data=_100km_buffer_raster_tif)
        TSM_extract_buffer_tif.save(Extract_by_Mask_2_)


        # Process: Raster Calculator (Raster Calculator) (ia)
        TSM_extract_buffer_reverse_value_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\TSM_extract_buffer_reverse_value.tif"
        Raster_Calculator = TSM_extract_buffer_reverse_value_tif
        TSM_extract_buffer_reverse_value_tif =  TSM_extract_buffer_tif *(-1)
        TSM_extract_buffer_reverse_value_tif.save(Raster_Calculator)


        # Process: Resample (Resample) (management)
        TSM_resample_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\TSM_resample.tif"
        with arcpy.EnvManager(scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            arcpy.management.Resample(in_raster=TSM_extract_buffer_reverse_value_tif, out_raster=TSM_resample_tif, cell_size="1000 1000", resampling_type="NEAREST")
            TSM_resample_tif = arcpy.Raster(TSM_resample_tif)

        # Process: Extract by Mask (4) (Extract by Mask) (sa)
        TSM_grid_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\TSM_grid.tif"
        Extract_by_Mask_4_ = TSM_grid_tif
        with arcpy.EnvManager(scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            TSM_grid_tif = arcpy.sa.ExtractByMask(in_raster=TSM_resample_tif, in_mask_data=grid_1km_v01_1_tif)
            TSM_grid_tif.save(Extract_by_Mask_4_)


        # Process: Raster Calculator (2) (Raster Calculator) (ia)
        TSM_grid_changeunit_to_mgm3_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\TSM_grid_changeunit_to_mgm3.tif"
        Raster_Calculator_2_ = TSM_grid_changeunit_to_mgm3_tif
        with arcpy.EnvManager(scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            TSM_grid_changeunit_to_mgm3_tif =  TSM_grid_tif /1000
            TSM_grid_changeunit_to_mgm3_tif.save(Raster_Calculator_2_)


        # Process: Raster Calculator (3) (Raster Calculator) (ia)
        turbidity_proxy_combine_TSM_chla_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\turbidity_proxy_combine_TSM_chla.tif"
        Raster_Calculator_3_ = turbidity_proxy_combine_TSM_chla_tif
        with arcpy.EnvManager(scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            turbidity_proxy_combine_TSM_chla_tif =  chla_grid_tif + TSM_grid_changeunit_to_mgm3_tif
            turbidity_proxy_combine_TSM_chla_tif.save(Raster_Calculator_3_)


        # Process: Extract by Mask (5) (Extract by Mask) (sa)
        Ib_industry_pol_v_0_1_0_1_km_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\Ib_industry_pol_v.0.1.0_1_km.tif"
        Extract_by_Mask_5_ = Ib_industry_pol_v_0_1_0_1_km_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", snapRaster="grid_1km_v01.1.tif", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            Ib_industry_pol_v_0_1_0_1_km_tif = arcpy.sa.ExtractByMask(in_raster=turbidity_proxy_combine_TSM_chla_tif, in_mask_data=grid_1km_v01_1_tif)
            Ib_industry_pol_v_0_1_0_1_km_tif.save(Extract_by_Mask_5_)


        # Process: Raster Calculator (4) (Raster Calculator) (ia)
        industry_pol_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industry_pol.tif"
        Raster_Calculator_4_ = industry_pol_tif
        with arcpy.EnvManager(cellSize=r"\\sgu.se\SGU\prod\proj\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", scratchWorkspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", 
                              snapRaster="grid_1km_v01.1.tif", workspace=r"\\sgu.se\sgu\prod\proj\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            industry_pol_tif =  ((Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif- Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif_2_.minimum)/( Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif_2_.maximum-Ib_industry_pol_v_0_1_0_1_km_tif_Ib_industry_pol_v_0_1_0_1_km_tif_2_.minimum))*100
            industry_pol_tif.save(Raster_Calculator_4_)


        # Process: Mosaic To New Raster (Mosaic To New Raster) (management)
        Ib_industry_pol_v_0_2_0_1km_tif = arcpy.management.MosaicToNewRaster(input_rasters=[industry_pol_tif, grid_1km_v01_1_tif], output_location=lb_industry, raster_dataset_name_with_extension="Ib_industry_pol_v.0.2.0_1km.tif", coordinate_system_for_the_raster="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", pixel_type="32_BIT_FLOAT", cellsize=None, number_of_bands=1, mosaic_method="FIRST", mosaic_colormap_mode="FIRST")[0]
        Ib_industry_pol_v_0_2_0_1km_tif = arcpy.Raster(Ib_industry_pol_v_0_2_0_1km_tif)

        # Process: Project (Project) (management)
        power_plant_wio_reproject_shp = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\power_plant_wio_reproject.shp"
        with arcpy.EnvManager(scratchWorkspace=r"M:\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb", workspace=r"M:\marin\swoc\work\wiosym\data\reg\pres\lb_industry\v01\proc\industrial_pol\Default.gdb"):
            arcpy.management.Project(in_dataset=power_plant_clip_wio1_shp, out_dataset=power_plant_wio_reproject_shp, out_coor_system="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", transform_method=[], in_coor_system="GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]", preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")

        # Process: Point Statistics (Point Statistics) (sa)
        industry_size_pointstate_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industry_size_pointstate.tif"
        Point_Statistics = industry_size_pointstate_tif
        industry_size_pointstate_tif = arcpy.sa.PointStatistics(in_point_features=power_plant_wio_reproject_shp, field="capacity_m", cell_size="1000", neighborhood="Circle 20000 MAP", statistics_type="MEAN")
        industry_size_pointstate_tif.save(Point_Statistics)


        # Process: Raster to Polygon (Raster to Polygon) (conversion)
        grid_1km_polygon_shp = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\grid_1km_polygon.shp"
        with arcpy.EnvManager(outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", outputMFlag="Disabled", outputZFlag="Disabled", 
                              snapRaster="grid_1km_v01.1.tif"):
            arcpy.conversion.RasterToPolygon(in_raster=grid_1km_v01_1_tif_3_, out_polygon_features=grid_1km_polygon_shp, simplify="SIMPLIFY", raster_field="Value", create_multipart_features="SINGLE_OUTER_PART", max_vertices_per_feature=None)

        # Process: Intersect (Intersect) (analysis)
        power_plant_intersect_grid_shp = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\power_plant_intersect_grid.shp"
        arcpy.analysis.Intersect(in_features=[[grid_1km_polygon_shp, ""]], out_feature_class=power_plant_intersect_grid_shp, join_attributes="ALL", cluster_tolerance="", output_type="POINT")

        # Process: Raster Calculator (5) (Raster Calculator) (ia)
        grid_0proxy_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\grid_0proxy.tif"
        Raster_Calculator_5_ = grid_0proxy_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            grid_0proxy_tif =  grid_1km_v01_1_tif *0
            grid_0proxy_tif.save(Raster_Calculator_5_)


        # Process: Point to Raster (2) (Point to Raster) (conversion)
        industrial_capacity_raster_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industrial_capacity_raster.tif"
        arcpy.conversion.PointToRaster(in_features=power_plant_wio_reproject_shp, value_field="capacity_m", out_rasterdataset=industrial_capacity_raster_tif, cell_assignment="SUM", priority_field="capacity_m", cellsize="10000", build_rat="BUILD")

        # Process: Resample (3) (Resample) (management)
        industrial_capacity_raster_R = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industrial_pol\\Default.gdb\\industrial_capacity_raster_R"
        with arcpy.EnvManager(outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            arcpy.management.Resample(in_raster=industrial_capacity_raster_tif, out_raster=industrial_capacity_raster_R, cell_size="10000 10000", resampling_type="NEAREST")
            industrial_capacity_raster_R = arcpy.Raster(industrial_capacity_raster_R)

        # Process: Mosaic To New Raster (2) (Mosaic To New Raster) (management)
        industrial_cap_mosaic_grid_tif = arcpy.management.MosaicToNewRaster(input_rasters=[grid_0proxy_tif, industrial_capacity_raster_R], output_location=industrial_pol_2_, raster_dataset_name_with_extension="industrial_cap_mosaic_grid.tif", coordinate_system_for_the_raster="PROJCS[\"WGS_1984_Cylindrical_Equal_Area\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",12.0],PARAMETER[\"standard_parallel_1\",-12.0],UNIT[\"Meter\",1.0]]", pixel_type="32_BIT_FLOAT", cellsize=None, number_of_bands=1, mosaic_method="SUM", mosaic_colormap_mode="FIRST")[0]
        industrial_cap_mosaic_grid_tif = arcpy.Raster(industrial_cap_mosaic_grid_tif)

        # Process: Extract by Mask (6) (Extract by Mask) (sa)
        industrail_mask_grid_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\v01\\proc\\industrail_mask_grid.tif"
        Extract_by_Mask_6_ = industrail_mask_grid_tif
        with arcpy.EnvManager(mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", snapRaster="grid_1km_v01.1.tif"):
            industrail_mask_grid_tif = arcpy.sa.ExtractByMask(in_raster=industrial_cap_mosaic_grid_tif, in_mask_data=grid_1km_v01_1_tif_3_)
            industrail_mask_grid_tif.save(Extract_by_Mask_6_)


        # Process: Focal Statistics (Focal Statistics) (ia)
        Ib_industry_pol_v_0_2_1_1km_tif = "M:\\marin\\swoc\\work\\wiosym\\data\\reg\\pres\\lb_industry\\Ib_industry_pol_v.0.2.1_1km.tif"
        Focal_Statistics = Ib_industry_pol_v_0_2_1_1km_tif
        with arcpy.EnvManager(cellSize=r"M:\marin\swoc\work\wiosym\data\reg\grid\grid\v01\grid_1km_v01.1.tif", mask="grid_1km_v01.1.tif", outputCoordinateSystem="PROJCS["WGS_1984_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",12.0],PARAMETER["standard_parallel_1",-12.0],UNIT["Meter",1.0]]", 
                              snapRaster="grid_1km_v01.1.tif"):
            Ib_industry_pol_v_0_2_1_1km_tif = arcpy.ia.FocalStatistics(in_raster=industrail_mask_grid_tif, neighborhood="Wedge 5 0 360 CELL", statistics_type="MEAN", ignore_nodata="DATA", percentile_value=90)
            Ib_industry_pol_v_0_2_1_1km_tif.save(Focal_Statistics)


if __name__ == '__main__':
    Ib_industry2()
