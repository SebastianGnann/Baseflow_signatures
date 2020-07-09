# Transform all wetland shapefiles to rasters and merge them to get a single raster file
# Note that the Wetlands have all been added to a layer "Wetland_States" first. 
# You can change the path if you want to load them directly.

import arcpy

# Environment settings are set in ArcGIS. 
# Output coordinate system is WGS1984 Web Mercator.
# Grid is snapped to dem_30s.
# Everything else is default.

# loop over all states or state regions
in_polygon_list = ["Wetlands_States\WY_Wetlands_West",
                   "Wetlands_States\WY_Wetlands_East",
                   "Wetlands_States\WI_Wetlands_South",
                   "Wetlands_States\WI_Wetlands_North",
                   "Wetlands_States\TX_Wetlands_West",
                   "Wetlands_States\TX_Wetlands_East",
                   "Wetlands_States\TX_Wetlands_Central",
                   "Wetlands_States\SD_Wetlands_West",
                   "Wetlands_States\SD_Wetlands_East",
                   "Wetlands_States\OR_Wetlands_West",
                   "Wetlands_States\OR_Wetlands_East",
                   "Wetlands_States\OK_Wetlands_West",
                   "Wetlands_States\OK_Wetlands_East",
                   "Wetlands_States\NV_Wetlands_South",
                   "Wetlands_States\NV_Wetlands_North",
                   "Wetlands_States\NE_Wetlands_West",
                   "Wetlands_States\NE_Wetlands_East",
                   "Wetlands_States\ND_Wetlands_West",
                   "Wetlands_States\ND_Wetlands_East",
                   "Wetlands_States\MT_Wetlands_West",
                   "Wetlands_States\MT_Wetlands_East",
                   "Wetlands_States\North_Wetlands_North_West",
                   "Wetlands_States\North_Wetlands_North_East",
                   "Wetlands_States\MN_Wetlands_South",
                   "Wetlands_States\MN_Wetlands_Central_West",
                   "Wetlands_States\MN_Wetlands_Central_East",
                   "Wetlands_States\KS_Wetlands_West",
                   "Wetlands_States\KS_Wetlands_East",
                   "Wetlands_States\IA_Wetlands_West",
                   "Wetlands_States\IA_Wetlands_East",
                   "Wetlands_States\CO_Wetlands_West",
                   "Wetlands_States\CO_Wetlands_East",
                   "Wetlands_States\IL_Wetlands",
                   "Wetlands_States\ME_Wetlands",
                   "Wetlands_States\MI_Wetlands",
                   "Wetlands_States\OH_Wetlands",
                   "Wetlands_States\KY_Wetlands",
                   "Wetlands_States\PA_Wetlands",
                   "Wetlands_States\DC_Wetlands",
                   "Wetlands_States\WV_Wetlands",
                   "Wetlands_States\IN_Wetlands",
                   "Wetlands_States\VT_Wetlands",
                   "Wetlands_States\NJ_Wetlands",
                   "Wetlands_States\MD_Wetlands",
                   "Wetlands_States\MA_Wetlands",
                   "Wetlands_States\NH_Wetlands",
                   "Wetlands_States\CT_Wetlands",
                   "Wetlands_States\DE_Wetlands",
                   "Wetlands_States\RI_Wetlands",
                   "Wetlands_States\WA_Wetlands",
                   "Wetlands_States\MO_Wetlands",
                   "Wetlands_States\LA_Wetlands",
                   "Wetlands_States\ID_Wetlands",
                   "Wetlands_States\MS_Wetlands",
                   "Wetlands_States\GA_Wetlands",
                   "Wetlands_States\UT_Wetlands",
                   "Wetlands_States\NM_Wetlands",
                   "Wetlands_States\FL_Wetlands",
                   "Wetlands_States\NC_Wetlands",
                   "Wetlands_States\AR_Wetlands",
                   "Wetlands_States\TN_Wetlands",
                   "Wetlands_States\VA_Wetlands",
                   "Wetlands_States\SC_Wetlands",
                   "Wetlands_States\AL_Wetlands",
                   "Wetlands_States\AZ_Wetlands",
                   "Wetlands_States\CA_Wetlands_North",
                   "Wetlands_States\CA_Wetlands_NorthCentral",
                   "Wetlands_States\CA_Wetlands_South",
                   "Wetlands_States\CA_Wetlands_SouthCentral"]

out_raster_list = ["C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wy1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wy2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wi1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wi2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/tx1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/tx2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/tx3_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/sd1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/sd2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/or1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/or2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ok1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ok2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nv1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nv2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ne1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ne2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nd1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nd2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mt1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mt2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mn1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mn2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mn3_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mn4_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mn5_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ks1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ks2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ia1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ia2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/co1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/co2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/il_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/me_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mi_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/oh_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ky_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/pa_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/dc_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wv_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/in_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/vt_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nj_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/md_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ma_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nh_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ct_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/de_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ri_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wa_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/mo_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/la_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/id_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ms_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ga_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ut_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nm_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/fl_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/nc_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ar_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/tn_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/va_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/sc_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/al_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/az_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ca1_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ca2_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ca3_class_30s",
                   "C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/ca4_class_30s"]

# out_ascii_list = ["C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands/wy1_class_30s.TXT"]

# 34 didn't work, MI

for i in range(1, len(in_polygon_list)):
	
	in_polygon = in_polygon_list[i]
	out_raster = out_raster_list[i]
	# out_ascii = out_raster + ".TXT"
	# out_ascii = out_ascii_list[i]

	# join table to layer, to get consistent classes
	arcpy.AddJoin_management(
		in_layer_or_view=in_polygon, 
		in_field="WETLAND_TY", 
		join_table="C:/Users/sg16200/Documents/GIS/CAMELS_geodatabase.gdb/wetlands", 
		join_field="WETLAND_TY", 
		join_type="KEEP_ALL")
		
	# transform polygon to grid
	arcpy.PolygonToRaster_conversion(
		in_features=in_polygon, 
		value_field="wetlands.CLASS", 
		out_rasterdataset=out_raster, 
		cell_assignment="MAXIMUM_COMBINED_AREA", 
		priority_field="NONE", 
		cellsize="C:\Users\sg16200\Documents\GIS\Data_transformed\hydrosheds\dem_30s")
		
	# transform grid to ascii
	# arcpy.RasterToASCII_conversion(
	# 	in_raster=out_raster, 
	# 	out_ascii_file=out_ascii)

# merge individual grids
arcpy.MosaicToNewRaster_management(
    input_rasters=out_raster_list, 
    output_location="C:/Users/sg16200/Documents/GIS/Data_transformed/",  #CAMELS_geodatabase.gdb
    raster_dataset_name_with_extension="wetlands_30s", 
    coordinate_system_for_the_raster="", 
    pixel_type="8_BIT_UNSIGNED", 
    cellsize="", 
    number_of_bands="1", 
    mosaic_method="LAST", 
    mosaic_colormap_mode="FIRST")

# transform grid to ascii
arcpy.RasterToASCII_conversion(
 	in_raster="C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands_30s", 
 	out_ascii_file="C:/Users/sg16200/Documents/GIS/Data_transformed/wetlands_30s.TXT")