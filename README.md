# Linking baseflow signatures to catchment attributes and hydrological processes

This repository contains code and data used for the submitted manuscript: Gnann et al. (2020) - Linking baseflow signatures to catchment attributes and hydrological processes.
Please email sebastian.gnann@bristol.ac.uk if you have questions.

- This repository doesn't contain any of the original datasets (links to obtain the datasets can be found in the table below)
- To create the Matlab struc file that contains CAMELS data you can use code from this repository: https://github.com/SebastianGnann/CAMELS_Matlab
- One additional open-source package is needed to create the plots: the Brewerplot package (https://github.com/DrosteEffect/BrewerMap).
- Matlab version: **Matlab R2018a**

## Overview of functions and folders:
- *baseflow_signatures_attribute_plots.m* loads data and creates scatter plots between baseflow signatures and catchment attributes shown in the paper.

- *baseflow_signatures_hydrograph_plots.m* loads data and creates hydrographs shown in the paper. 

- *baseflow_signatures_water_balance.m* estimates regional groundwater flow via the water balance (see Supplement).

- *calculate_signatures.m* calculates baseflow signatures for CAMELS catchments and saves them in a mat-file.

- *calculate_attributes.m* calculates catchment attributes using CAMELS shapefiles and raster datasets created with ArcGIS.

- *ArcGIS*: contains ArcGIS files and Python code used to create raster datasets.

- *Functions*: functions used to plot, calculate signatures, etc.

- *Results*: folder where results (e.g. catchment attributes) are saved.

- *Images*: folder where plots are saved.

- *new_CAMELS_attributes_and_signatures.xlsx* contains new catchment attributes and hydrological signatures.


## Acknowledgements

Colours are based on www.ColorBrewer.org, by Cynthia A. Brewer, Penn State.

The function *read_ascii_grid.m* was provided by Gemma Coxon.


## Links to data sources

**DATASET** | **URL**
--- | ---
CAMELS | https://ral.ucar.edu/solutions/products/camels
**CAMELS datasets** |
STATSGO | http://www.soilinfo.psu.edu/index.cgi?soil_data&index.html
GLiM |  https://www.geo.uni-hamburg.de/en/geologie/forschung/geochemie/glim.html
GLHYMPS |  https://dataverse.scholarsportal.info/dataset.xhtml?persistentId=doi:10.5683/SP2/DLGXYO
**Additional datasets** |
HydroSHEDS | https://hydrosheds.org/downloads
Generalized Glacial Limit Lines | https://purl.stanford.edu/vz874sc7648
USGS Geological Map | https://pubs.er.usgs.gov/publication/ds1052
Principal Aquifers of the U.S. | https://water.usgs.gov/ogw/aquifer/map.html
National Wetlands Inventory | https://www.fws.gov/wetlands/index.html
MGS Sinkhole Points  | https://apps5.mo.gov/geostrat/
TWDB Major Aquifers  | http://www.twdb.texas.gov/mapping/gisdata.asp
MODIS | http://www.ntsg.umt.edu/project/modis/mod16.php
GLEAM | https://www.gleam.eu/
**Additional datasets not used in the paper** |
GLHYMPS2 | https://dataverse.scholarsportal.info/dataset.xhtml?persistentId=doi:10.5683/SP2/TTJNIU
Global Thickness of Permeable Layers | https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1304
Global Lakes and Wetlands Database | https://www.worldwildlife.org/publications/global-lakes-and-wetlands-database-lakes-and-wetlands-grid-level-3
Randolph Glacier Inventory | https://www.glims.org/RGI/rgi60_dl.html


