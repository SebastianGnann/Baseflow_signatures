%% calculate_attributes
%   
%   Calculates catchment attributes using CAMELS shapefiles and raster
%   datasets created with ArcGIS.
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

clc
% clear all
% close all

%% add directories for functions to path

if exist('./BrewerMap') == 7
    addpath(genpath('./BrewerMap'));
else
    error('BrewerMap toolbox needed. Can be downloaded from https://github.com/DrosteEffect/BrewerMap and should be in a folder named BrewerMap in the same directory.')
end

if exist('./CAMELS_Matlab') == 7
    addpath(genpath('./CAMELS_Matlab')); % CAMELS data
else
    error('CAMELS data processing repository needed. Can be downloaded from https://github.com/SebastianGnann/CAMELS_Matlab and should be in a folder named CAMELS_Matlab in the same directory.')
end

% if exist('./topotoolbox') == 7
%     addpath(genpath('./topotoolbox'));
% end

% if exist('./arcgridwrite') == 7
%     addpath(genpath('./arcgridwrite'));
% end

%% load catchments
if exist('data_CAMELS_struc') % check if data is already loaded
else
    data_CAMELS_struc = load('./CAMELS_Matlab/Data/CAMELS_data.mat');
end
CAMELS_data = data_CAMELS_struc.CAMELS_data;
clear data_CAMELS_struc

%% load files
% load shapefiles
path = "C:\Users\sg16200\Documents\GIS"; % LOCAL path with raster files etc.

% pre-processing in ArcGIS:
%   - coordinate system: WGS84 Web Mercator
%   - catchment boundaries as well as corresponding rivers are shapefiles
%   - all other files (data) are raster files (all snapped to the same
%   grid) with approx. 1km size
%   - for some small catchments there is no river/data available

% load catchment shapefiles
CAMELS_catchments_filename = strcat(path,"\CAMELS\CAMELS_catchments_WGS1984_projected.shp");
CAMELS_catchments = shaperead(CAMELS_catchments_filename);

% load river shapefiles
GloRiC_filename = strcat(path,"\Data_transformed\hydrosheds\riv_CAMELS_dissolved.shp");
GloRiC = shaperead(GloRiC_filename);

%% loop over attributes
for attribute_name = ["sinkholes"]
    %         "Edwards","Pelletier","elevation",...
    %         "permeability", "porosity", "permeability2", "porosity2",...
    %         "glaciers","Pre_Wisconsin","Wisconsin","gleam","modis",...
    %         "wetlands","glwd","mean_geo_age","min_geo_age","sinkholes"
    switch attribute_name
        case 'elevation'
            raster_filename = strcat(path,"\Data_transformed\dem_30s.txt"); % Hydrosheds
        case 'flow_accumulation'
            raster_filename = strcat(path,"\Data_transformed\acc_30s.txt"); % Hydrosheds
        case 'flow_direction'
            raster_filename = strcat(path,"\Data_transformed\dir_30s.txt"); % Hydrosheds
        case 'permeability'
            raster_filename = strcat(path,"\Data_transformed\perm1_30s.txt"); % GLHYMPS 1
        case 'porosity'
            raster_filename = strcat(path,"\Data_transformed\poro1_30s.txt"); % GLHYMPS 1
        case 'permeability2'
            raster_filename = strcat(path,"\Data_transformed\perm2_30s.txt"); % GLHYMPS 2
        case 'porosity2'
            raster_filename = strcat(path,"\Data_transformed\poro2_30s.txt"); % GLHYMPS 2
        case 'mean_geo_age'
            raster_filename = strcat(path,"\Data_transformed\mean_age_30s.txt"); % USGS geology map
        case 'min_geo_age'
            raster_filename = strcat(path,"\Data_transformed\min_age_30s.txt"); % USGS geology map
        case 'glwd'
            raster_filename = strcat(path,"\Data_transformed\glwd3_30s.txt"); % GLWD3
        case 'wetlands'
            raster_filename = strcat(path,"\Data_transformed\test\wetlands_30s.txt"); % wetlands US
        case 'sinkholes'
            raster_filename = strcat(path,"\Data_transformed\sinkholes_30s.txt"); % Missouri sinkhole map
        case 'modis'
            raster_filename = strcat(path,"\Data_transformed\modis_m_30s.txt"); % MODIS ET
        case 'gleam'
            raster_filename = strcat(path,"\Data_transformed\gleam_30s.txt"); % GLEAM ET
        case 'glaciers'
            raster_filename = strcat(path,"\Data_transformed\glaciers_30s.txt"); % Randolph glacier database
        case 'Wisconsin'
            raster_filename = strcat(path,"\Data_transformed\glac_Wi_30s.txt"); % Winsconsin glaciation
        case 'Pre_Wisconsin'
            raster_filename = strcat(path,"\Data_transformed\glac_PWi_30s.txt"); % Pre-Winsconsin glaciation
        case 'Pelletier'
            raster_filename = strcat(path,"\Data_transformed\pelletier_30s.txt"); % Upland/Lowland classification
        case 'Edwards'
            raster_filename = strcat(path,"\Data_transformed\Edwards_30s.txt"); % Edwards-Trinity aquifer system
        case 'xxx'
            raster_filename = strcat(path,"\Data_transformed\xxx.txt"); % xxx
        otherwise
            error('Attribute does not exist.')
    end
    
    % read ascii files
    [data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = read_ascii_grid(raster_filename); % from Gemma
    % create coordinates (2d)
    x_vec = xllcorner + cellsize.*(0:ncols-1);% + cellsize/2; % starts at lower left corner
    y_vec = yllcorner + cellsize.*(nrows-1:-1:0);% + cellsize/2;
    [X,Y] = meshgrid(x_vec,y_vec);
    data(data==-9999) = NaN;
    
    % data transformation
    switch attribute_name
        case 'permeability'
            data = 10.^data;
        case 'permeability2'
            data = 10.^(data./100);
        case 'porosity2'
            data = data./100;
        case 'glwd'
            % recategorise lakes and wetlands
            % class 1: lakes
            data(data==1) = -1; % Lake
            data(data==2) = -1; % Reservoir
            % class 2: rivers
            data(data==3) = -2; % River
            % class 3: other
            data(data==4) = -3; % Freshwater Marsh, Floodplain
            data(data==5) = -3; % Swamp Forest, Flooded Forest
            data(data==6) = -3; % Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)
            data(data==7) = -3; % Pan, Brackish/Saline Wetland
            data(data==8) = -3; % Bog, Fen, Mire (Peatland)
            data(data==9) = -3; % Intermittent Wetland/Lake
            % class 4: major wetland
            data(data==10) = -4; % 50-100% Wetland
            % class 5: minor wetland
            data(data==11) = -5; % 25-50% Wetland
            % class 6: trace wetland
            data(data==12) = -6; % Wetland Complex (0-25% Wetland)
            data = abs(data);
        case 'wetlands'
            % recategorise lakes and wetlands
            % class 1: estuarine_fraction
            data(data==1) = -1; % estuarine and marine deepwater
            data(data==2) = -1; % estuarine and marine wetland
            % class 2: freshwater_fraction
            data(data==3) = -2; % freshwater emergent wetland
            data(data==4) = -2; % freshwater forested/shrub wetland
            data(data==5) = -2; % freshwater pond
            % class 3: lake_fraction
            data(data==6) = -3; % lake
            % class 4: other
            data(data==7) = -4; % other
            data(data==8) = -4; % riverine
            data = abs(data);
        case 'modis'
            data = data./10;
            data(data>6000) = NaN;
        case 'glaciers'
            data (~isnan(data)) = 1; % glacier
        case 'Edwards'
            % recategorise Edwards-Trinity aquifer system
            data(data==5) = -1; % Edwards Trinity
            data(data==8) = -2; % Trinity
            data(data==9) = -3; % Edwards BFZ
            % class 4: other
            data(data==1) = -4; % Seymour
            data(data==2) = -4; % Ogallala
            data(data==3) = -4; % Hueco Bolson
            data(data==4) = -4; % Gulf coast
            data(data==6) = -4; % Pecos valley
            data(data==7) = -4; % Carrizo
            data = abs(data);
    end
    
    % convert to vectors with all coordinates (representing the center)
    x_flat = reshape(X,nrows*ncols,1) + cellsize/2; % get vectors with all coordinates
    y_flat = reshape(Y,nrows*ncols,1) + cellsize/2;
    data_flat = reshape(data,nrows*ncols,1);
    
    %% check which cells are inside catchment and calculate indices
    % preallocate arrays
    n_catchments = length(CAMELS_catchments);
    % special statistics
    switch attribute_name
        case 'elevation'
            hypsometric_integral = NaN(n_catchments,1);
            relief = NaN(n_catchments,1);
            attribute_mean = NaN(n_catchments,1);
            attribute_stdev = NaN(n_catchments,1);
            %         case 'mean_age'
            %             young_fraction = NaN(n_catchments,1);
            %         case 'min_age'
            %             young_fraction = NaN(n_catchments,1);
        case 'glwd'
            lake_fraction = NaN(n_catchments,1);
            river_fraction = NaN(n_catchments,1);
            other_fraction = NaN(n_catchments,1);
            major_wetland_fraction = NaN(n_catchments,1);
            minor_wetland_fraction = NaN(n_catchments,1);
            trace_wetland_fraction = NaN(n_catchments,1);
        case 'wetlands'
            estuarine_fraction = NaN(n_catchments,1);
            freshwater_fraction = NaN(n_catchments,1);
            lake_fraction = NaN(n_catchments,1);
            other_fraction = NaN(n_catchments,1);
        case 'glaciers'
            glacier_fraction = NaN(n_catchments,1);
        case 'Wisconsin'
            Wisconsin_fraction = NaN(n_catchments,1);
        case 'Pre_Wisconsin'
            Pre_Wisconsin_fraction = NaN(n_catchments,1);
        case 'Pelletier'
            upland_fraction = NaN(n_catchments,1);
            lowland_fraction = NaN(n_catchments,1);
            lake_fraction = NaN(n_catchments,1);
            ice_fraction = NaN(n_catchments,1);
        case 'Edwards'
            edwards_trinity_frac = NaN(n_catchments,1);
            trinity_frac = NaN(n_catchments,1);
            edwards_bfz_frac = NaN(n_catchments,1);
            other_frac = NaN(n_catchments,1);
        otherwise
            % basic statistics
            attribute_mean = NaN(n_catchments,1);
            attribute_median = NaN(n_catchments,1);
            attribute_sum = NaN(n_catchments,1);
            attribute_geomean = NaN(n_catchments,1);
            attribute_harmmean = NaN(n_catchments,1);
            attribute_stdev = NaN(n_catchments,1);
            %             attribute_fraction = NaN(n_catchments,1);
    end
    attribute_isnan = NaN(n_catchments,1);
    
    % loop over all catchments
    fprintf('\n%s\n',attribute_name)
    
    for i = 1:n_catchments
        
        if mod(i,10) == 0 % check progress
            fprintf('%.0f/%.0f\n',i,n_catchments)
        end
        
        % extract catchment shapes and check which area is inside the catchment
        catchment_tmp = CAMELS_catchments(i);
        
        % check which cells are inside
        x_flat_tmp = x_flat;
        y_flat_tmp = y_flat;
        data_flat_tmp = data_flat;
        [in,on] = inpolygon(x_flat,y_flat,catchment_tmp.X,catchment_tmp.Y);
        % approximate since cellsize is about 1sqkm. Could make grid finer
        % (e.g. 100 cells) and then count how many are in polygon
        x_flat_tmp(~in) = [];
        y_flat_tmp(~in) = [];
        data_flat_tmp(~in) = [];
        
        % plot catchment
        %{
        plotCatchmentShapefile(catchment_tmp,GloRiC,...
            x_flat_tmp,y_flat_tmp,data_flat_tmp,attribute_name,i,false);
        %}
        
        % calculate statistics
        attribute_isnan(i) = sum(isnan(data_flat_tmp))./length(data_flat_tmp);
        
        % calculate special statistics
        switch attribute_name
            case 'elevation'
                hypsometric_integral(i) = ...
                    (nanmean(data_flat_tmp)-min(data_flat_tmp))./...
                    (max(data_flat_tmp)-min(data_flat_tmp));
                relief(i) = max(data_flat_tmp)-min(data_flat_tmp);
                attribute_mean(i) = nanmean(data_flat_tmp);
                attribute_stdev(i) = nanstd(data_flat_tmp);
                %             case 'mean_age'
                %                 young_fraction(i) = sum(data_flat_tmp<2)./length(data_flat_tmp);
                %             case 'min_age'
                %                 young_fraction(i) = sum(data_flat_tmp<2)./length(data_flat_tmp);
            case 'glwd'
                lake_fraction(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
                river_fraction(i) = sum(data_flat_tmp==2)./length(data_flat_tmp);
                other_fraction(i) = sum(data_flat_tmp==3)./length(data_flat_tmp);
                major_wetland_fraction(i) = sum(data_flat_tmp==4)./length(data_flat_tmp);
                minor_wetland_fraction(i) = sum(data_flat_tmp==5)./length(data_flat_tmp);
                trace_wetland_fraction(i) = sum(data_flat_tmp==6)./length(data_flat_tmp);
            case 'wetlands'
                estuarine_fraction(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
                freshwater_fraction(i) = sum(data_flat_tmp==2)./length(data_flat_tmp);
                lake_fraction(i) = sum(data_flat_tmp==3)./length(data_flat_tmp);
                other_fraction(i) = sum(data_flat_tmp==4)./length(data_flat_tmp);
            case 'glaciers'
                glacier_fraction(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
            case 'Wisconsin'
                Wisconsin_fraction(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
            case 'Pre_Wisconsin'
                Pre_Wisconsin_fraction(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
            case 'Pelletier'
                upland_fraction(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
                lowland_fraction(i) = sum(data_flat_tmp==2)./length(data_flat_tmp);
                lake_fraction(i) = sum(data_flat_tmp==3)./length(data_flat_tmp);
                ice_fraction(i) = sum(data_flat_tmp==4)./length(data_flat_tmp);
            case 'Edwards'
                edwards_trinity_frac(i) = sum(data_flat_tmp==1)./length(data_flat_tmp);
                trinity_frac(i) = sum(data_flat_tmp==2)./length(data_flat_tmp);
                edwards_bfz_frac(i) = sum(data_flat_tmp==3)./length(data_flat_tmp);
                other_frac(i) = sum(data_flat_tmp==4)./length(data_flat_tmp);
            otherwise
                % calculate basic statistics
                if attribute_isnan(i) == 1
                    % if every entry is NaN, then attribute is NaN
                else
                    attribute_mean(i) = nanmean(data_flat_tmp);
                    attribute_median(i) = nanmedian(data_flat_tmp);
                    attribute_sum(i) = nansum(data_flat_tmp);
                    attribute_geomean(i) = geomean(data_flat_tmp(~isnan(data_flat_tmp)));
                    attribute_harmmean(i) = harmmean(data_flat_tmp(~isnan(data_flat_tmp)));
                    attribute_stdev(i) = nanstd(data_flat_tmp);
                end
        end
        
    end
    
    %% save results
    data = {};
    data.attribute_isnan = attribute_isnan;
    
    % special statistics
    switch attribute_name
        case 'elevation'
            data.hypsometric_integral = hypsometric_integral;
            data.relief = relief;
            data.attribute_mean = attribute_mean;
            data.attribute_stdev = attribute_stdev;
            %         case 'mean_age'
            %             data.young_fraction = young_fraction;
            %         case 'min_age'
            %             data.young_fraction = young_fraction;
        case 'glwd'
            data.lake_fraction = lake_fraction;
            data.river_fraction = river_fraction;
            data.other_fraction = other_fraction;
            data.major_wetland_fraction = major_wetland_fraction;
            data.minor_wetland_fraction = minor_wetland_fraction;
            data.trace_wetland_fraction = trace_wetland_fraction;
        case 'wetlands'
            data.estuarine_fraction = estuarine_fraction;
            data.freshwater_fraction = freshwater_fraction;
            data.lake_fraction = lake_fraction;
            data.other_fraction = other_fraction;
        case 'glaciers'
            data.glacier_fraction = glacier_fraction;
        case 'Wisconsin'
            data.Wisconsin_fraction = Wisconsin_fraction;
        case 'Pre_Wisconsin'
            data.Pre_Wisconsin_fraction = Pre_Wisconsin_fraction;
        case 'Pelletier'
            data.Upland_fraction = upland_fraction;
            data.Lowland_fraction = lowland_fraction;
            data.Lake_fraction = lake_fraction;
            data.Ice_fraction = ice_fraction;
        case 'Edwards'
            data.edwards_trinity_frac = edwards_trinity_frac;
            data.trinity_frac = trinity_frac;
            data.edwards_bfz_frac = edwards_bfz_frac;
            data.other_frac = other_frac;
        otherwise
            % basic statistics
            data.attribute_mean = attribute_mean;
            data.attribute_median = attribute_median;
            data.attribute_sum = attribute_sum;
            data.attribute_geomean = attribute_geomean;
            data.attribute_harmmean = attribute_harmmean;
            data.attribute_stdev = attribute_stdev;
    end
    
    % save struc
    save(strcat('./Baseflow_signatures/Results/',attribute_name,'.mat'),'data') 
    
end
