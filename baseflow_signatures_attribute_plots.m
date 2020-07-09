%% baseflow_signatures_attribute_plots
%
%   Loads data and creates scatter plots between baseflow signatures and 
%   catchment attributes shown in the paper.
%   (Contains additional attributes not used in the paper.)
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

clc

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

% figure path
fig_path = '.\Baseflow_signatures\Images';

%% load catchment data
if exist('data_CAMELS_struc') % check if data is already loaded
else
    data_CAMELS_struc = load('./CAMELS_Matlab/Data/CAMELS_data.mat');
end
CAMELS_data = data_CAMELS_struc.CAMELS_data;
clear data_CAMELS_struc

%% load additional signatures
if exist('signatures_CAMELS_struc') % check if data is already loaded
else
    signatures_CAMELS_struc = load('./Baseflow_signatures/Results/CAMELS_signatures.mat');
end
CAMELS_signatures = signatures_CAMELS_struc.CAMELS_signatures;
clear signatures_CAMELS_struc

%% load additional data
% sinkholes
sinkholes = load('sinkholes.mat'); sinkholes = sinkholes.data.attribute_sum;
% formerly glaciated areas
Wisconsin = load('Wisconsin.mat');  Wisconsin_fraction = Wisconsin.data.Wisconsin_fraction;
Pre_Wisconsin = load('Pre_Wisconsin.mat'); Pre_Wisconsin_fraction = Pre_Wisconsin.data.Pre_Wisconsin_fraction;
% geological age
mean_age = load('mean_geo_age.mat'); mean_age = mean_age.data.attribute_mean;
min_age = load('min_geo_age.mat'); min_age = min_age.data.attribute_mean;
% surface water bodies
wetlands = load('wetlands.mat');
freshwater_fraction = wetlands.data.freshwater_fraction;
lake_fraction = wetlands.data.lake_fraction;
estuarine_fraction = wetlands.data.estuarine_fraction;
other_fraction = wetlands.data.other_fraction;
% glaciers
glaciers = load('glaciers.mat'); glacier_fraction = glaciers.data.glacier_fraction;
% geology
porosity = load('porosity.mat'); porosity = porosity.data.attribute_mean;
permeability = load('permeability.mat'); permeability = permeability.data.attribute_geomean;
% porosity2 = load('porosity2.mat'); porosity2 = porosity2.data.attribute_mean;
% permeability2 = load('permeability2.mat'); permeability2 = permeability2.data.attribute_geomean;
% topography
elevation = load('elevation.mat');
hypsometric_integral = elevation.data.hypsometric_integral;
relief = elevation.data.relief;
% Pelletier classification
Pelletier = load('Pelletier.mat');
lowland_fraction = Pelletier.data.Lowland_fraction;
upland_fraction = Pelletier.data.Upland_fraction;
lake_fraction_pelletier = Pelletier.data.Lake_fraction;
ice_fraction = Pelletier.data.Ice_fraction;
% ET
ET_modis = load('modis.mat'); ET_modis = ET_modis.data.attribute_mean;
ET_gleam = load('gleam.mat'); ET_gleam = ET_gleam.data.attribute_mean;
% Edwards-Trinity
Edwards = load('Edwards.mat');
isn = Edwards.data.attribute_isnan;
trinity_frac = Edwards.data.trinity_frac; trinity_frac(isn==1) = NaN;
edwards_trinity_frac = Edwards.data.edwards_trinity_frac; edwards_trinity_frac(isn==1) = NaN;
edwards_bfz_frac = Edwards.data.edwards_bfz_frac; edwards_bfz_frac(isn==1) = NaN;
other_frac = Edwards.data.other_frac; other_frac(isn==1) = NaN;
% translate categorical values (strings) into numeric ones
% Pelletier and other classes
lowland = lowland_fraction>0.5;
upland = upland_fraction>0.5;
snow = CAMELS_data.frac_snow>0.3;
% lakes = (lake_fraction)>0.01;
lakes = (lake_fraction+freshwater_fraction)>0.01;
% freshwater = (freshwater_fraction)>0.1;
% lakes2 = ((freshwater_fraction+lake_fraction)>0.01 & data_CAMELS.clay_frac>20);
glaciers = glacier_fraction>0;

%% glacial deposits
plotTrivariate(-1,-1,-1,...
    'x_name','f{clay} [-]','y_name','f{sand} [-]','z_name','BFI5 [-]',...
    'ID',CAMELS_data.gauge_id,...
    'colour_scheme','YlGNBu','flip_colour_scheme',false,'nr_colours',10,...
    'x_limits',[0 .5],'y_limits',[0 1],...
    'z_limits',[0 1],'z_lower_limit_open',false,'z_upper_limit_open',false,...
    'figure_title','','figure_name','soil_texture',...
    'save_figure',false,'figure_path',fig_path,'figure_type','-dmeta')
% p4=scatter(data_CAMELS.clay_frac./100,data_CAMELS.sand_frac./100,20,...
%     signatures_CAMELS.BFI5,'filled','MarkerFaceColor',[0.8 0.8 0.8]); 
p0=scatter(CAMELS_data.clay_frac(Wisconsin_fraction>0.5)./100,CAMELS_data.sand_frac(Wisconsin_fraction>0.5)./100,20,...
    CAMELS_signatures.BFI5(Wisconsin_fraction>0.5),'filled'); 
p1=scatter(CAMELS_data.clay_frac(Pre_Wisconsin_fraction>0.5&Wisconsin_fraction<0.5)./100,CAMELS_data.sand_frac(Pre_Wisconsin_fraction>0.5&Wisconsin_fraction<0.5)./100,20,...
    CAMELS_signatures.BFI5(Pre_Wisconsin_fraction>0.5&Wisconsin_fraction<0.5),'^','filled'); 
p2=scatter(-1,-1,20,0.6,'filled'); % fake legend
p3=scatter(-1,-1,20,0.3,'^','filled');
legend([p2,p3],{'Wisconsin','Pre Wisconsin'},'box','off','location','northeast','FontSize',7,'ItemTokenSize',[7,7])
% set(gca,'YScale','log')
saveFig(gcf,strcat('trivariate_','soil_texture'),fig_path,'-dmeta')
disp(corr(CAMELS_signatures.BFI5(Pre_Wisconsin_fraction>0.5),CAMELS_data.clay_frac(Pre_Wisconsin_fraction>0.5),'Type','Spearman','Rows','Complete'))
disp(corr(CAMELS_signatures.BFI5(Pre_Wisconsin_fraction>0.5),CAMELS_data.sand_frac(Pre_Wisconsin_fraction>0.5),'Type','Spearman','Rows','Complete'))

%% North Carolina
catchments = [115:117 119 122:136 211:213 239:244 247 250:252]; % catchment indices

plotTrivariate(CAMELS_data.clay_frac(catchments)./100,...
    CAMELS_data.sand_frac(catchments)./100,...
    CAMELS_signatures.BFI5(catchments),...
    'x_name','f{clay} [-]','y_name','f{sand} [-]','z_name','BFI5 [-]',...
    'ID',CAMELS_data.gauge_id,...
    'colour_scheme','YlGNBu','flip_colour_scheme',false,'nr_colours',10,...
    'x_limits',[0 .5],'y_limits',[0 1],...
    'z_limits',[0 1],'z_lower_limit_open',false,'z_upper_limit_open',false,...
    'figure_title','','figure_name','soil_texture',...
    'save_figure',false,'figure_path',fig_path,'figure_type','-dmeta')
disp(corr(CAMELS_signatures.BFI5(catchments),CAMELS_data.clay_frac(catchments),'Type','Spearman','Rows','Complete'))
disp(corr(CAMELS_signatures.BFI5(catchments),CAMELS_data.sand_frac(catchments),'Type','Spearman','Rows','Complete'))

%% Oregon age plots
catchments = [644:649 651:652 654:657]; % Oregon central cascades westwards draining

plotTrivariate(mean_age(catchments),CAMELS_signatures.BFI90(catchments),CAMELS_data.frac_snow(catchments),...
    'x_name','Mean age [Ma]','y_name','BFI{90} [-]','z_name','f{snow} [-]','ID',CAMELS_data.gauge_id(catchments),...
    'colour_scheme','Greys','flip_colour_scheme',false,'nr_colours',10,...
    'y_limits',[0.0 .5],'x_limits',[0 20],...
    'z_limits',[0.0 .5],'z_lower_limit_open',false,'z_upper_limit_open',false,...
    'figure_title','','figure_name','Mean_Age',...
    'save_figure',false,'figure_path',fig_path,'figure_type','-dmeta')
% set(gca,'YScale','log');
saveFig(gcf,strcat('trivariate_','Mean_Age'),fig_path)
disp(corr(mean_age(catchments),CAMELS_signatures.BFI90(catchments),'Type','Spearman'))

%% Ozarks sinkholes plot
sinkholes(Pre_Wisconsin_fraction>0) = NaN; % remove catchments overlain by glacial deposits

plotBivariate(sinkholes./CAMELS_data.area_gages2,CAMELS_signatures.BFI5,...
    'x_name','# Sinkholes / km2 ','y_name','BFI{5} [-]','ID',CAMELS_data.gauge_id,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'y_limits',[0 1],'x_limits',[0 0.5],...
    'figure_title','','figure_name','Sinkholes',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dmeta')
disp(corr(sinkholes./CAMELS_data.area_gages2,CAMELS_signatures.BFI5,'Type','Spearman','rows','complete'))

%% Edwards Trinity
catchments = [470 472:473 475:479 466:467 460:462 456:457 455 ]'; %459

BFI90 = CAMELS_signatures.BFI90(catchments);
% BFI90(BFI90==0) = 1e-4; % to be able to show data on log plot
plotTrivariate(edwards_trinity_frac(catchments),BFI90,...
    CAMELS_data.runoff_ratio(catchments),...
    'x_name','Edwards-Trinity frac. [-]','y_name','BFI{90} [-]','z_name','Q/P [-]','ID',CAMELS_data.gauge_id(catchments),...
    'colour_scheme','Oranges','flip_colour_scheme',true,'nr_colours',10,...
    'y_limit',[0 0.5],'x_limit',[0 1],...
    'z_limits',[0.0 .25],'z_lower_limit_open',false,'z_upper_limit_open',false,...
    'figure_title','','figure_name','Edwards',...
    'save_figure',false,'figure_path',fig_path,'figure_type','-dmeta')
% set(gca,'YScale','log');
% set(gca,'XDir','reverse')
% yticks([1e-4 1e-3 1e-2 1e-1 1])
% yticklabels({'0','10^{-3}','10^{-2}','10^{-1}','10^{-0}'})
% annotation(gcf,'line',[0.1475 0.1475],[0.23 0.263],'color',[1 1 1],'linewidth',3);
% annotation(gcf,'line',[0.135 0.157],[0.253 0.273]);
% annotation(gcf,'line',[0.135 0.157],[0.22 0.24]);
saveFig(gcf,strcat('trivariate_','Edwards'),fig_path)
disp(corr(edwards_trinity_frac(catchments),BFI90,'Type','Spearman'))

%% Lakes and wetlands
colour_mat = brewermap(12,'Paired');
% correlation for all catchments
disp(corr(CAMELS_signatures.BFI5,CAMELS_signatures.recession_beta,'Type','Spearman','Rows','Complete'))

% subsurface
plotBivariate(-1,-1,...
    'x_name','BFI{5} [-]','y_name','\betam [-]','ID',1,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'x_limits',[0 1],'y_limits',[0 7],...
    'figure_title','(c) Subsurface','figure_name','release','save_figure',false,'figure_path',fig_path)
p0=scatter(CAMELS_signatures.BFI5,CAMELS_signatures.recession_beta,20,'filled','MarkerFaceColor',[0.8 0.8 0.8]); %colour_mat(7,:)
p1=scatter(CAMELS_signatures.BFI5(~lakes&~snow),CAMELS_signatures.recession_beta(~lakes&~snow),20,'filled','MarkerFaceColor',[0.4 0.4 0.4]); %colour_mat(12,:)
% p2=scatter(signatures_CAMELS.BFI5(upland&~lakes&~snow&strcmp(data_CAMELS.geol_1st_class,'Siliciclastic sedimentary rocks')),signatures_CAMELS.recession_beta(upland&~lakes&~snow&strcmp(data_CAMELS.geol_1st_class,'Siliciclastic sedimentary rocks')),20,'filled','MarkerFaceColor',colour_mat(7,:)); %colour_mat(12,:)
% p3=scatter(signatures_CAMELS.BFI5(upland&~lakes&~snow&strcmp(data_CAMELS.geol_1st_class,'Metamorphics')),signatures_CAMELS.recession_beta(upland&~lakes&~snow&strcmp(data_CAMELS.geol_1st_class,'Metamorphics')),20,'filled','MarkerFaceColor',colour_mat(4,:)); %colour_mat(12,:)
% legend([p3 p2 p1],{'Metamorphic rock','Siliciclastic sedimentary rock','Other'},'box','off','location','northwest','FontSize',7,'ItemTokenSize',[7,7])
% set(gca,'YScale','log')
saveFig(gcf,strcat('bivariate_','recession_subsurface'),fig_path,'-dmeta')
disp(corr(CAMELS_signatures.BFI5(~lakes&~snow),CAMELS_signatures.recession_beta(~lakes&~snow),'Type','Spearman','Rows','Complete'))

% lakes
catchments = [288:294 297 363:365 367 368 145:158 160:162]; % Prairie Pothole Region and Florida
plotBivariate(-1,-1,...
    'x_name','BFI{5} [-]','y_name','\betam [-]','ID',1,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'x_limits',[0 1],'y_limits',[0 7],...
    'figure_title','(a) Surface Water','figure_name','release','save_figure',false,'figure_path',fig_path)
p0=scatter(CAMELS_signatures.BFI5,CAMELS_signatures.recession_beta,20,'filled','MarkerFaceColor',[0.8 0.8 0.8]); %colour_mat(7,:)
p1=scatter(CAMELS_signatures.BFI5(lakes),CAMELS_signatures.recession_beta(lakes),20,'filled','MarkerFaceColor',colour_mat(2,:));
p2=scatter(CAMELS_signatures.BFI5(catchments),CAMELS_signatures.recession_beta(catchments),25,'MarkerEdgeColor',[0 0 0],'linewidth',1.0);
% set(gca,'YScale','log')
saveFig(gcf,strcat('bivariate_','recession_lakes'),fig_path,'-dmeta')
disp(corr(CAMELS_signatures.BFI5(lakes),CAMELS_signatures.recession_beta(lakes),'Type','Spearman','Rows','Complete'))

% snow
plotBivariate(-1,-1,...
    'x_name','BFI{5} [-]','y_name','\betam [-]','ID',1,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'x_limits',[0 1],'y_limits',[0 7],...
    'figure_title','(b) Snow','figure_name','release','save_figure',false,'figure_path',fig_path)
p0=scatter(CAMELS_signatures.BFI5,CAMELS_signatures.recession_beta,20,'filled','MarkerFaceColor',[0.8 0.8 0.8]); %colour_mat(7,:)
p1=scatter(CAMELS_signatures.BFI5(snow&~lakes),CAMELS_signatures.recession_beta(snow&~lakes),20,'filled','MarkerFaceColor',colour_mat(10,:));%colour_mat(10,:));
% set(gca,'YScale','log')
saveFig(gcf,strcat('bivariate_','recession_snow'),fig_path,'-dmeta')
disp(corr(CAMELS_signatures.BFI5(snow&~lakes),CAMELS_signatures.recession_beta(snow&~lakes),'Type','Spearman','Rows','Complete'))

%% GLHYMPS (to check catchment averages)
plotBivariate(porosity,CAMELS_data.geol_porosity,...
    'x_name','Porosity [-]','y_name','Porosity CAMELS [-]','ID',CAMELS_data.gauge_id,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'figure_title','(a)','figure_name','porosity_comparison','save_figure',true,'figure_path',fig_path)

plotMapUS(CAMELS_data.gauge_lat,CAMELS_data.gauge_lon,porosity,...
    'attribute_name','Porosity [-]','ID',CAMELS_data.gauge_id,...
    'colour_scheme','PRGn','flip_colour_scheme',false,...
    'c_limits',[0 .25],'c_lower_limit_open',false,'c_upper_limit_open',false,...
    'figure_title','','figure_name','Porosity')

plotBivariate(log10(permeability),CAMELS_data.geol_permeability,...
    'x_name','Permeability [m^2, log10]','y_name','Permeability CAMELS [m^2, log10]','ID',CAMELS_data.gauge_id,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'figure_title','(b)','figure_name','permeability_comparison','save_figure',true,'figure_path',fig_path)

plotMapUS(CAMELS_data.gauge_lat,CAMELS_data.gauge_lon,log10(permeability),...
    'attribute_name','Permeability [m^2, log10]','ID',CAMELS_data.gauge_id,...
    'colour_scheme','RdYlBu','flip_colour_scheme',false,...
    'c_limits',[-17 -11],'c_lower_limit_open',false,'c_upper_limit_open',false,...
    'figure_title','','figure_name','Permeability')

%% maps of baseflow signatures
plotMapUS(CAMELS_data.gauge_lat,CAMELS_data.gauge_lon,CAMELS_signatures.BFI5,...
    'attribute_name','BFI5 [-]','ID',CAMELS_data.gauge_id,...
    'colour_scheme','YlGnBu','flip_colour_scheme',false,...
    'c_limits',[0 1],'c_lower_limit_open',false,'c_upper_limit_open',false,...
    'figure_title','(a)','figure_name','BFI5','save_figure',true,'figure_path',fig_path)

plotMapUS(CAMELS_data.gauge_lat,CAMELS_data.gauge_lon,CAMELS_signatures.BFI90,...
    'attribute_name','BFI{90} [-]','ID',CAMELS_data.gauge_id,...
    'colour_scheme','RdPu','flip_colour_scheme',false,...
    'c_limits',[0 .5],'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','(b)','figure_name','BFI90','save_figure',true,'figure_path',fig_path)

plotMapUS(CAMELS_data.gauge_lat,CAMELS_data.gauge_lon,CAMELS_signatures.recession_beta,...
    'attribute_name','\betam [-]','ID',CAMELS_data.gauge_id,...
    'colour_scheme','PuOr','flip_colour_scheme',false,...
    'c_limits',[0 4],'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','(c)','figure_name','beta','save_figure',true,'figure_path',fig_path)

%% Correlation plots
plotBivariate(CAMELS_signatures.BFI5,CAMELS_signatures.BFI90,...
    'x_name','BFI5 [-]','y_name','BFI90 [-]','ID',CAMELS_data.gauge_id,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'figure_title','(a)','figure_name','BFI5_BFI90',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dmeta')
disp(corr(CAMELS_signatures.BFI5,CAMELS_signatures.BFI90,'Type','Spearman','Rows','Complete'))

plotBivariate(CAMELS_signatures.BFI5,CAMELS_signatures.recession_beta,...
    'x_name','BFI5 [-]','y_name','\betam [-]','ID',CAMELS_data.gauge_id,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'figure_title','(b)','figure_name','BFI5_betam',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dmeta')
disp(corr(CAMELS_signatures.BFI5,CAMELS_signatures.recession_beta,'Type','Spearman','Rows','Complete'))

plotBivariate(CAMELS_signatures.BFI90,CAMELS_signatures.recession_beta,...
    'x_name','BFI90 [-]','y_name','\betam [-]','ID',CAMELS_data.gauge_id,...
    'show_corr',false,'show_fit',false,'show_hist',false,...
    'figure_title','(c)','figure_name','BFI90_betam',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dmeta')
disp(corr(CAMELS_signatures.BFI90,CAMELS_signatures.recession_beta,'Type','Spearman','Rows','Complete'))

%% print correlations between attributes and BFI etc.
% catchments = [1:671]';
additional_attributes = [
    mean_age, freshwater_fraction+lake_fraction, ...
    sinkholes, edwards_trinity_frac]; % matrix with additional attributes
dispCorrelationsCAMELS(CAMELS_data,CAMELS_signatures.BFI5,'BFI5',...
    catchments,additional_attributes)
dispCorrelationsCAMELS(CAMELS_data,CAMELS_signatures.BFI90,'BFI90',...
    catchments,additional_attributes)
dispCorrelationsCAMELS(CAMELS_data,CAMELS_signatures.recession_beta,'\betam',...
    catchments,additional_attributes)
