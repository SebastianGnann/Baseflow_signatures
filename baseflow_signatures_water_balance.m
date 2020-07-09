%% baseflow_signatures_water_balance
%
%   Estimates regional groundwater flow via the water balance 
%   (see Supplement).
%
% ---
%
% Sebastian Gnann (2020)
% sebastian.gnann@bristol.ac.uk

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

%% load catchments
if exist('data_CAMELS_struc') % check if data is already loaded
else
    data_CAMELS_struc = load('./CAMELS_Matlab/Data/CAMELS_data.mat');
end
CAMELS_data = data_CAMELS_struc.CAMELS_data;
clear data_CAMELS_struc

%% calculate mean annual data
% get mean from 2000 to 2013 (last year in CAMELS)
Q_m = NaN(size(CAMELS_data.gauge_id));
P_m = NaN(size(CAMELS_data.gauge_id));
PET_m = NaN(size(CAMELS_data.gauge_id));
n_years = NaN(size(CAMELS_data.gauge_id));

for i = 1:length(CAMELS_data.gauge_id)
    
    t = datetime(CAMELS_data.P{i}(:,1),'ConvertFrom','datenum'); 
    P = CAMELS_data.P{i}(:,2);
    PET = CAMELS_data.PET{i}(:,2);
    Q = CAMELS_data.Q{i}(:,2);

    [P_tmp, ~, year_P] = aggregateTimeSeries([datenum(t), P], 1);
    [PET_tmp, ~, year_PET] = aggregateTimeSeries([datenum(t), PET], 1);
    [Q_tmp, ~, year_Q] = aggregateTimeSeries([datenum(t), Q], 1);

    P_m(i) = nanmean(P_tmp(year_P>1999 & year_P<2014)).*365;
    PET_m(i) = nanmean(PET_tmp(year_PET>1999 & year_PET<2014)).*365;
    Q_m(i) = nanmean(Q_tmp(year_Q>1999 & year_Q<2014)).*365;
    n_years = length(year_P(year_P>1999 & year_P<2014)); % check
    
end

%% plots
ET_modis = load('modis.mat'); ET_modis = ET_modis.data.attribute_mean;
ET_gleam = load('gleam.mat'); ET_gleam = ET_gleam.data.attribute_mean;
ET_Budyko = P_m.*BudykoCurve(PET_m./P_m);
ET_water_balance = (1-Q_m./P_m).*P_m;

% scatter plots
i = 0;
label = ["(a)","(b)","(c)"];
for ET_str = ["ET_modis", "ET_gleam", "ET_Budyko"]
    i = i+1;
    eval(strcat('ET = ',ET_str,';'));    
    plotTrivariate(ET./P_m,Q_m./P_m,(ET+Q_m)./P_m,...
        'x_name','E_a/P','y_name','Q/P','z_name','(E_a+Q)/P','ID',CAMELS_data.gauge_id,...
        'colour_scheme','RdBu','flip_colour_scheme',false,...
        'x_limits',[0 1],'y_limits',[0 1],...
        'z_limits',[0.75 1.25],'z_lower_limit_open',true,'z_upper_limit_open',true,...
        'figure_title',char(label(i)),'figure_name',char(ET_str),'save_figure',true,'figure_path',fig_path)
    hold on; plot(0:.1:1,1:-.1:0,'k');
end

% maps
i = 0;
for ET_str = ["ET_modis", "ET_gleam", "ET_Budyko"]
    i = i+1;
    eval(strcat('ET = ',ET_str,';'));
    wb = (P_m - (ET+Q_m))./P_m;
    plotMapUS([CAMELS_data.gauge_lat]',[CAMELS_data.gauge_lon]',wb,...
        'attribute_name','Q_{gw}/P','ID',[CAMELS_data.gauge_id]',...
        'colour_scheme','RdBu','flip_colour_scheme',true,'nr_colours',7,...
        'c_limits',[-0.25 0.25],'c_lower_limit_open',true,'c_upper_limit_open',true,...
        'figure_title',char(label(i)),'figure_name',char(ET_str),...
        'save_figure',true,'figure_path',fig_path,'figure_type','-dmeta')
%     set(gca,'ColorScale','log');
%     saveFig(gcf,['map_US_',char(ET_str)],fig_path)
end

% Budyko and related plots
i = 0;
for ET_str = ["ET_modis", "ET_gleam", "ET_water_balance"]
    i = i+1;
    eval(strcat('ET = ',ET_str,';')); 
    disp(corr(PET_m./P_m,ET,'type','Spearman'))
    plotTrivariate(PET_m./P_m,ET,CAMELS_data.frac_snow,...
        'x_name','PET/P','y_name','E_a [mm/year]','z_name','frac snow','ID',CAMELS_data.gauge_id,...
        'colour_scheme','RdBu','flip_colour_scheme',false,...
        'x_limits',[0.0 4],'y_limits',[0 4*365],...
        'z_limits',[0 1],'z_lower_limit_open',false,'z_upper_limit_open',false,...
        'figure_title',char(label(i)),'figure_name',char(ET_str),'save_figure',false,'figure_path',fig_path)
%     hold on; plot(0:.1:1,0:.1:1,'k');
%     plot(1:.1:5,ones(size(1:.1:5)),'k');
%     plot(PET_m./P_m,ET_Budyko./P_m,'.k');
%     set(gca,'xScale','log')
end
