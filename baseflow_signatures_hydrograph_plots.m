%% baseflow_signatures_hydrograph_plots
%
%   Loads data and creates hydrographs shown in the paper.
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk( 2020)

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

%% catchment analysis
% loop over regions

% additional topography data (used for perceptual models)
elevation = load('elevation.mat');
relief = elevation.data.relief;

region_list = ["EX","NC","GL","LA","OR","OZ","TX"];
clc
for k = 1:length(region_list)
    region = region_list(k);
    switch region
        case 'EX'
            catchments = [250]; Qmax = 3; xl = [3*10^-1 10^2/5]; yl = [2*10^-3 10^2/2]; % Example
        case 'GL'
            catchments = [318 267]; Qmax = 2.5; xl = [10^-2 10^1]; yl = [10^-4 10^2]; % Glacial Areas
        case 'NC'
            catchments = [123 128]; Qmax = 1.5; xl = [10^-2 10^1]; yl = [10^-4 10^2]; % North Carolina
        case 'OR'
            catchments = [656 644]; Qmax = 15; xl = [10^-1 10^2]; yl = [10^-2 10^2]; % Oregon
        case 'OZ'
            catchments = [394 404]; Qmax = 2.0; xl = [10^-1 10^1]; yl = [10^-3 10^2]; % Ozarks
        case 'TX'
            catchments = [475 461]; Qmax = 1.0; xl = [10^-4 10^2]; yl = [10^-4 10^1];  % Texas
        case 'TX2'
            catchments = [474 479]; Qmax = 1.0; xl = [10^-4 10^2]; yl = [10^-4 10^1];  % Texas
        case 'LA'
            catchments = [289 147]; Qmax = 2; xl = [10^-3 10^1]; yl = [10^-4 10^1]; % Lakes / North Dakota / Florida
    end
    
    counter = 0;
    for i = catchments %1:length(data_CAMELS.gauge_id)
%         counter = counter + 1;
        
        fig1 = figure('Name',num2str(i),'NumberTitle','off',...
            'pos',[100 100 350 180],...
            'Renderer','painters'); % OpenGL leads to strange plots
        
        % plot hydrograph
        t = datetime(CAMELS_data.P{i}(:,1),'ConvertFrom','datenum');
        Q = CAMELS_data.Q{i}(:,2);
        %         subplot(4,2,1:4)
        plotHydrograph(CAMELS_data,i,false,'',Qmax); %sig_x_percentile(Q,t,90)
        
        % plot recessions
        % subplot(4,2,[5 7])
        % plotRecessions(data_CAMELS,i,false,'',xl,yl); %
        
        % plot catchment attributes
        plotCatchmentAttributes(CAMELS_data,CAMELS_signatures,relief,i);
        
%         if counter == 1
%             title(strcat(num2str(data_CAMELS.gauge_id(i)))) 
%         else
%             title(strcat(num2str(data_CAMELS.gauge_id(i))))
%         end
        
        % save figure
        saveFig(fig1,strcat('hydrograph_',num2str(i)),fig_path,'-dmeta')
        
    end
    
end
