%% calculate_signatures
%
%   Calculates baseflow signatures for CAMELS catchments and saves them in 
%   a mat-file.
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

%% load catchments
if exist('data_CAMELS_struc') % check if data is already loaded
else
    data_CAMELS_struc = load('./CAMELS_Matlab/Data/CAMELS_data.mat');
end
CAMELS_data = data_CAMELS_struc.CAMELS_data;
clear data_CAMELS_struc

%% loop over all catchments
n_CAMELS = length(CAMELS_data.gauge_id);
BFI5 = NaN(n_CAMELS,1); % UKIH method
BFI30 = NaN(n_CAMELS,1); % UKIH method
BFI60 = NaN(n_CAMELS,1); % UKIH method
BFI90 = NaN(n_CAMELS,1); % UKIH method
BFI120 = NaN(n_CAMELS,1); % UKIH method
BFI180 = NaN(n_CAMELS,1); % UKIH method
BFI_LH = NaN(n_CAMELS,1); % Lyne-Hollick filter
Q5 = NaN(n_CAMELS,1); %
Qmean = NaN(n_CAMELS,1); %
recession_beta = NaN(n_CAMELS,1); % recession analysis using Exponential Time Stepping dQdt calculation
recession_alpha = NaN(n_CAMELS,1);
recession_beta_spring = NaN(n_CAMELS,1); % seasonal recession analysis
recession_alpha_spring = NaN(n_CAMELS,1);
recession_beta_summer = NaN(n_CAMELS,1);
recession_alpha_summer = NaN(n_CAMELS,1);
recession_beta_autumn = NaN(n_CAMELS,1);
recession_alpha_autumn = NaN(n_CAMELS,1);
recession_beta_winter = NaN(n_CAMELS,1);
recession_alpha_winter = NaN(n_CAMELS,1);
recession_beta_BN = NaN(n_CAMELS,1); % recession analysis using Brutsaert-Nieber dQdt calculation
recession_alpha_BN = NaN(n_CAMELS,1);

for i = 1:n_CAMELS
    
    if mod(i,10) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,n_CAMELS)
    end
    
    % extract data
    t = datetime(CAMELS_data.P{i}(:,1),'ConvertFrom','datenum');
    P = CAMELS_data.P{i}(:,2);
    PET = CAMELS_data.PET{i}(:,2);
    Q = CAMELS_data.Q{i}(:,2);
    ID = CAMELS_data.gauge_id(i);
    
    % calculate signatures
    BFI5(i) = sig_BFI(Q,t,'method','UKIH');
    BFI30(i) = sig_BFI(Q,t,'method','UKIH','parameters',30);
    BFI60(i) = sig_BFI(Q,t,'method','UKIH','parameters',60);
    BFI90(i) = sig_BFI(Q,t,'method','UKIH','parameters',90);
    BFI120(i) = sig_BFI(Q,t,'method','UKIH','parameters',120);
    BFI180(i) = sig_BFI(Q,t,'method','UKIH','parameters',180);
    BFI_LH(i) = sig_BFI(Q,t,'method','Lyne_Hollick');
    
    Q5(i) = sig_x_percentile(Q,t,5);
    Qmean(i) = sig_Q_mean(Q,t);
    
    [para_mat, recession_month] = sig_RecessionAnalysis(Q,t,'plot_results',false);
    [~, ind] = min(abs(para_mat(:,2) - median(para_mat(:,2)))); % find recession according to median exponent
    recession_alpha(i) = para_mat(ind,1);
    recession_beta(i) = para_mat(ind,2);
    
    % seasonal recessions
    spring = (recession_month==4|recession_month==5|recession_month==6);
    if sum(spring) < 1
    else
        para_mat_spring = para_mat(spring,:);
        [~, ind] = min(abs(para_mat_spring(:,2) - median(para_mat_spring(:,2))));
        recession_alpha_spring(i) = para_mat_spring(ind,1);
        recession_beta_spring(i) = para_mat_spring(ind,2);
    end
    
    summer = (recession_month==7|recession_month==8|recession_month==9);
    if sum(summer) < 1
    else
        para_mat_summer = para_mat(summer,:);
        [~, ind] = min(abs(para_mat_summer(:,2) - median(para_mat_summer(:,2))));
        recession_alpha_summer(i) = para_mat_summer(ind,1);
        recession_beta_summer(i) = para_mat_summer(ind,2);
    end
    
    autumn = (recession_month==10|recession_month==11|recession_month==12);
    if sum(autumn) < 1
    else
        para_mat_autumn = para_mat(autumn,:);
        [~, ind] = min(abs(para_mat_autumn(:,2) - median(para_mat_autumn(:,2))));
        recession_alpha_autumn(i) = para_mat_autumn(ind,1);
        recession_beta_autumn(i) = para_mat_autumn(ind,2);
    end
    
    winter = (recession_month==1|recession_month==2|recession_month==3);
    if sum(winter) < 1
    else
        para_mat_winter = para_mat(winter,:);
        [~, ind] = min(abs(para_mat_winter(:,2) - median(para_mat_winter(:,2))));
        recession_alpha_winter(i) = para_mat_winter(ind,1);
        recession_beta_winter(i) = para_mat_winter(ind,2);
    end
    
    [para_mat] = sig_RecessionAnalysis(Q,t,'dQdt_method','BN');
    [~, ind] = min(abs(para_mat(:,2) - median(para_mat(:,2)))); % find recession according to median exponent
    recession_alpha_BN(i) = para_mat(ind,1);
    recession_beta_BN(i) = para_mat(ind,2);
    
end

CAMELS_signatures.BFI5 = BFI5;
CAMELS_signatures.BFI30 = BFI30;
CAMELS_signatures.BFI60 = BFI60;
CAMELS_signatures.BFI90 = BFI90;
CAMELS_signatures.BFI120 = BFI120;
CAMELS_signatures.BFI180 = BFI180;
CAMELS_signatures.BFI_LH = BFI_LH;
CAMELS_signatures.Qmean = Qmean;
CAMELS_signatures.Q5 = Q5;
CAMELS_signatures.Q5n = Q5./Qmean;
CAMELS_signatures.recession_beta = recession_beta;
CAMELS_signatures.recession_alpha = recession_alpha;
CAMELS_signatures.recession_beta_spring = recession_beta_spring;
CAMELS_signatures.recession_alpha_spring = recession_alpha_spring;
CAMELS_signatures.recession_beta_summer = recession_beta_summer;
CAMELS_signatures.recession_alpha_summer = recession_alpha_summer;
CAMELS_signatures.recession_beta_autumn = recession_beta_autumn;
CAMELS_signatures.recession_alpha_autumn = recession_alpha_autumn;
CAMELS_signatures.recession_beta_winter = recession_beta_winter;
CAMELS_signatures.recession_alpha_winter = recession_alpha_winter;
CAMELS_signatures.recession_beta_BN = recession_beta_BN;
CAMELS_signatures.recession_alpha_BN = recession_alpha_BN;

save('./Baseflow_signatures/Results/CAMELS_signatures.mat','CAMELS_signatures')
