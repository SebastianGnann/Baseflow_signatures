function [] = plotDFI(data_CAMELS,i,save_figure)
%plotDFI Plots DFI for various parameter values.
% Optional:
%   - ...
%
% INPUT
%   catchment_data: structure with CAMELS data
%   i: Matlab index of catchment
%   save_figure: save plot true/false
%
%   OUTPUT
%   plot
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

% check input parameters
if nargin < 3
    error('Not enough input arguments.')
end

%% extract data
t = datetime(data_CAMELS.P{i}(:,1),'ConvertFrom','datenum');
P = data_CAMELS.P{i}(:,2);
PET = data_CAMELS.PET{i}(:,2);
Q = data_CAMELS.Q{i}(:,2);
ID = data_CAMELS.gauge_id(i);

%% estimate baseflow
% UKIH method
max_N = 365;
day_blocks = 1:max_N;
DFI = NaN(size(day_blocks));

for j = 1:length(day_blocks)
    DFI(j) = sig_BFI(Q,t,'method','UKIH','parameters',[day_blocks(j)]);
end

%% plot
colour_mat = brewermap(6,'PuBuGn');
colour_mat(1,:) = [];
hold on
plot(day_blocks(1:5),DFI(1:5),'linewidth',2,'color',colour_mat(4,:))
plot(day_blocks(5:90),DFI(5:90),'linewidth',2,'color',colour_mat(3,:))
plot(day_blocks(90:end),DFI(90:end),'linewidth',2,'color',colour_mat(1,:))
leg = legend({'<=5d','<=90d','<=365d'},'location','best','box','off');
leg.FontSize = 7;
xlim([0 100])
ylim([0 1])
xlabel('filter parameter [d]')
ylabel('DFI')
% set(gca,'XScale','log')
% set(gca,'YScale','log')

v = get(gca,'Position');
set(gca,'Position',[v(1) v(2)*2 v(3) v(4)*0.8])

%% save fig
if save_figure
    fig = gcf;
    set(fig,'Units','Inches');
    position = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[position(3),position(4)]);
    fig_name = strcat('DFI_',figure_name);
    path_name = './Baseflow_and_geology/Images';
    path = strcat(path_name,'\',fig_name);
    print(fig,path,'-dpdf','-r500');
end

end
