function [] = plotAnnualFlows(data_CAMELS,i,save_figure)
%plotAnnualFlows Plots annual base flow and quick flow vs annual precip.
%
%   INPUT
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

Q_b = util_UKIH_Method(Q, 'n_days', 5);
Q_f = Q - Q_b;

% plot annual baseflow against annual P
[P_annual, ~, ~] = aggregateTimeSeries([datenum(t), P], true);
% [PET_annual, ~, ~] = aggregateTimeSeries([datenum(t), PET], true);
% [Q_annual, ~, ~] = aggregateTimeSeries([datenum(t), Q], true);
[Q_b_annual, ~, ~] = aggregateTimeSeries([datenum(t), Q_b], true);
[Q_f_annual, ~, ~] = aggregateTimeSeries([datenum(t), Q_f], true);

hold on;
isn = isnan(Q_b_annual);
P_annual(isn) = [];
Q_b_annual(isn) = [];
Q_f_annual(isn) = [];
P_annual = P_annual.*365;
Q_b_annual = Q_b_annual.*365;
Q_f_annual = Q_f_annual.*365;

plot(P_annual,Q_b_annual,'b .')
plot(P_annual,Q_f_annual,'r .')
xlabel('Annual P [mm]');
ylabel('Annual Q_b and Q_f [mm]');
% xlim([0 ceil(1.1*max(P_a))]); ylim([0 1.1*max(Q_b_a)])
n = length(P_annual);
SSxy = sum(P_annual.*Q_b_annual) - sum(P_annual)*sum(Q_b_annual)/n;
SSxx = sum(P_annual.^2) - sum(P_annual)^2/n;
b = SSxy/SSxx;
a = mean(Q_b_annual) - b*mean(P_annual);
SSxy = sum(P_annual.*Q_f_annual) - sum(P_annual)*sum(Q_f_annual)/n;
SSxx = sum(P_annual.^2) - sum(P_annual)^2/n;
b2 = SSxy/SSxx;
a2 = mean(Q_f_annual) - b2*mean(P_annual);
% plot([0:ceil(1.1*max(P_a))],(a+b.*[0:ceil(1.1*max(P_a))]),'k-')
plot([floor(0.9*min(P_annual)):ceil(1.1*max(P_annual))],(a+b.*[floor(0.9*min(P_annual)):ceil(1.1*max(P_annual))]),'b-')
plot([floor(0.9*min(P_annual)):ceil(1.1*max(P_annual))],(a2+b2.*[floor(0.9*min(P_annual)):ceil(1.1*max(P_annual))]),'r-')
%     leg = legend([p1; p2],{str1; str2},'location','best','box','off');
%     leg.FontSize = 7;
str1 = (sprintf('Q_b=%.0f+%.2f P \n',a,b));
str2 = (sprintf('Q_f=%.0f+%.2f P \n',a2,b2));

% dim = [0.525 0.847 0.108 0.122];
% annotation('textbox',dim,'String',[str1,str2],...
%     'FitBoxToText','on','FontSize',7,'EdgeColor','None');

leg = legend({str1,str2},'location','best','box','off');
leg.FontSize = 7;

% Q_min_monthly(isn) = [];
% Q_max_monthly(isn) = [];
% plot(P_a,Q_min_monthly.*365.25,'.','color',[0.5 0.5 0.5])
% plot(P_a,Q_max_monthly.*365.25,'.','color',[0.2 0.2 0.2])

v = get(gca,'Position');
set(gca,'Position',[v(1) v(2)*2 v(3) v(4)*0.8])

colour_mat = brewermap(12,'Spectral');
colour_mat(6:7,:) = [];

%% save fig
if save_figure
    fig = gcf;
    set(fig,'Units','Inches');
    position = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[position(3),position(4)]);
    fig_name = strcat('annual_flows_',figure_name);
    path_name = './Baseflow_and_geology/Images';
    path = strcat(path_name,'\',fig_name);
    print(fig,path,'-dpdf','-r500');
end

end
