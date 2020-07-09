function [] = plotRecessions(data_CAMELS,i,save_figure,title_str,xl,yl)
%plotRecessions Plots recessions.
% Optional:
%   - ...
%
% INPUT
%   catchment_data: structure with CAMELS data
%   i: Matlab index of catchment
%   save_figure: save plot true/false
%   title_str: title of plot
%   xl: xlimits
%   yl: ylimits
%
%   OUTPUT
%   plot
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

% check input parameters
if nargin < 4
    error('Not enough input arguments.')
end

%% extract data
t = datetime(data_CAMELS.P{i}(:,1),'ConvertFrom','datenum');
P = data_CAMELS.P{i}(:,2);
PET = data_CAMELS.PET{i}(:,2);
Q = data_CAMELS.Q{i}(:,2);
ID = data_CAMELS.gauge_id(i);

%% recession analysis
[para_mat] = sig_RecessionAnalysis(Q, t, 'plot_results', true,...
    'fitIndividual',true,'fitting_type','linear');
% figure; plot(para_mat(:,1),para_mat(:,2),'o')

[~, ind] = min(abs(para_mat(:,2) - median(para_mat(:,2)))); % find recession according to median exponent
str = (sprintf('Median recession:\n-dQ/dt = %.3f Q^{%.1f}',para_mat(ind,1),para_mat(ind,2)));

v = get(gca,'Position');
% set(gca,'Position',[v(1) v(2)*2 v(3) v(4)*0.8])
% set(gca,'Position',[v(1)*1.05 v(2)*2 v(3) v(4)*0.8])
set(gca,'Position',[v(1) v(2) v(3) v(4).*0.8])
leg = get(gca,'Legend');
leg.Visible = 'off';
% set(leg,'Position',[leg.Position(1)*0.95 leg.Position(2)*1.05 leg.Position(3) leg.Position(4)])
% leg.ItemTokenSize = [0,0];
% dim =  [0.737 0.735 0.124 0.174];
dim =  [0.129505720823799 0.313727687150444 0.248283758097983 0.0972540062133999];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'FontSize',6,'EdgeColor','None');
% axis equal
if nargin < 5
else
    xlim(xl)
    ylim(yl)
end
xlabel('Q [mm/d]')
ylabel('-dQ/dt [mm/d^2]')

% dim =  [0.819 0.9 0.142 0.128];
% annotation('textbox',dim,'String',title_str,'FitBoxToText','on',...
%     'FontSize',10,'FontWeight','bold','EdgeColor','None');

% legend
colour_mat = brewermap(11,'Spectral');
p0=plot(0.00001:.1:100,0.00001:.1:100,'--','color',[.5 .5 .5]);
p1=plot(0,0,'.','Color',colour_mat(5,:)); 
p2=plot(0,0,'.','Color',colour_mat(2,:));  
p3=plot(0,0,'.','Color',colour_mat(9,:)); 
p4=plot(0,0,'.','Color',colour_mat(11,:)); 
% leg = legend([p1 p2 p3 p4],{'MAM','JJA','SON','DJF'},'box','off',...
%     'Position',[0.882 0.215 0.073 0.254]);
leg = legend([p1 p2 p3 p4],{'MAM','JJA','SON','DJF'},'box','off',...
    'Position',[0.367 0.103 0.113 0.144]);
leg.FontSize = 6;
leg.ItemTokenSize = [7,7];


%% save fig
if save_figure
    fig = gcf;
    set(fig,'Units','Inches');
    position = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[position(3),position(4)]);
    fig_name = strcat('catchment_data',figure_name);
    path_name = './Baseflow_and_geology/Images';
    path = strcat(path_name,'\',fig_name);
    print(fig,path,'-dpdf','-r500');
end

end
