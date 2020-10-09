function [] = plotHydrograph(data_CAMELS,i,save_figure,title_str,max_flow)
%plotHydrograph Plots hydrograph with estimated baseflow components.
% Optional:
%   - ...
%
% INPUT
%   catchment_data: structure with CAMELS data
%   i: Matlab index of catchment
%   save_figure: save plot true/false
%   title_str: title of plot
%   max_flow: maximum flow for plot (y-axis)
%
%   OUTPUT
%   plot
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

% check input parameters
if nargin < 5
    error('Not enough input arguments.')
end

%% extract data
t = datetime(data_CAMELS.P{i}(:,1),'ConvertFrom','datenum');
P = data_CAMELS.P{i}(:,2);
PET = data_CAMELS.PET{i}(:,2);
Q = data_CAMELS.Q{i}(:,2);
ID = data_CAMELS.gauge_id(i);

if nargin < 4
    max_flow = max(Q);
end

%% estimate baseflow
% UKIH method
B5 = util_UKIH_Method(Q,'n_days',5);
B5(isnan(B5)) = 0;
B90 = util_UKIH_Method(Q,'n_days',90);
B90(isnan(B90)) = 0;

%% plot
colour_mat = brewermap(6,'PuBuGn');
colour_mat(1,:) = [];

hold on
% yyaxis left
% plot(t,Q,'-','color',[.7 .7 .7],'Linewidth',0.5)
try
    Q_tmp = Q; Q_tmp(isnan(Q)) = nanmedian(Q);
    B5(isnan(B5)) = nanmedian(B5);
    B90(isnan(B90)) = nanmedian(B90);
    t2 = [t', fliplr(t')];
    p1 = fill(t2,[B5', fliplr(Q_tmp')],'-','FaceColor',colour_mat(5,:),'EdgeColor','none');
    p2 = fill(t2,[B90', fliplr(B5')],'-','FaceColor',colour_mat(3,:),'EdgeColor','none');
    p3 = fill(t2,[zeros(size(B90))', fliplr(B90')],'-','FaceColor',colour_mat(1,:),'EdgeColor','none');
catch
    %     disp('Problem')
end
set(gca,'ycolor','k')
ylabel('Flow [mm/d]')
% set(gca,'YScale','log')

ylim([0 max_flow])

% yyaxis right
% p4 = plot(t,movmean(P,1),'color',0.4.*[1 1 1],'Linewidth',0.5);
% p4.Color(4) = 0.3;
% p5 = plot(t,movmean(PET,1),'-','color',0.4.*[1 1 1],'Linewidth',0.5);
% p5.Color(4) = 0.6;
% % plot(t,movmean(-PET,1),'-','color',0.6.*[1 1 1],'Linewidth',1.0)
% set(gca,'Ydir','reverse')
% set(gca,'ycolor','k')
% ca = gca; yright = ca.YAxis(2);
% yright.Visible = 'off';

% ticks = datenum(1980:0.5:2010,1,1);
% set(gca, 'xtick', ticks);
% datetick('x', 'mm-yy', 'keepticks');
% legend_P = strcat('P (fs=',num2str(round(data_CAMELS.frac_snow(i),2)),')');
% legend_P = strcat('P');

xlabel('Date')
xlim([datetime(1999,10,1) datetime(2002,9,30)])
% datetick('x','mmm-YY') %,'keepticks'
set(gca,'Layer','top')

% leg = legend([p1 p2 p3 p4 p5],{'Q','Q{b,5}','Q{b,90}',legend_P,'PET'},...
%     'Position',[0.779 0.520 0.122 0.310],'box','off');
leg = legend([p1 p2 p3],{'Q','Q{b,5}','Q{b,90}'},...
    'Position',[0.779 0.716 0.137 0.228],'box','off');

leg.FontSize = 7;
leg.ItemTokenSize = [7,7];

v = get(gca,'Position');
% set(gca,'Position',[v(1).*0.7 v(2)*2 v(3) v(4)*0.8])
set(gca,'Position',[v(1) v(2) v(3)*0.8 v(4)])

% dim =  [0.286 0.9 0.142 0.128];
% annotation('textbox',dim,'String',title_str,'FitBoxToText','on',...
%     'FontSize',10,'FontWeight','bold','EdgeColor','None');

%% save fig
if save_figure
    fig = gcf;
    set(fig,'Units','Inches');
    position = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[position(3),position(4)]);
    fig_name = strcat('hydrograph_',figure_name);
    path_name = './Baseflow_and_geology/Images';
    path = strcat(path_name,'\',fig_name);
    print(fig,path,'-dpdf','-r0');
end

end
