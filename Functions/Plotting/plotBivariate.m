function [] = plotBivariate(x,y,varargin)
%PLOTBIVARIATE Plots scatter plot of two variables.
%   Optional:
%       - plot correlation
%       - plot fit
%       - plot histograms
%
%   INPUT
%   x: variable 1
%   y: variable 2
%   x_name: string with name of x
%   y_name: string with name of y
%   ID: if data cursor shall show ID input the respective array
%   x_limit: axis limits for x-axis
%   y_limit: axis limits for y-axis
%   show_corr: boolean to specify whether to show correlation
%   show_fit: boolean to specify whether to plot linear regression line
%   show_hist: boolean to specify whether to plot histogram
%   figure_title: title of plot, e.g. '(a)'
%   figure_name: name for saving, e.g. UK_P_Q
%   save_figure: save plot true/false
%   figure_path: path to folder where figure should be saved
%   figure_type: figure type, e.g. -dpdf or -dmeta
%
%   OUTPUT
%   plot
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;

addRequired(ip, 'xdata', ...
    @(x) isnumeric(x) && (size(x,1)==1 || size(x,2)==1))
addRequired(ip, 'ydata', ...
    @(y) isnumeric(y) && (size(y,1)==1 || size(y,2)==1))

addParameter(ip, 'x_name', [], @ischar)
addParameter(ip, 'y_name', [], @ischar)
addParameter(ip, 'ID', NaN(size(x)), @isnumeric)
addParameter(ip, 'x_limits', [min(x) max(x)], @(z) isnumeric(z) && length(z)==2)
addParameter(ip, 'y_limits', [min(y) max(y)], @(z) isnumeric(z) && length(z)==2)
addParameter(ip, 'show_corr', false, @islogical)
addParameter(ip, 'show_fit', false, @islogical)
addParameter(ip, 'show_hist', false, @islogical)
addParameter(ip, 'figure_title', '', @ischar)
addParameter(ip, 'figure_name', 'no_name', @ischar)
addParameter(ip, 'save_figure', false, @islogical)
addParameter(ip, 'figure_path', '', @ischar)
addParameter(ip, 'figure_type', '-dpdf', @ischar)

parse(ip, x, y, varargin{:})

x_name = ip.Results.x_name;
y_name = ip.Results.y_name;
ID = ip.Results.ID;
x_limits = ip.Results.x_limits;
y_limits = ip.Results.y_limits;
show_corr = ip.Results.show_corr;
show_fit = ip.Results.show_fit;
show_hist = ip.Results.show_hist;
figure_title = ip.Results.figure_title;
figure_name = ip.Results.figure_name;
save_figure = ip.Results.save_figure;
figure_path = ip.Results.figure_path;
figure_type = ip.Results.figure_type;

%% plotting
% get rid of NaN points
xNaN = isnan(x); x(xNaN) = []; y(xNaN) = [];
yNaN = isnan(y); x(yNaN) = []; y(yNaN) = [];
index = [1:length(ID)]';
ID(xNaN) = []; ID(yNaN) = [];
index(xNaN) = []; index(yNaN) = [];

% nice colours
colour_mat = brewermap(10,'YlGnBu');
if show_hist
    fig = figure('Name',figure_name,'NumberTitle','off','pos',[10 10 500 500]);
else
    fig = figure('Name',figure_name,'NumberTitle','off','pos',[10 10 250 230]);
end

hold on
% grid on
title(figure_title)

% plot
if show_hist
    scatterhist(x,y,'Color','k','NBins',20)
    hold on
else
    scatter(x,y,30,'k') 
end

Pearson_cor = corr(x,y,'Type','Pearson','rows','complete');
Kendall_cor = corr(x,y,'Type','Kendall','rows','complete');
Spearman_cor = corr(x,y,'Type','Spearman','rows','complete');

% do linear regression if requested
if show_fit
    P = polyfit(x,y,1);
    yfit = polyval(P,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    
    set(0,'CurrentFigure',fig(1));
    plot(x,yfit,'-','color',colour_mat(end/2,:),'linewidth',3);
    equation = sprintf(...
        'y = %.2f x + %.2f \nR^2 = %.2f \nPearson = %.2f \nSpearman = %.2f',...
        P(1),P(2),rsq,Pearson_cor,Spearman_cor); %,Kendall_cor \nKendall = %.2f 
    if show_corr
        if show_hist
            annotation('textbox',[0.025, 0.035, 0.25, 0.2],...
                'string',equation,'EdgeColor','none')
        else
            TextLocation(equation,'Location',[0.52 0.31 0.46 0.17]);
        end
    end
elseif show_corr
    set(0,'CurrentFigure',fig(1));
    equation = sprintf(...
        'Pearson = %.2f \nSpearman = %.2f',...
        Pearson_cor,Spearman_cor); %,Kendall_cor \nKendall = %.2f 
    TextLocation(equation,'Location',[0.52 0.18 0.46 0.17]);
end

xlim(x_limits)
ylim(y_limits)
% axis equal
xlabel(x_name)
ylabel(y_name)

% update cursor
dcm_obj = datacursormode(figure(fig(1)));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,index})

%% save fig
if save_figure
    fig_name = strcat('bivariate_',figure_name);  
    saveFig(fig,fig_name,figure_path,figure_type)
end

end

function hOut = TextLocation(textString,varargin)

l = legend(textString,varargin{:});
t = annotation('textbox');
t.String = textString;
t.Position = l.Position;
delete(l);
t.LineStyle = 'None';
t.FontSize = 8;

if nargout
    hOut = t;
end
end