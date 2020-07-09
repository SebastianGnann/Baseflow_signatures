function [] = plotTrivariate(x,y,z,varargin)
%plotTrivariate Plots scatter plot of two variables with the dots coloured
%   according to a third variable.
%
%   INPUT
%   x: variable 1
%   y: variable 2
%   z: variable 3
%   OPTIONAL
%   x_name: string with name of x
%   y_name: string with name of y
%   z_name: string with name of z
%   ID: if data cursor shall show ID input the respective array
%   x_limit: axis limits for x-axis
%   y_limit: axis limits for y-axis
%   z_limits: limits of colour axis, e.g. [0 1]
%   z_lower_limit_open: is the lower limit open?
%   z_upper_limit_open: is the upper limit open?
% 	colour_scheme: name of colour scheme
%   flip_colour_scheme: flip colour scheme?
%   nr_colours: nr of colours used for colourscale
%   figure_title: title of plot, e.g. '(a)'
%   figure_name: name for saving, e.g. UK_P_Q
%   save_figure: save plot true/false
%   figure_path: path to folder where figure should be saved
%   figure_type: figure type, e.g. -dpdf or -dmeta
%
%   OUTPUT
%   plot and saved figure
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

if nargin < 3
    error('Not enough input arguments.')
end

ip = inputParser;

addRequired(ip, 'x', ...
    @(x) isnumeric(x) && (size(x,1)==1 || size(x,2)==1))
addRequired(ip, 'y', ...
    @(y) isnumeric(y) && (size(y,1)==1 || size(y,2)==1))
addRequired(ip, 'z', ...
    @(z) isnumeric(z) && (size(z,1)==1 || size(z,2)==1))

addParameter(ip, 'x_name', [], @ischar)
addParameter(ip, 'y_name', [], @ischar)
addParameter(ip, 'z_name', [], @ischar)
addParameter(ip, 'ID', NaN(size(x)), @isnumeric)
addParameter(ip, 'x_limits', [min(x) max(x)], @(x_limits) isnumeric(x_limits) && length(x_limits)==2)
addParameter(ip, 'y_limits', [min(y) max(y)], @(y_limits) isnumeric(y_limits) && length(y_limits)==2)
addParameter(ip, 'z_limits', [min(z) max(z)], @(z_limits) isnumeric(z_limits) && length(z_limits)==2)
addParameter(ip, 'z_lower_limit_open', false, @islogical)
addParameter(ip, 'z_upper_limit_open', false, @islogical)
addParameter(ip, 'nr_colours', 10, @isnumeric)
addParameter(ip, 'colour_scheme', 'Spectral', @ischar)
addParameter(ip, 'flip_colour_scheme', false, @islogical)
addParameter(ip, 'figure_title', '', @ischar)
addParameter(ip, 'figure_name', 'no_name', @ischar)
addParameter(ip, 'save_figure', false, @islogical)
addParameter(ip, 'figure_path', '', @ischar)
addParameter(ip, 'figure_type', '-dpdf', @ischar)

parse(ip, x, y, z, varargin{:})

x_name = ip.Results.x_name;
y_name = ip.Results.y_name;
z_name = ip.Results.z_name;
ID = ip.Results.ID;
x_limits = ip.Results.x_limits;
y_limits = ip.Results.y_limits;
z_limits = ip.Results.z_limits;
z_lower_limit_open = ip.Results.z_lower_limit_open;
z_upper_limit_open = ip.Results.z_upper_limit_open;
nr_colours = ip.Results.nr_colours;
colour_scheme = ip.Results.colour_scheme;
flip_colour_scheme = ip.Results.flip_colour_scheme;
figure_title = ip.Results.figure_title;
figure_name = ip.Results.figure_name;
save_figure = ip.Results.save_figure;
figure_path = ip.Results.figure_path;
figure_type = ip.Results.figure_type;

%% plotting
% create colormap
if flip_colour_scheme
    colour_mat = flip(brewermap(nr_colours,colour_scheme));
else
    colour_mat = brewermap(nr_colours,colour_scheme);
end

% get rid of NaN points
xNaN = isnan(x); x(xNaN) = []; y(xNaN) = []; z(xNaN) = [];
yNaN = isnan(y); x(yNaN) = []; y(yNaN) = []; z(yNaN) = [];
zNaN = isnan(z); x(zNaN) = []; y(zNaN) = []; z(zNaN) = [];
index = [1:length(ID)]';
ID(xNaN) = []; ID(yNaN) = []; ID(zNaN) = [];
index(xNaN) = []; index(yNaN) = []; index(zNaN) = [];

fig = figure('Name',figure_name,'NumberTitle','off','pos',[10 10 300 230]);
hold on
% grid on

% plot
colormap(colour_mat)
scatter(x,y,25,z,'filled')
dcm_obj = datacursormode(figure(fig));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,index})
xlabel(x_name)
ylabel(y_name)
% axis equal
xlim(x_limits)
ylim(y_limits)
title(figure_title)
c = colorbar;
title(c,z_name)
x1=get(gca,'position');
x=get(c,'Position');
x(3)=10/400;
set(c,'Position',x)
set(gca,'position',x1)
caxis(z_limits)
if z_lower_limit_open
    c.TickLabels{1} = ['<' c.TickLabels{1}];
end
if z_upper_limit_open
    c.TickLabels{end} = ['>' c.TickLabels{end}];
end

%% save fig
if save_figure
    fig_name = strcat('trivariate','_',figure_name);
    saveFig(fig,fig_name,figure_path,figure_type)
end

end

