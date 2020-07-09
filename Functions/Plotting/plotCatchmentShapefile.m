function [] = plotCatchmentShapefile(...
    catchment_tmp,GloRiC,x_flat_tmp,y_flat_tmp,data_flat_tmp,...
    attribute_name,i,save_figure)
%plotCatchmentShapefile Plots shapefile of catchment with streams and
%   attributes.
%
%   INPUT
%   catchment_tmp: catchment shapefile
%   GloRiC: river shapefile
%   x_flat_tmp: x data
%   y_flat_tmp: y data
%   data_flat_tmp: data
%   attribute_median: name of attribute
%   i: iteration
%   save_figure: save plot true/false
%
%   OUTPUT
%   plot
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

% check input parameters
if nargin < 7
    error('Not enough input arguments.')
end
if nargin < 8
    save_figure = false;
end

ID_tmp = catchment_tmp.hru_id;
ID_GloRiC = [GloRiC.hru_id]';
try
    [ism, index] = ismember(ID_tmp,ID_GloRiC);
    river_tmp = GloRiC(index);
catch
    disp('No river data.')
end

colour_mat = (brewermap(10,'Spectral'));

fig = figure();
hold on
xlabel('Easting [km]')
ylabel('Northing [km]')
title(sprintf('hru id: %.0f / i: %.0f', catchment_tmp.hru_id, i))
%     h = pcolor(X_tmp,Y_tmp+cellsize,data_tmp);
%     set(h, 'EdgeColor', 'none')
plot_cellsize = 1e5/(length(x_flat_tmp)+length(y_flat_tmp));
h = scatter(x_flat_tmp./1000,y_flat_tmp./1000,plot_cellsize,(data_flat_tmp),'s','filled');%,'markeredgecolor','k');
axis equal
colormap(colour_mat)
c = colorbar;
%     caxis([0 1000])
title(c,attribute_name, 'Interpreter', 'none')
plot(catchment_tmp.X./1000,catchment_tmp.Y./1000,'k','linewidth',1.5)
try
    plot(river_tmp.X./1000,river_tmp.Y./1000,'b','linewidth',1.5)
catch
end

%% save fig
if save_figure
    set(fig,'Units','Inches');
    position = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[position(3),position(4)]);
    fig_name = strcat('shapefile_',attribute_name);
    path_name = './Baseflow_and_geology/Images';
    path = strcat(path_name,'\',fig_name);
    print(fig,path,'-dpdf','-r500');
end

end