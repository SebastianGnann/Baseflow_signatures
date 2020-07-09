function [data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata]  = read_ascii_grid(filename)
%read_ascii_grid - from Gemma Coxon

fid = fopen(filename, 'r');

tmp = textscan(fid, '%s%f', 1);
ncols = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
nrows = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
xllcorner = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
yllcorner = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
cellsize = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
nodata = tmp{1,2};

formatString = repmat('%f', 1, ncols);
data = textscan(fid, formatString, nrows);
data = cell2mat(data);

fclose(fid);

end