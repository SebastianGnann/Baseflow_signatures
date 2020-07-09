function [B] = util_LyneHollickFilter(Q, varargin)
%util_LyneHollickFilter Calculates baseflow using the Lyne-Hollick filter.
%   Calculates baseflow using the Lyne and Hollick recursive digital filter
%   (Lyne and Hollick, 1979).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   filter_parameter: filter parameter (default = 0.925)
%   nr_passes: number of passes (default = 3; forwards, backwards, forwards)
%   threshold_type: how to threshold resulting time series (default = at
%   the end of all passes (end); other options are after each pass (pass), 
%   or after each timestep (timestep), or no thresholding (none))
%
%   OUTPUT
%   B: baseflow [mm/timestep]
%
%   References
%   Lyne, V. and Hollick, M., 1979. Stochastic time-variable
%   rainfall-runoff modelling. In Institute of Engineers Australia National
%   Conference (Vol. 1979, pp. 89-93). Barton, Australia: Institute of
%   Engineers Australia.
%   Su, C.H., Costelloe, J.F., Peterson, T.J. and Western, A.W., 2016. On
%   the structural limitations of recursive digital filters for base flow
%   estimation. Water Resources Research, 52(6), pp.4745-4764.
%   Ladson, A.R., Brown, R., Neal, B. and Nathan, R., 2013. A standard
%   approach to baseflow separation using the Lyne and Hollick filter.
%   Australasian Journal of Water Resources, 17(1), pp.25-34.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.


% check input parameters
if nargin < 1
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1)) 

% optional input arguments
addParameter(ip, 'filter_parameter', 0.925, @isnumeric)
addParameter(ip, 'nr_passes', 3, @isnumeric)
addParameter(ip, 'threshold_type', 'end', @ischar)

parse(ip, Q, varargin{:})
filter_parameter = ip.Results.filter_parameter;
nr_passes = ip.Results.nr_passes;
threshold_type = ip.Results.threshold_type;

if filter_parameter>1 && filter_parameter<=0
    error('Filter parameter must be between 0 and 1.')
end

if floor(nr_passes)~=nr_passes && nr_passes<1
    error('Number of filter passes must be an integer larger than zero.')
end

% Baseflow separation is problematic with NaN values. Therefore, we set NaN
% values to median, apply the filter, and then set baseflow to NaN where
% streamflow is NaN. If there are a lot of NaN values, we encourage the
% user to either interpolate these values or to calculate the signature for
% each block individually and then calculate a weighted average.
Q_tmp = Q;
Q_tmp(isnan(Q)) = median(Q,'omitnan');

% calculate baseflow by applying RDF several times
B = LyneHollickFilter(Q_tmp, filter_parameter, threshold_type);
for nr = 2:nr_passes
    B = LyneHollickFilter(flip(B), filter_parameter, threshold_type);
end

% set baseflow to NaN where streamflow is NaN
B(isnan(Q)) = NaN;

% constrain baseflow not to be higher than streamflow (see also Ladson et
% al. (2013)) at the end
if strcmp(threshold_type,'none')
else
    B(B>Q) = Q(B>Q);
end

end

function B = LyneHollickFilter(Q, filter_parameter, threshold_type)
%LyneHollickFilter Helper function that runs the Lyne-Hollick filter.

% define thresholding method
threshold_timestep = false;
threshold_pass = false;
switch threshold_type
    case 'end'
    case 'timestep'
        threshold_timestep = true;
    case 'pass'
        threshold_pass = true;
    case 'none'
    otherwise
        error('Not a valid thresholding method. Choose either end, timestep, pass, or none.')
end

n = length(Q);
F = NaN(n,1);
F(1) = Q(1) - min(Q); % initial condition, see Su et al. (2016)

if threshold_timestep
    for i=2:1:n
        F(i) = filter_parameter*F(i-1) + ((1+filter_parameter)/2)*(Q(i) - Q(i-1));
        if F(i)<0 % constrain after each timestep
            F(i) = 0;
        end
    end
else
    for i=2:1:n
        F(i) = filter_parameter*F(i-1) + ((1+filter_parameter)/2)*(Q(i) - Q(i-1));
    end
end

if threshold_pass
    F(F<0) = 0; % constrain after each filter pass
end

% calculate baseflow
B = Q - F;

end
