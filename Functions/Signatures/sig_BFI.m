function [BFI] = sig_BFI(Q, t, varargin)
%sig_BFI Calculates baseflow index (BFI).
%   Calculates BFI, that is the ratio between baseflow (volume) and
%   total streamflow (volume), with the following options:
%       - choose different baseflow separation methods (Lyne and Hollick, 
%       1979; UK Institute of Hydrology, 1980)
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%	OPTIONAL
%   method: which baseflow separation method should be employed
%   ('Lyne_Hollick','UKIH')
%   parameters: specify filter parameters ([filter_parameter nr_passes] for
%   Lyne Hollick and [n_days] for UKIH)
%
%   OUTPUT
%   BFI: baseflow index [-]
%
%   References
%   Lyne, V. and Hollick, M., 1979. Stochastic time-variable
%   rainfall-runoff modelling. In Institute of Engineers Australia National
%   Conference (Vol. 1979, pp. 89-93). Barton, Australia: Institute of
%   Engineers Australia.
%   Institute of Hydrology (Great Britain), 1980. Low Flow Studies Reports.
%   Institute of Hydrology.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1)) 
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1)) 

% optional input arguments
addParameter(ip, 'method', 'Lyne_Hollick', @ischar) % Which method? Default: Lyne_Hollick
addParameter(ip, 'parameters', [], @isnumeric) % Which parameter values?
addParameter(ip, 'threshold_type', [], @ischar) % How to threshold Lyne-Hollick filter?

parse(ip, Q, t, varargin{:})
method = ip.Results.method;
parameters = ip.Results.parameters;
threshold_type = ip.Results.threshold_type;

% data checks
[~,~,t] = util_DataCheck(Q, t);

% calculate signature

% pad time series to compensate for warm up effect (Ladson et al., 2013)
if length(Q)>60
    Q_padded = [Q(30:-1:1); Q; Q(end-29:end)];
else
    Q_padded = Q;
    warning('Very short time series. Baseflow separation might be unreliable.')
end
% obtain baseflow
switch method
    
    case 'Lyne_Hollick'
        if isempty(parameters)
            parameters = [0.925, 3];
        elseif length(parameters) == 1
            parameters(2) = 3;
        elseif length(parameters) > 2
            error('Too many filter parameters.')
        end
        
        if isempty(threshold_type)
            B = util_LyneHollickFilter(Q_padded, ...
                'filter_parameter', parameters(1), 'nr_passes', parameters(2));
        else
            B = util_LyneHollickFilter(Q_padded, ...
                'filter_parameter', parameters(1), 'nr_passes', parameters(2),...
                'threshold_type',threshold_type);
        end
        
    case 'UKIH'
        if isempty(parameters)
            parameters = 5;
        elseif length(parameters) > 1
            error('Too many filter parameters.')
        end
        B = util_UKIH_Method(Q_padded, 'n_days', parameters(1));
        
    otherwise
        error('Please choose one of the available baseflow separation methods.')
end

% remove padding
if length(Q)>60
    B(1:30) = [];
    B(end-29:end) = [];
else
end

% calculate BFI
BFI = sum(B,'omitnan')./sum(Q,'omitnan');

% check if 0<BFI<=1
if BFI<0 || BFI>1
    BFI = NaN;
end

end
