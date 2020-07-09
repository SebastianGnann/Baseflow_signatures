function [Q_mean] = sig_Q_mean(Q, t)
%sig_Q_mean Mean of flow time series with variable timesteps.
%   Calculates mean flow with the following options:
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%	OPTIONAL
%
%   OUTPUT
%   Q_mean: mean flow [mm/timestep]
%
%   References
%   https://en.wikipedia.org/wiki/Arithmetic_mean
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

parse(ip, Q, t)

% data checks
[~,~,t] = util_DataCheck(Q, t);

% calculate signature
Q_mean = mean(Q,'omitnan'); % always ignoring NaNs for now
    
end
