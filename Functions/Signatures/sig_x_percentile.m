function [Q_x] = sig_x_percentile(Q, t, x, varargin)
%sig_x_percentile Calculate x-th flow percentile of streamflow.
%   Following Addor et al. (2018) Q95 is a high flow measure, i.e. the 95%
%   NON-exceedance probability. 
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   x: x-th percentile(s) (e.g. 95 for Q95), can also be a vector
%	OPTIONAL
%
%   OUTPUT
%   Q_x: x-th flow percentile (flow that is not reached for x of the time)
%
%   References
%   Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and
%   Clark, M.P., 2018. A ranking of hydrological signatures based on their 
%   predictability in space. Water Resources Research, 54(11), pp.8792-8812.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 3
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1)) 
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1)) 
addRequired(ip, 'x', @(x) isnumeric(x) && (size(x,1)==1 || size(x,2)==1)) 

% optional input arguments

parse(ip, Q, t, x, varargin{:})

% data checks
[~,~,t] = util_DataCheck(Q, t);

if any(x>100) || any(x<0)
    error('x must be between 0 and 100.')
end

% calculate signature
p = 1 - x./100; % transform to get exceedance probability

% get ranks as a proxy for exceedance probabilities   
Q_sorted = sort(Q);
% Q_ranked = tiedrank(Q_sorted);
Q_ranked = [1:length(Q)]'; % give unique (random) rank to every measurement
FDC = 1 - Q_ranked./length(Q_ranked); % flow duration curve

% find x-th flow percentile
indices = 1:length(FDC);
bound_x = NaN(size(p));
for i = 1:length(p)
    % if flow is highly ephemeral, FDC might not be well defined
    if isempty(max(indices(FDC >= p(i))))        
    else
        bound_x(i) = max(indices(FDC >= p(i)));
    end
end

Q_x = NaN(size(p)); 
Q_x(~isnan(bound_x)) = Q_sorted(bound_x(~isnan(bound_x)));
    
end
