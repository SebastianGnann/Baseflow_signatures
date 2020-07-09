function [a, b] = util_fitPowerLaw(Q_rec, dQdt_rec, varargin)
%util_fitPowerLaw Fits a power law (to recession segments).
%   Can either be done by linear fitting in loglog-space or non-linear
%   fitting.
%   dQ/dt = - a Q ^ b
%
%   INPUT
%   Q_rec: streamflow (recession periods only)
%   dQdt_rec: corresponding flow rate gradients
%
%   OPTIONAL
%   fitting_type: specifies if fitting procedure is linear or non-linear
%   weights: weights to fit curve
%
%   OUTPUT
%   a: scaling parameter
%   b: parameter of non-linearity
%
%   References
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 2    
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q_rec', @(Q_rec) isnumeric(Q_rec) && (size(Q_rec,1)==1 || size(Q_rec,2)==1))
addRequired(ip, 'dQdt_rec', @(dQdt_rec) isnumeric(dQdt_rec) && (size(dQdt_rec,1)==1 || size(dQdt_rec,2)==1)) 

% optional input arguments
addParameter(ip, 'fitting_type', 'linear', @ischar) %
addParameter(ip, 'weights', NaN(size(Q_rec)), @isnumeric) % 

parse(ip, Q_rec, dQdt_rec, varargin{:})
fitting_type = ip.Results.fitting_type;
weights = ip.Results.weights;

if strcmp(fitting_type, 'linear') % linear regression in log log space
        
    p0 = [0.1, 1.0];
    fit_lin = fit(log(Q_rec),log(-dQdt_rec),'a + b.*x','Start',p0,'Weight',weights);
    p = coeffvalues(fit_lin);
    a = exp(p(1));
    b = p(2);       
    
elseif strcmp(fitting_type, 'nonlinear') % nonlinear regression in lin space
    
    powerFcn = @(p,x) p(1).*x.^p(2);
    p0 = [0.1, 1.0];
    fit_nonlin = fitnlm(Q_rec,-dQdt_rec,powerFcn,p0,'Weight',weights);
    a = fit_nonlin.Coefficients.Estimate(1);
    b = fit_nonlin.Coefficients.Estimate(2);   
        
else
    error('Invalid fitting type.')
end

end


