function [dQdt, Qm, flow_section, R2] = util_dQdt(Q, t, flow_section, varargin)
%util_dQdt Calculates flow rate gradient with the following options:
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datenum]
%   flow_section: n-by-2 array where n is the number of recession segments
%
%	OPTIONAL
%   method: method for dQdt calculation
%
%   OUTPUT
%   dQdt: flow rate gradient
%   Qm: corresponding flow
%   flow_section: updated flow_section array (some recession points have to
%   be removed due to approx. of derivative)
%   R2: R^2 from exponential time stepping method
%
%   References
%	Brutsaert, W. and Nieber, J.L., 1977. Regionalized drought flow 
%   hydrographs from a mature glaciated plateau. Water Resources Research, 
%   13(3), pp.637-643.
%   Roques, C., Rupp, D.E. and Selker, J.S., 2017. Improved streamflow 
%   recession parameter estimation with attention to calculation of? dQ/dt.
%   Advances in water resources, 108, pp.29-43.
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
addRequired(ip, 'flow_section', @(flow_section) isnumeric(flow_section) && size(flow_section,2)==2) 

% optional input arguments
addParameter(ip, 'method', 'ETS', @ischar) 

parse(ip, Q, t, flow_section, varargin{:})
method = ip.Results.method;

dQdt = NaN(size(Q));
Qm = NaN(size(Q));
m = NaN(size(Q));
R2 = ones(size(Q)); % weights

switch method    
        
    case 'BN' % Brutsaert and Nieber (1979)
        for j = 1:size(flow_section,1)
            rec = [flow_section(j,1):flow_section(j,2)]'; % get recession
            dQdt(rec(2:end)) = (Q(rec(2:end))-Q(rec(1:end-1)))./days(t(2)-t(1));
            Qm(rec(2:end)) = (Q(rec(2:end)) + Q(rec(1:end-1)))./2;            
            flow_section(j,1) = flow_section(j,1)+1; % shorten recession
        end
       
    case 'backwards' % similar to Brutsaert and Nieber (1979), but we keep measured Q, see also Thomas et al. (2015)
        for j = 1:size(flow_section,1)
            rec = [flow_section(j,1):flow_section(j,2)]'; % get recession
            dQdt(rec(2:end)) = (Q(rec(2:end))-Q(rec(1:end-1)))./days(t(2)-t(1));
            Qm(rec(2:end)) = Q(rec(2:end));
            flow_section(j,1) = flow_section(j,1)+1; % shorten recession
        end
        
    case 'ETS' % exponential time stepping following Roques et al. (2017)        
        for j = 1:size(flow_section,1)
            rec = [flow_section(j,1):flow_section(j,2)]'; % get recession
            n = 0.2*(length(rec)); % n = 10% of recession led to good results according to Roques et al. (2017) % 20% as in Roques examples
            gamma = fitExponential(Q(rec), t(rec)); % get gamma
            m(rec) = 1 + ceil(n.*exp(-1./(gamma.*[1:length(rec)]))); % multiply by time step?
            
            i = rec(1);
            while i+m(i) <= rec(end)
                Qm(i) = mean(Q(i:i+m(i)));
                [~, dQdt(i), R2(i)] = fitLinear(...
                    datenum(t(i:i+m(i))),Q(i:i+m(i)));
                i = i+1;
            end 
            
            flow_section(j,2) = flow_section(j,2) - m(i);%(m(i)+1); % shorten recession
            
        end
        
    otherwise
        error('Differentiation method not available.')
end

% dQdt has to be negative and weights have to be non-negative
dQdt(dQdt>=0) = NaN;
Qm(isnan(dQdt)) = NaN;
R2(R2<=0) = 0; 
R2(isnan(R2)) = 0;

end

function [gamma] = fitExponential(Q, t)
%fitExponential Fits an exponential function to recession segments.
%   Q = Q0*exp(-gamma*t)
t = [1:length(t)]';
ExponentialObjective = @(gamma) Q(1).*exp(-gamma.*t) - Q;
gamma0 = [0.1];
options = optimoptions(@lsqnonlin,'Display','off');
gamma = lsqnonlin(ExponentialObjective, gamma0, [1e-9], [100], options);

end

function [a, b, R2] = fitLinear(x,y)
%fitLinear fits linear function and returns parameters and residuals.
%   y = a + b*x;

n = length(x);
SSxy = sum(x.*y) - sum(x)*sum(y)/n;
SSxx = sum(x.^2) - sum(x)^2/n;
b = SSxy/SSxx;
a = mean(y) - b*mean(x);
y_hat = a + b*x;
R2 = 1 - sum((y - y_hat).^2)./sum((y - mean(y)).^2);

end


