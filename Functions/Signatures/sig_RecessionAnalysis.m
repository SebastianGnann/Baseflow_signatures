function [para_mat, recession_month] = sig_RecessionAnalysis(Q, t, varargin)
%sig_RecessionAnalysis Fits power law function to individual recession 
%   segments and returns recession parameters (see Brutsaert and Nieber,
%   1977; Roques et al., 2017; and Jachens et al., 2020)
%   dQ/dt = - a Q ^ b
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   recession_length: min. length of recession segments (days), default = 5
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default=0
%   start_of_recession: define start of recession when baseflow filter
%   rejoins the curve (reference) "baseflow" or after peak "peak"
%   filter_par: smoothing parameter of Lyne Hollick Filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   plot_results: whether to plot results, default = 0
%   fitIndividual: fit each individual recession segment
%   fitting_type: fit non-linear or linear curve ('nonlinear','linear')
%   reservoir), etc.
%   method: method for dQ/dt calculation
%
%   OUTPUT
%   para_mat: matrix with parameters alpha, beta (=1 for exponential fit)
%   for each recession segment
%   recession_month: approx. month of recession
%
%   References
%	Brutsaert, W. and Nieber, J.L., 1977. Regionalized drought flow 
%   hydrographs from a mature glaciated plateau. Water Resources Research, 
%   13(3), pp.637-643.
%   Roques, C., Rupp, D.E. and Selker, J.S., 2017. Improved streamflow 
%   recession parameter estimation with attention to calculation of? dQ/dt.
%   Advances in water resources, 108, pp.29-43.
%   Jachens, E.R., Rupp, D.E., Roques, C. and Selker, J.S., 2020. Recession 
%   analysis revisited: Impacts of climate on parameter estimation. 
%   Hydrology and Earth System Sciences, 24(3), pp.1159-1170.
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
addParameter(ip, 'ignoreNaN', 'y', @ischar) % ignore NaN values y/n?
addParameter(ip, 'recession_length', 5, @isnumeric) % length of decreasing flow section (amount of timesteps) to be declared a recession
addParameter(ip, 'n_start', 1, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'peak', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of Lyne Hollick Filter to determine start of recession (higher = later recession start)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)
addParameter(ip, 'fitIndividual', true, @islogical) % fit individual recessions or point cloud
addParameter(ip, 'fitting_type', 'linear', @ischar) % nonlinear or linear fit
addParameter(ip, 'dQdt_method', 'ETS', @ischar) % how to calculate dQ/dt

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
filter_par = ip.Results.filter_par;
plot_results = ip.Results.plot_results;
fitIndividual = ip.Results.fitIndividual;
fitting_type = ip.Results.fitting_type;
dQdt_method = ip.Results.dQdt_method;

% data checks
[~,~,t] = util_DataCheck(Q, t);

% get recession segments
[flow_section] = util_RecessionSegments(Q, t, ...
    'recession_length',recession_length,'eps',eps,...
    'filter_par',filter_par,'plot_results',false,...
    'start_of_recession',start_of_recession,'n_start',n_start);

% get flow rate gradient and corresponding flows
[dQdt, Qm, flow_section, R2] = ...
    util_dQdt(Q, t, flow_section, 'method', dQdt_method);

% get recession month
date_tmp = datevec(floor(mean(flow_section,2)));
recession_month = date_tmp(:,2);

if ~fitIndividual
    rec = ~isnan(Qm);
    [para_mat(1), para_mat(2)] = ...
        util_fitPowerLaw(Qm(rec), dQdt(rec),...
        'fitting_type', fitting_type, 'weights', R2(rec));
    
else
    para_mat = NaN(length(flow_section),2);
    for i = 1:size(flow_section,1)
        rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
        [para_mat(i,1), para_mat(i,2)] = ...
            util_fitPowerLaw(Qm(rec), dQdt(rec), ...
            'fitting_type', fitting_type, 'weights', R2(rec));
    end
end

if plot_results
    %     figure('pos',[10 10 300 250])
    hold on
    colour_mat_seasons = [...
        0 0 1;  0 0 1;...
        0 1 0; 0 1 0; 0 1 0;...
        1 0 0; 1 0 0; 1 0 0;...
        1 1 0; 1 1 0; 1 1 0; ...
        0 0 1];
    p1=plot(0,0,'.','Color',[0 1 0]); 
    p2=plot(0,0,'.','Color',[1 0 0]);  
    p3=plot(0,0,'.','Color',[1 1 0]); 
    p4=plot(0,0,'.','Color',[0 0 1]); 
    for i = 1:size(flow_section,1)
        rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
        Q_tmp = Qm(rec);
        dQdt_tmp = dQdt(rec);
        date_vec = datevec(t(rec));
        ind = floor(median(date_vec(:,2))); % get approx. month
        plot(Q_tmp,-dQdt_tmp,'.','color',colour_mat_seasons(ind,:),'linewidth',2)
        if fitIndividual
            plot(Q_tmp,para_mat(i,1).*Q_tmp.^para_mat(i,2),'color',colour_mat_seasons(ind,:))
        end
    end
    legend([p1 p2 p3 p4],{'MAM','JJA','SON','DJF'},'box','off','Location','best');
    
    if ~fitIndividual
        rec = ~isnan(Qm);
        plot(sort(Qm(rec)),para_mat(1).*sort(Qm(rec)).^para_mat(2),...
            '-','color',[0.3 0.3 0.3],'linewidth',2);
        str = (sprintf('-dQ/dt = %.2f Q^{%.1f} \n',para_mat(1),para_mat(2)));
        title(str);
    else
        [~, ind] = min(abs(para_mat(:,2) - median(para_mat(:,2)))); % find recession according to median exponent
        str = (sprintf('-dQ/dt = %.3f Q^{%.1f}',para_mat(ind,1),para_mat(ind,2)));
        title(str);
    end
    xlabel('Q') % [mm/d]
    ylabel('-dQ/dt') % [mm/d^2]
    set(gca,'XScale','log')
    set(gca,'YScale','log')
end

end

