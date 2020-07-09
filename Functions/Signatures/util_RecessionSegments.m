function [flow_section] = util_RecessionSegments(Q, t, varargin)
%util_RecessionSegments Identify all individual baseflow recession segments
% (see e.g. Safeeq et al., 2013, Jachens et al., 2020).
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%
%	OPTIONAL
%   recession_length: min. length of recessions (days),default = 15
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default=0
%   start_of_recession: define start of recession when baseflow filter
%   rejoins the curve (reference) "baseflow" or after peak "peak"
%   filter_par: smoothing parameter of Lyne Hollick Filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   plot_results: whether to plot results, default = 0
%
%   OUTPUT
%   flow_section: n-by-2 array where n is the number of recession segments;
%   columns are the indices into the flow array of the start and end of 
%   the recession segments
%
%   References 
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
addParameter(ip, 'recession_length', 15, @isnumeric) % length of decreasing flow in days to be declared a recession
addParameter(ip, 'n_start', 0, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'baseflow', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of Lyne Hollick Filter to determine start of recession (higher = later recession start)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
plot_results = ip.Results.plot_results;
filter_par = ip.Results.filter_par;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
n_start = ip.Results.n_start;

% Identify all individual recession segments with length > recession_length days.
% how many decreasing timesteps depends on length of timestep
len_decrease = recession_length/days(t(2)-t(1));
% find timesteps with decreasing flow
decreasing_flow = Q(2:end)<(Q(1:end-1)+eps); %SG: changed from <= to <
% start on a non-decreasing point
start_point = find(decreasing_flow==0,1);
decreasing_flow = decreasing_flow(start_point:end);
% find start and end of decreasing sections
flow_change = find(decreasing_flow(1:end-1) ~= decreasing_flow(2:end));
% reshape into x by 2 array (columns = start, end of decrease)
flow_change = flow_change(1:(2*floor(size(flow_change,1)./2)));
flow_change = reshape(flow_change,2,[]).';
% find sections
flow_section = flow_change((flow_change(:,2)-flow_change(:,1))>=len_decrease+n_start,:);
flow_section = flow_section+start_point;
flow_section(:,1) = flow_section(:,1)+n_start; % move start point n days
if numel(flow_section)==0
    error('No long enough recession periods, consider setting eps parameter > 0')
end

if plot_results
    [dQdt,Qm,flow_section] = util_dQdt(Q, t, flow_section);
    figure('Position',[300 300 1000 400]); hold on;
    plot(t,Q);
    for i = 1:size(flow_section,1)
        h1=plot(t(flow_section(i,1):flow_section(i,2)),Qm(flow_section(i,1):flow_section(i,2)),'r-');
        plot(t(flow_section(i,1):flow_section(i,2)),dQdt(flow_section(i,1):flow_section(i,2)),'g');
    end
    title('Selected recession segments')
%     datetick('x')
    ylabel('Flow')
end

switch start_of_recession
    case 'peak'
    case 'baseflow'
        % The beginning of recession = point where baseflow filter rejoins the curve
        % (alternative needs basin area)
        % Use Lyne Hollick Filter with default a, #passes = 1 (doesn't work with > 1
        % as then baseflow is never equal to flow)
        [B] = util_LyneHollickFilter(Q,'filter_parameter',filter_par,'nr_passes',1);
        
        % Find sections identified as baseflow in complete series
        isbaseflow = B == Q;
        for i = 1:size(flow_section,1)
            % Find which points in the decreasing section count as baseflow only
            isb_section = isbaseflow(flow_section(i,1):flow_section(i,2));
            if numel(isb_section)==0
                flow_section(i,:) = NaN;
            else
                % Find the first such point
                isb_start = find(isb_section==1,1,'first');
                % Limit section to only those points
                flow_section(i,1)=flow_section(i,1)+isb_start-1;
            end
        end
        % if flow_section less than 4 points, remove
        flow_section = flow_section(~isnan(flow_section(:,1)),:);
        if numel(flow_section)==0
            error('No long enough baseflow recession periods, consider increasing filter_par parameter')
        end
        
        %Add the baseflow recession sections to the plot
        if plot_results
            %Add baseflow to plot
            plot(B, 'k--')
            for i = 1:size(flow_section,1)
                h2=plot(t(flow_section(i,1):flow_section(i,2)),Q(flow_section(i,1):flow_section(i,2)),'g-');
            end
            legend([h1 h2],{'Complete recessions', 'Baseflow recessions'})
        end
    otherwise
        error('Incorrect option for start of recession.')
end

end