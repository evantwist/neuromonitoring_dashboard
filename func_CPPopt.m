%% Function for CPPopt analysis
function [CPPopt] = func_CPPopt(CPP, PRx)
% Input:
%       CPP: CPP time series (optionally filtered to remove noise),
%       can also be calculated as the difference between synchronized MAP
%       and ICP time series. Comes with a time axis t.
%       PRx: PRx as calculated by func_PRx (PRx.data)
%
% Output: 
%       CPPopt_trend: the weighted average CPP over all time windows,
%       calculated as the mean of the mean CPPopt and the CPPopt with the
%       lowest PRx (i.e. intact autoregulation), updated per minute.
%       CPPopt_trend_display: contains the CPPopt trend, but in case of
%       impaired autoregulation (PRx > 0.2) the CPPopt value is rejected
%       and replaced by a missing value.
%       range_cpp: optimal range of CPP determined as minimum PRx (ll_prx) 
%       to minimum PRx + 0.2 (ul_prx).
%
% Inspired by Donnely et al. and Depreitere et al.
% Adjusted from Bart Formsma and Tahisa Robles
% Erasmus MC Sophia Childrens Hospital, Pediatric Intensive Care Unit in
% collaboration with Delft University of Technology, 3mE 

% Disclaimer: not validated for real time clinical use
% Eris van Twist, May 2024.

close all
clear all
clc

CPPopt = struct();
%% Preprocessing
% Fixed parameters
Ts = 10;                                                                    % Seconds per sample
window_ave_CPP = 5*60/(2*Ts);                                               % Averaging window for CPP (5 min)
window_ave_PRx = 1*60/(2*Ts);                                               % Averaging window for PRx (1 min)
bin_width = 5;                                                              % CPP bins in mmHg
window = [1*60 2*60 4*60 6*60 8*60];                                        % CPPopt windows of 1, 2, 4, 6 and 8 hours in minutes

% Calculate mean CPP per minute
for bb = window_ave_CPP+1:(60/Ts):(length(CPP)-window_ave_CPP)            
    CPP_av1(bb) = mean(CPP(bb-window_ave_CPP:bb+window_ave_CPP),'omitnan');
end

% Samples per minute
CPP_av = CPP_av1(window_ave_CPP+1:(60/Ts):(length(CPP)-window_ave_CPP)); 

% Calculate mean PRx per minute
for cc = window_ave_CPP+1:(60/Ts):(length(CPP)-window_ave_CPP)
    PRx_av1(cc) = mean(PRx(cc-window_ave_PRx:cc+window_ave_PRx),'omitnan');
end

% Samples per minute
PRx_av = PRx_av1(window_ave_CPP+1:(60/Ts):(length(CPP)-window_ave_CPP));

%% Calculations of CPPopt per minute
% Initialize vectors 
CPPopt_trend = [];
CPPopt_loop = [];
CPPopt_trend_display = [];

% Add new cell for every minute containing the bin edges where CA is intact
range_cpp = cell(length((8*60+1):1:length(CPP_av)), 1); 

%% Loop over windows to calculate CPPopt and optimal CPP range
% Calculations are done per minute, start with 8 hour window
for i = (8*60+1):1:length(CPP_av)
    CPPopt_loop = [];   
    for k = 1:length(window)

        %% Extract window data
        nn = window(k);
        CPPopt_loop = [];
        window_CPP = CPP_av(i-nn:i);
        window_PRx = PRx_av(i-nn:i);

        %% Binning
        % Local minima and maxima for bin edges
        min_loop = floor(min(window_CPP)) - mod((floor(min(window_CPP))),5); 
        max_loop = ceil(max(window_CPP)) + mod((ceil(max(window_CPP))),5); 

        % Determine bin edges
        bin_edges_loop = min_loop:bin_width:max_loop;
        bin_indices_loop = cell(numel(bin_edges_loop)-1,1);

        % Divide measured CPP values into bins 
        for n = 1:length(window_CPP)
            CPP_value = window_CPP(n);

            for j = 1:numel(bin_edges_loop)-1
                if CPP_value >= bin_edges_loop(j) && CPP_value < bin_edges_loop (j+1)
                    bin_indices_loop{j} = [bin_indices_loop{j},n];
                end
            end 
        end 
        
        %% Calculations over bins
        % Calculate mean and standard deviation of PRx for each bin 
        bin_means_loop = NaN(size(bin_edges_loop));
        bin_std_loop = NaN(size(bin_edges_loop));

        for f = 1:numel(bin_edges_loop)-1
            indices_loop = bin_indices_loop{f};
            bin_data_loop = window_PRx(indices_loop);
            bin_means_loop(f) = mean(bin_data_loop,'omitnan');
            bin_std_loop(f) = std(bin_data_loop,'omitnan');
        end

        % Find bins that contain at least 1% of measured data
        valid_bin_means_loop = bin_means_loop;

        for n = 1:length(bin_indices_loop)
            if numel(bin_indices_loop{n}) <= 0.01*numel(window_CPP)
                valid_bin_means_loop(n) = NaN;
            end
        end

        % Only use >= 1% bins for curve fitting to remove outliers
        valid_indices_loop = ~isnan(valid_bin_means_loop);
        valid_bin_edges_loop = bin_edges_loop(valid_indices_loop);
        valid_bin_means_loop = valid_bin_means_loop(valid_indices_loop);

        %% Curve fitting
        % Fit second order polynomial curve over remaining bins 
        coefficients_loop = polyfit(valid_bin_edges_loop,valid_bin_means_loop',2);
        fitted_curve_loop = polyval(coefficients_loop,valid_bin_edges_loop);
        fitted_curve_func_loop = @(x) polyval(coefficients_loop,x);
        
        %% Calculate CPPopt in window
        % Find local minimum and maximum, where CPPopt is the local minimum
        min_curve = min(valid_bin_edges_loop);
        max_curve = max(valid_bin_edges_loop);

        % If the minimum corresponds with PRx > 0.2 the CPPopt is replaced
        % by NaN due to impaired autoregulation
        if ~isempty (valid_bin_edges_loop)
            if ~isnan(min_curve) && ~isnan(max_curve) 
                min_x_loop{k} = fminbnd(fitted_curve_func_loop,min_curve,max_curve);
                min_y_loop{k} = polyval(coefficients_loop,min_x_loop{k});
            else
                min_x_loop{k} = NaN;
                min_y_loop{k} = NaN;
            end
            
            if min_y_loop{k} < 0.2
                CPPopt_loop = [CPPopt_loop; min_x_loop; min_y_loop];        
            else 
                CPPopt_loop = [CPPopt_loop; NaN; NaN];
            end              
        end

        %% Calculate optimal CPP range
        % Calculate optimal range as min PRx + 0.2 in the 8 hour window,
        % allowing up to PRx = 0.25 as threshold for intact autoregulation.
        if nn == 480                                                        % 8 hour window in minutes       
            if min(valid_bin_means_loop) <= 0.05                           
                ll_prx = min(valid_bin_means_loop);
                ul_prx = ll_prx + 0.2;
            else 
                ll_prx = NaN;
                ul_prx = NaN;
            end

            if  ~isnan(ll_prx) &&  ~isnan(ul_prx)
                range_cpp{i,1} = [valid_bin_edges_loop(valid_bin_means_loop <= ul_prx)]';
            else
                range_cpp{i,1} = NaN;
            end

        end
    end
    %% Determine overall CPPopt
    % Extract resultant PRx and CPP of windows
    if isempty(CPPopt_loop) 
        CPPopt_loop = NaN;
    elseif iscell(CPPopt_loop)
        loop_result_cpp = cell2mat(CPPopt_loop(1,:));
        loop_result_prx = cell2mat(CPPopt_loop(2,:));
    else 
        loop_result_cpp = CPPopt_loop(1,:);
        loop_result_prx = CPPopt_loop(2,:); 
    end 
    
    % Calculate mean PRx
    mean_loop_result_prx = nanmean(loop_result_prx);
    
    % Find lowest PRx (best autoregulation) to find 'best' CPP
    index_lowest_prx = find(loop_result_prx == min(loop_result_prx), 1);
    cpp_with_lowest_prx = loop_result_cpp(index_lowest_prx);

    % The CPPopt that will be displayed is the weighted average between the
    % mean CPP overall and the 'best CPP'. However, CPPopt values with PRx
    % < 0.2 i.e. impaired autoregulation are rejected.
    weighted_av = nanmean([nanmean(loop_result_cpp), cpp_with_lowest_prx]);

    if mean_loop_result_prx < 0.2 && ~isnan(weighted_av)
        loop_display = nanmean([nanmean(loop_result_cpp), cpp_with_lowest_prx]);
    else
        % Here you can enter static targets as fail-safe mechanisms, e.g.
        % age-based targets which may be used in the case that dynamic
        % targets cannot be determined due to impaired autoregulation. 
    end

    % Trend of CPPopt with rejection of impaired autoregulation
    CPPopt_trend_display = [CPPopt_trend_display loop_display]; 
    
    % Trend of weighted average CPPopt over all windows
    CPPopt_trend = [CPPopt_trend weighted_av]; 

end

%% Save to output
CPPopt.trend = CPPopt_trend;
CPPopt.trend_display = CPPopt_trend_display;
CPPopt.range = range_cpp;

%% Visualize result
% t10 = t(1:10:length(t))                                                     % Downsampled time axis of input CPP data (for visualization only)
% start_calc = 8*60*6;                                                        % Time calculations in seconds (8 h window)
% t_av_indices = start_calc + (0:length(CPPopt_trend)-1)*(60/Ts);             % Indices of the CPPopt trend (with 1 min increments)
% t_av = t10(t_av_indices);                                                   % Get corresponding timestamps
% 
% figure, hold on ;
% plot(t10,CPP,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
% plot(t_av,CPPopt_trend,'Linewidth',2,'Color','blue');
% title ('CPPopt calculated using 5 time windows','FontSize',18)
% xlabel ('Time','FontSize',14)
% ylabel ('Pressure [mmHg]','FontSize',14)
% legend ('Measured CPP','CPPopt','FontSize',12)

end
