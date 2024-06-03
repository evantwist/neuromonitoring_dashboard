%% Function for analysis of ICP on an insult-basis
function [ICP_insults] = func_ICP(ICP, outcome)
% Input:
%       ICP_data: ICP time series (optionally filtered to remove noise)
%       outcome_data: to make the correlation plot, individual ICP data
%       needs to be correlated to individual outcome scores (e.g. PCPC,
%       GCS) at predefined follow-up.
% 
% Ideally use a mat-file or struct containing ICP and outcome data of all
% patients, or refer to a directory (e.g. files = dir('*.mat')) to perform
% a loop. Here it is illustrated for a single patient-file as an example.
%
% Output:
%       Color-coded figure displaying the correlation between the average
%       number of ICP insults (denoted by their intensity in mmHg and
%       duration in min) and outcome (denoted by outcome score, e.g. PCPC).
%       
% Inspired by Guiza et al. Adjusted from Tahisa Robles
% Erasmus MC Sophia Childrens Hospital, Pediatric Intensive Care Unit in
% collaboration with Delft University of Technology, 3mE 
%
% Disclaimer: not validated for real time clinical use
% Eris van Twist, May 2024

clear all
close all
clc 

ICP_insults = struct();

%% Determine ICP insults per patient
%files = refer to a directory/struct/.mat-file
Ts = 10;      % Seconds per sample

for i = 1:length(files)
    
    % Preallocate for speed [360 min x 60 mmHg range]
    intensity = array2table(zeros(360,60));                                  

    % Load file
    load(files(i).name);

    % Find ICP insults from 5 min duration and 10 mmHg intensity
    for row = 5:1:360
        for col = 10:1:60
            % Duration threshold in samples
            threshold_t  = row*60/Ts;
            
            % Intensity threshold
            threshold_p = ICP.filt >= col;
            episode_start = find(diff([0; threshold_p]) == 1);
            episode_end = find(diff([threshold_p; 0]) ==-1);
            
            % Exceed threshold duration
            valid_insults = (episode_end - episode_start +1)>= threshold_t;
    
            % Number of insults where ICP >= threshold for at least >= duration threshold
            intensity(row, col) = table(sum(valid_insults));
        end
    end 

    % Save table for each patient 
    %filename = fullfile(outputdir,['intensity_' files(i).name]);
    %save(filename,'intensity')
    insults{i} = intensity;
end 


%% Average ICP insults across outcome groups
% This section depends on the outcome data available and the scale used to
% score outcome. In our case using the PCPC score it may look something
% like the example below.

% outputdir = refer to directory or file where outcomes are stored
% outcome_data = refer to directory containing outcomes of patients

% Loop over outcome score (PCPC scores)
for i = 1:6 
    % Fetch patients in that outcome group
    pids = outcomes.PID(outcomes.outcome == i);

    % Preallocate matrix
    IntensityData = cell(1, length(pids));

    % Load all intensity tables for the current outcome group 
    for j = 1:length(pids)
        fileName = strcat('intensity',num2str(pids(j)));
        intensity = load(fileName);
        IntensityData{j} = intensity;
    end 
    
    % Preallocate table
    MeanTable = zeros(360,60);

    for idx = 1:numel(IntensityData)
        currentTable = table2array(IntensityData{idx}.intensity);
        MeanTable = MeanTable + currentTable; 
    end 

    MeanTable = MeanTable/(numel(IntensityData));
  
    % Save mean average number of ICP insults per outcome group
    filename = fullfile(outputFolder,['MeanIntensityOutcome_' num2str(i)]);
    save(filename,'MeanTable')

end 


%% Calculate correlation between average ICP insults and outcome

correlations = table();

% PCPC scores (5 was not in cohort)
scores = [1, 2, 3, 4 ,6];

% Loop over insult duration [min] and intensities [mmHg]
for row = 1:360
    for col = 1:60
        % Average number of ICP insults
        x = [MeanIntensityOutcome_1(row,col), MeanIntensityOutcome_2(row,col), ...
            MeanIntensityOutcome_3(row,col), MeanIntensityOutcome_4(row,col), ...
            MeanIntensityOutcome_6(row,col)];
        % If insults are not present, return NaN
        if sum (x,'omitnan') == 0 
            correlations(row,col) = table(NaN);
        % Otherwise calculate correlation coefficient
        elseif sum (x,'omitnan') > 0          
            correlation = corrcoef(x, scores);
            correlations(row,col) = table(correlation(1,2));
        end 
    end 
end

%% Visualize results
% Folder publicatie > resultaten intensity tabellen 

rownames = arrayfun(@(x) num2str(x),1:size(correlations,1),'UniformOutput',false);
colnames = arrayfun(@(x) num2str(x),1:size(correlations,2),'UniformOutput',false);
correlations.Properties.RowNames = rownames;
correlations.Properties.VariableNames = colnames;

% Remove ICP 40-60 mmHg (non-physiological)
correlations(:,40:60) = [];

% Remove empy rows and columns (ICP < 10 mmHg or duration < 5 min)
correlations(1:4,:) = []; 
correlations(:,1:9) = [];

correlations = correlations{:,:};


figure;
colormap(jet); 
flippedCorrelations = flipud(correlations);
imagesc(flippedCorrelations);
caxis([min(min(correlations(~isnan(correlations)))),max(max(correlations(~isnan(correlations))))]);
colorbar;
xlabel('Insult intensity [mmHg]')
ylabel('Insult duration [minutes]')
xticks(1:5:40);
xticklabels(10:5:40);
yticks(1:30:360)
yticklabels(linspace(360,30,12))
title('')
set(gca,'FontSize',12)

%%

% %% Calculate average number of ICP insults
% % files = refer to a directory/struct/.mat-file
% bin_edges = 1:1:100;
% 
% % Initialize count vectors
% total_counts_g = table();
% total_counts_p = table();
% total_counts_m = table();
% 
% % Initialize count vectors normalized for time 
% counts_g_time = table();
% counts_p_time = table();
% counts_m_time = table();
% 
% % Calculate length of ICP monitoring period 
% final_day = arrayfun(@(x) find(~isnan(ICP.filt(:,x)), 1, 'last'), 1:size(ICP.filt, 2));     % Non-nan final index of ICP
% length_hosp = days(ICP.filt(final_day) - ICP.filt(1));                         
%     
% % Good outcome
% if outcome < 4
%     [counts_g, ~] = histcounts(ICP, bin_edges);
%     total_counts_g = [total_counts_g; array2table(counts_g)];
%     time_norm = counts_g./length_hosp;
%     counts_g_time = [counts_g_time; array2table(time_norm)];
% end 
% 
% % Poor outcome
% if outcome > 3 && outcome < 6
%     [counts_p,~] = histcounts(ICP, bin_edges);
%     newrow_p = array2table(counts_p);
%     time_norm_p = newrow_p./length_hosp;
%     total_counts_p = [total_counts_p;newrow_p];
%     counts_p_time=[counts_p_time;time_norm_p];
% end 
% 
% % Mortality
% if outcome == 6 
%     [counts_m,~] = histcounts(ICP, bin_edges);
%     newrow_m = array2table(counts_m);
%     total_counts_m = [total_counts_m;newrow_m];
%     time_norm_m = newrow_m./length_hosp;
%     counts_m_time = [counts_m_time;time_norm_m];
% end 
% 
% 
% % Calculate mean count per minute for group
% counts_g_av = (mean(total_counts_g{:,:},1))/60; 
% counts_p_av = (mean(total_counts_p{:,:},1))/60; 
% counts_m_av = (mean(total_counts_m{:,:},1))/60; 
% 
% % Normalized for time
% counts_g_time_av = (mean(counts_g_time{:,:},1))/60; 
% counts_p_time_av = (mean(counts_p_time{:,:},1))/60; 
% counts_m_time_av = (mean(counts_m_time{:,:},1))/60; 
% 
% %% Visualize results
% 
% % Non-normalized counts
% figure, hold on;
% plot(bin_edges(1:end-1),counts_g_av,'Color','green','LineWidth',2)
% plot(bin_edges(1:end-1),counts_p_av,'Color',[1,0.5,0],'LineWidth',2)
% plot(bin_edges(1:end-1),counts_m_av,'Color','red','LineWidth',2)
% legend ('Good outcome (PCPC 1-2)','Poor Outcome (PCPC 3-5)','Mortality (PCPC 6)')
% title ('ICP distribution without time normalization ')
% xlabel ('ICP')
% ylabel ('Time [minutes]')
% 
% 
% % Counts normalized for time
% figure, hold on;
% plot(bin_edges(1:end-1),counts_g_time_av,'Color','green','LineWidth',2)
% plot(bin_edges(1:end-1),counts_p_time_av,'Color',[1,0.5,0],'LineWidth',2)
% plot(bin_edges(1:end-1),counts_m_time_av,'Color','red','LineWidth',2)
% legend ('Good outcome (PCPC 1-2)','Poor Outcome (PCPC 3-5)','Mortality (PCPC 6)')
% title ('ICP distribution with time normalization ')
% xlabel ('ICP')
% ylabel ('Time [minutes]')
% 
% 
% %% Segment ICP data into ICP insults 
% intensityThresholds=10:40;
% durationThresholds=30:2160;
% 
% icpInsults=cell(length(intensityThresholds),length(durationThresholds));
% 
% for i = 1:length(intensityThresholds)
%     for j = 1:length(durationThresholds)
%         intensityThreshold=intensityThresholds(i);
%         durationThreshold=durationThresholds(j);
% 
%         currentInsult=[];
%         insultStartTime=0;
% 
%         %Loop through ICP data
%         for k=1:length(ICP.filt)
%             icp=ICP.filt(k);
%             t=time.t10;
% 
%             if icp > intensityThreshold
%                 if isempty(currentInsult)
%                     insultStartTime=t;
%                 end 
% 
%                 currentInsult = [currentInsult,icp];
%             else
%                 if ~isempty(currentInsult) && all(t - insultStartTime >= durationThreshold)
%                     icpInsults{i,j}=[icpInsults{i,j};currentInsult];
%                 end 
%                 currentInsult =[];
%             end 
%         end 
% 
%         if ~isempty (currentInsult) && all(t - insultStartTime >= durationThreshold)
%             icpInsults{i,j}=[icpInsults{i,j};currentInsult];
%         end 
%     end 
% end 
% end