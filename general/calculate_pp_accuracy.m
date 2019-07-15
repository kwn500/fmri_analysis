%% ------------- CALCULATE PARTICIPANT RESPONSE ACCURACY --------------- %%

% Loads in participant trial-by-trial response accuracies (hits, misses,
% false alarms and correct rejections), and calculates percentage hits,
% overall percentage correct (hits and correct rejections) and d prime for
% each run, then averages across runs to give one value for each of these
% for each participants. Outputs data to a .csv file.

% KWN 14/06/2018

clear; close; % general housekeeping.
addpath(genpath('/scratch/groups/Projects/P1323/colour/Code/functions')); % add the function directory.

%%  ---------------- SPECIFY PARTICIPANT INFORMATION ------------------- %%

% participants{1} = 'R2268';
participants{1} = 'R2590';
participants{2} = 'R2904';
participants{3} = 'R3111';
participants{4} = 'R3517';
participants{5} = 'R3773';
participants{6} = 'R4059';
participants{7} = 'R4065';
participants{8} = 'R4127';
% participants{10} = 'R4200';
participants{9} = 'R4244';
participants{10} = 'R4831';
participants{11} = 'R4928';

%% -------------------- CALCULATE RESPONSE ACCURACY -------------------- %%

for i = 1:length(participants) % for each participant in turn:
    
    participant = participants{i}; % specify the individual participant to process.
    
    % get the total number of runs that participant completed.
    runs = dir(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/logs/run*', participant));
    
    for x = 1:length({runs.name}) % for each of these runs in turn:
        
        % get the name of the specific run directory.
        all_runs = {runs.name}; run = all_runs{x};
        
        % load in the participant and run-specific edited response .csv file.
        file = dir(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/logs/%s/*%s_%s_data_edit.csv', participant, run, participant, run));
        
        % load in the data columns of interest.
        [feature, marked_response] = csvimport(strcat(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/logs/%s/', participant, run),...
            file.name), 'columns', {'feature','marked_response'});
        
        data = [feature, marked_response]; % combine the data columns into a singular array. 
        
        % remove data from any trial with a passive attentional focus. 
        passive_remove = strcmp(data(:,1),'passive');
        data(passive_remove,:) = []; 
        
        % remove data from any trial in which the participant did not respond. 
        noresp_remove = strcmp(data(:,2),'no_response');
        data(noresp_remove,:) = [];
        
        % count the number of trials in which each type of response was
        % made (hit, miss, false alarm and correct rejection). 
        hit = strcmp(data(:,2), 'Hit'); hit = sum(hit);
        miss = strcmp(data(:,2), 'Miss'); miss = sum(miss);
        fa = strcmp(data(:,2), 'False_Alarm'); fa = sum(fa);
        cr = strcmp(data(:,2), 'Correct_Rejection'); cr = sum(cr);
        
        signal = hit+miss; % calculate the total number of signal trials. 
        noise = fa+cr; % calculate the total number of noise trials. 
        
        percent_hit = (hit/signal)*100; % calculate the percentage of hits participants made. 
        percent_fa = (fa/noise)*100; % calculate the percentage of false alarms participants made. 
        
        % to calculate loglinear percentage hits and false alarms, we add
        % 0.5 to the number of hits and false alarms, and 1 to the number
        % of signal and noise trials (prevents infinite d prime values). 
        loglinear_percent_hit = (hit+0.5)/(signal+1);
        loglinear_percent_fa = (fa+0.5)/(noise+1);
        
        % calculate d prime. 
        dp = norminv(loglinear_percent_hit) - norminv(loglinear_percent_fa);
        
        % calculate the percentage of overall correct responses (hits and
        % correct rejections). 
        percent_correct = ((hit+cr)/(signal+noise))*100;
        
        % add these run-specific values to storage arrays. 
        total_hit(x) = percent_hit;
        total_dp(x) = dp;
        total_correct(x) = percent_correct;
    end % repeat for the next run. 
    
    % average across all runs for a singular participant and store these 
    % mean values in participant-specific locations in output storage arrays. 
    pp_hits(i) = mean(total_hit);
    pp_dp(i) = mean(total_dp);
    pp_correct(i) = mean(total_correct);
    
    clear total_hit total_dp total_correct % refresh some variables for the next loop iteration. 
end

% concatenate this data from each participant into a singular array.
data = [participants; num2cell(pp_hits); num2cell(pp_dp); num2cell(pp_correct)]';

% convert this cell array to a data table with relevant headings.
data_table = cell2table(data, 'VariableNames', {'participant','percenthit','dp','percentcorrect'});

% save this data to a .csv file.
writetable(data_table,'/scratch/groups/Projects/P1323/colour/fMRI/general_output/response_accuracy.csv');

%% -------------------------------------------------------------------- %%