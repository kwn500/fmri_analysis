%% ------------------- CALCULATE FMRI DPRIME SCORES -------------------- %%

% loads in participant- and run- specific .mat files from original
% achromatic experiment. calculates loglinear d prime for analysis with
% the rmse distance scores from connectivity analysis. 

% KWN 9/7/2019

clear;clc;

%% ---------------------- CALCULATE DPRIME SCORES ---------------------- %%

% specify participants to analyse.
participants = {'R2268', 'R2548', 'R2590','R2904', 'R3111', 'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

for thispp = 1:length(participants) % for each participant:
    participant = participants{thispp}; % extract the participant r number. 
    
    % list the run-specific files. 
    files = dir(sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/logs/performance/*.mat', participant));
    filenames = {files.name};
    
    for thisrun = 1:length(filenames) % for each run:
        
        % load the run-specific data file. 
        load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/logs/performance/%s', participant, filenames{thisrun}));
        
        % extract the number of specific responses (minus one as these have
        % been incremented by one in our old experimental code). 
        hit = hit_counter-1;
        miss = miss_counter-1;
        fa = FA_counter-1;
        cr = CR_counter-1;
        nr = no_response_counter-1;
        
        [signal, noise] = deal(0); % initalise signal and noise counters. 
        
        for thistrial = 1:length(changes) % for each trial:
            
            % if stimulus changes matches attention condition:
            if strfind(change_type{thistrial},char(condition{thistrial})) >= 1
                signal = signal + 1; % record as a signal trial.
            else
                noise = noise + 1; % otherwise record as a noise trial. 
            end
        end
           
        % increment signal and noise by one and hit and false alarm by 0.5
        % for loglinear dprime calculation. 
        signal = signal + 1;
        noise = noise + 1;
        hit = hit + 0.5;
        fa = fa + 0.5;
        
        % convert signal and noise to percentages. 
        signalpercent = (hit/signal)*100;
        noisepercent = (fa/noise)*100;
        
        % calculate loglinear dprime. 
        logdprime = norminv(signalpercent/100)-norminv(noisepercent/100);
        
        logdprime_runs(thisrun) = logdprime; % store run-specific loglinear dprime value. 
    end
    
    % calculate and store the mean loglinear dprime across runs. 
    logdprime_participants(thispp) = mean(logdprime_runs);
    clear logdprime_runs % refresh for next iteration. 
end

save('logdprime_participants.mat', 'logdprime_participants'); % save data in .mat file. 

%% --------------- CORRELATE DPRIME WITH CONNECTIVITY DATA ------------- %%

% load in correlation coefficients from connectivity analysis and the fMRI
% loglinear dprime scores.
load('/Users/kirstie/Documents/analysis/original/analysis_output/connectivity-fisher.mat');
load('/Users/kirstie/Documents/analysis/original/analysis_output/logdprime_participants.mat');

% extract the fisher-transformed correlation coefficient data.
data = connectivity.indconditioncorrs.fisherdata;

ROI = 'V1'; % specify ROI to analyse. 

for i = 1:3 % for each of the attention conditions:
    for thispp = 1:size(data,2) % for each participant in turn:
        
        % calculate the rmse distance between the attention-specific data
        % across partial correlations with the ROI of interest and the same
        % data for the passive viewing condition. 
        if isequal(ROI, 'V1')
            rmsedistance(i,thispp,:) = sum((squeeze(data(i,thispp,2:end,1))' - squeeze(data(4,thispp,2:end,1))').^2);
        elseif isequal(ROI, 'V4')
            rmsedistance(i,thispp,:) = sum((squeeze(data(i,thispp,3,[1,2,4:6]))' - squeeze(data(4,thispp,3,[1,2,4:6]))').^2);
        elseif isequal(ROI, 'LO1')
            rmsedistance(i,thispp,:) = sum((squeeze(data(i,thispp,4,[1,2,3,5,6]))' - squeeze(data(4,thispp,4,[1,2,3,5,6]))').^2);
        elseif isequal(ROI, 'LO2')
            rmsedistance(i,thispp,:) = sum((squeeze(data(i,thispp,5,[1,2,3,4,6]))' - squeeze(data(4,thispp,5,[1,2,3,4,6]))').^2);
        elseif isequal(ROI, 'IPS0')
            rmsedistance(i,thispp,:) = sum((squeeze(data(i,thispp,6,1:5))' - squeeze(data(4,thispp,6,1:5))').^2);
        elseif isequal(ROI, 'all')
            rmsedistance(i,thispp,:) = sum((squeeze(data(i,thispp,:,:))' - squeeze(data(4,thispp,:,:))').^2);
        end
    end
    
    % calculate the correlation between these rmse distances and the
    % loglinear dprime score for this attention condition.
    [r(i), p(i)] = corr(rmsedistance(i,:)', logdprime_participants');
end

% save the resulting correlation and p-value matrices. 
save('logdprime_participants_correlation.mat', 'r', 'p'); % save data in .mat file. 

%% --------------------------------------------------------------------- %%