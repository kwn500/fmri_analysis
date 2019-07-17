%% ------- PROCESS NATURALISTIC EXPERIMENT PSYCHOPY RESPONSE DATA ------ %% 

% reads in participant-specific response files, and we create FSL
% parfiles dictating when participants were 'attending' (keypress) versus
% non attending (no keypress). 

% 17/7/19 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %% 

clear; clc; % general housekeeping. 

participants = {'R2548', 'R2904', 'R3111', 'R3517', 'R3773', 'R4059', 'R4065', 'R4127',...
    'R4244', 'R4829', 'R4831', 'R4833', 'R4890', 'R4928', 'R5006'}; % specify participants to analyse.

%% ------------------------- CREATE EVENT FILES ------------------------ %%

for thispp = 1:length(participants); % for each participant in turn:
    participant = participants{thispp};
    
    % extract runs collected for this participant.
    runs = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/psychopy_data/runs/run_*', participant));
    runs = {runs.name};
    
    for thisRun = 1:length(runs); % for each run in turn:
        run_num = runs{thisRun};
        
        % load in event .txt file.
        data = dlmread(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/psychopy_data/runs/%s/%s_log_%s.txt',...
            participant, run_num, participant, run_num(end)));
        
        % extract frame and keypress data individually.
        frames = data(:,1);
        keypresses = data(:,2);
        
        % rescale the frame timing, to start at 0, as currently starts at 8
        % seconds (accounting for the four dummy volumes we included). 
        frames = frames- frames(1);
        
        % identify the frames in which a keypress is evident.
        keypress_index = find(keypresses == 1)';
        D = diff([0,diff(keypress_index)==1,0]);
        first = keypress_index(D>0);
        last = keypress_index(D<0);
        
        % extract the frame timings of the start and end of each keypress. 
        for thisevent = 1:length(first);
            starts(thisevent) = frames(first(thisevent));
            ends(thisevent) = frames(last(thisevent));
        end
        
        % round down the frame timings to the nearest TR (2 seconds). 
        keypress_timingsTRs = 2*floor(starts/2);
        
        % repeat the same process for the no-keypress periods.
        nokeypress_index = find(keypresses == 0)';
        D = diff([0,diff(nokeypress_index)==1,0]);
        first = nokeypress_index(D>0);
        last = nokeypress_index(D<0);
        
        for thisevent = 1:length(first);
            starts(thisevent) = frames(first(thisevent));
            ends(thisevent) = frames(last(thisevent));
        end

        nokeypress_timingsTRs = 2*floor(starts/2);

        % combine the TRs and event codes to produce mrvista parfiles (2 =
        % keypress, 1 = no keypress). 
        vista_keypress = [keypress_timingsTRs', ones(length(keypress_timingsTRs),1)*2];
        vista_nokeypress = [nokeypress_timingsTRs', ones(length(nokeypress_timingsTRs),1)];
        
        % combine and sort the keypress and no keypress data for the
        % mrvista parfile. 
        vista_data = [vista_keypress; vista_nokeypress];
        vista_data = sortrows(vista_data,1);
        
        % calculate the duration of each keypress and no keypress period.
        timing1 = vista_data(:,1); timing1(end+1) = 300;
        timing2 = vista_data(:,1); timing2 = [0;timing2];
        durations = timing1-timing2; durations(1) = [];
        
        % remove any keypress durations less than a TR (here, this is 0 as
        % the timings already been rounded to the nearest TR).
        short_timings = find(durations == 0);
        vista_data(short_timings,:) = [];
        durations(short_timings) = [];

        % format the keypress data for writing to FSL event files- timing
        % onsets, durations, and a column of ones. 
        fsl_data = [vista_data,durations];
        fsl_keypresses = find(fsl_data(:,2) == 2);
        fsl_keypress_data = fsl_data(fsl_keypresses,:);
        fsl_keypress_data_formatted = [fsl_keypress_data(:,1), fsl_keypress_data(:,3), ones(size(fsl_keypress_data,1),1)];
        
        % repeat the same process for the no keypress data. 
        fsl_nokeypresses = find(fsl_data(:,2) == 1);
        fsl_nokeypress_data = fsl_data(fsl_nokeypresses,:);
        fsl_nokeypress_data_formatted = [fsl_nokeypress_data(:,1), fsl_nokeypress_data(:,3), ones(size(fsl_nokeypress_data,1),1)];
        
        % save the keypress and no keypress data to FSL .txt files.
        data_directory = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/psychopy_data/processed/', participant);
        [~,~] = mkdir(data_directory);
        fsl_directory = strcat(data_directory, 'FSL'); [~,~] = mkdir(fsl_directory);
        vista_directory = strcat(data_directory, 'vista'); [~,~] = mkdir(vista_directory);
        
        dlmwrite(strcat(fsl_directory, sprintf('/%s_run_%s_keypress_log.txt', participant, run_num(end))),...
            fsl_keypress_data_formatted,'delimiter','\t','precision','%.2f');
        dlmwrite(strcat(fsl_directory, sprintf('/%s_run_%s_nokeypress_log.txt', participant, run_num(end))),...
            fsl_nokeypress_data_formatted,'delimiter','\t','precision','%.2f');
        
        % save the combined keypress and no keypress data to a singular
        % mrvista parfile.
        dlmwrite(strcat(vista_directory, sprintf('/%s_run_%s_response_log.txt', participant, run_num(end))),...
            vista_data, 'delimiter', '\t', 'precision', '%.2f');
    end % continue to next run.
end % continue to next participant. 

%% --------------------------------------------------------------------- %%