%% ---------------------- CREATE SPM DESIGN FILES ---------------------- %%

% loads in mr-vista parfiles and restructures timing and condition
% information to produce spm-toolbox compatible design files. 

% 10/7/2019 KWN

clear; clc;

%% ---------------------- RE-FORMAT DESIGN FILES ----------------------- %%

% specify participants to analyse.
participants{1} = 'R2268';
participants{2} = 'R2548';
participants{3} = 'R2590';
participants{4} = 'R2904';
participants{5} = 'R3111';
participants{6} = 'R3455';
participants{7} = 'R3517';
participants{8} = 'R3773';
participants{9} = 'R3932';
participants{10} = 'R4065';
participants{11} = 'R4496';

for thisPP = 1:length(participants) % for each participant in turn:
    
    % identify the participant-specific parfiles.
    files = dir(sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/mrvista/block/Stimuli/Parfiles/*.par',...
        participants{thisPP}));
    parfiles = {files.name};
    
    for thisRun = 1:length(files) % for each run in turn:
        
        % load in the run-specific parfile.
        [timing, event] = textread([sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/mrvista/block/Stimuli/Parfiles/',...
            participants{thisPP}),parfiles{thisRun}], '%d %d');
        
        % extract the timings of visual stimulation (passive) and visual
        % attention (orientation, contrast and shape blocks).
        names = {'Visual Stimulation', 'Visual Attention'};
        stimind = find(event == 4);
        attnoind = find(event == 1);
        attncind = find(event == 2);
        attnsind = find(event == 3);
        
        % combine the attention data in a single array. 
        attnind = sort([attnoind; attncind; attnsind]);
        
        % extract the onsets of each condition from the TR data. 
        stimonsets = timing(stimind);
        attnonsets = timing(attnind);
        
        % specify the duration of each block (TRs). 
        durations{1} = 15;
        durations{2} = 15;
        
        % store the run-specific onset data.
        stimrunonsets{thisRun} = stimonsets;
        attnrunonsets{thisRun} = attnonsets;
    end
    
    % increment the number of TRs across scans. 
    timeadjustment = [0, 393, 393*2, 393*3, 393*4, 393*5, 393*6, 393*7];
    [stimcatonsets, attncatonsets] = deal([]);
    
    for thisRun = 1:length(stimrunonsets) % for each run in turn:
        
        % increment the onsets by the run number so that TRs are continuous
        % across all scans. 
        stimcatonsets = [stimcatonsets; stimrunonsets{thisRun}+timeadjustment(thisRun)];
        attncatonsets = [attncatonsets; attnrunonsets{thisRun}+timeadjustment(thisRun)];
    end
    
    onsets{1} = stimcatonsets;
    onsets{2} = attncatonsets;
    
    % save this edited design matrix to a .mat file.
    output_directory = sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/spm/design/',...
        participants{thisPP}); [~,~] = mkdir(output_directory);
    save(strcat(output_directory, 'stimvsattn.mat'), 'names', 'durations', 'onsets');
    
    clear names durations onsets attncatonsets stimcatonsets stimrunonsets attnrunonsets
end

%% --------------------------------------------------------------------- %%