%% ------- PROCESS NATURALISTIC EXPERIMENT PSYCHOPY RESPONSE DATA ------ %% 

% reads in participant-specific response files, and we create FSL
% parfiles dictating when participants were 'attending' (keypress) versus
% non attending (no keypress). 

% 16/7/19 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %% 

clear; clc; % general housekeeping. 

Rnumber = 3773; % specify participant r number. 
Run = 'A'; % specify run to analyse. 

% for now, index into data storage directory (needs updating/not indexing). 
cd('/home/k/kwn500/Desktop'); 

%% ------------------------- CREATE EVENT FILES ------------------------ %%

% load in event .txt file.
data = dlmread(sprintf('R%d_Run_%d_log.txt', Rnumber, Run));

% extract frame and keypress data individually.
frames = data(:,1); 
keypresses = data(:,2);

% identify the frames in which a keypress is evident.
keypress_index = find(keypresses == 1);

for i = 1:length(keypress_index) % for each keypress frame.
    
    % if this is the first instance of the keypress event, extract the
    % frame number. 
    if i == 1 || keypress_index(i)-1 ~= keypress_index(i-1) 
        keypress_starts(i) = frames(keypress_index(i));
    end
    % else if this is the last instance of the keypress event, extract
    % the frame number. 
    if i == length(keypress_index) || keypress_index(i) +1 ~= keypress_index(i+1)
        keypress_ends(i) = frames(keypress_index(i));
    end
end

% remove the empty cells from the keypress start and end arrays. 
keypress_starts(keypress_starts==0) = []; keypress_starts = keypress_starts-frames(1);
keypress_ends(keypress_ends==0) = []; keypress_ends = keypress_ends-frames(1);

% calculate the duration of the keypresses. 
keypress_durations = keypress_ends-keypress_starts;

% convert to FSL format (rounding the timings to the nearest TR). 
fsl_keypress_data = [round(keypress_starts*2)/2; round(keypress_durations*2)/2; ones(1,length(keypress_starts))]';

% repeat the same process for the no-keypress periods. 
nokeypress_index = find(keypresses == 0);
for i = 1:length(nokeypress_index)
    if i == 1 || nokeypress_index(i)-1 ~= nokeypress_index(i-1)
        nokeypress_starts(i) = frames(nokeypress_index(i));
    end
    if i == length(nokeypress_index) || nokeypress_index(i) + 1 ~= nokeypress_index(i+1)
        nokeypress_ends(i) = frames(nokeypress_index(i));
    end
end

nokeypress_starts(nokeypress_starts==0)=[]; nokeypress_starts = nokeypress_starts-frames(1);
nokeypress_ends(nokeypress_ends==0)=[]; nokeypress_ends = nokeypress_ends-frames(1);
nokeypress_durations = nokeypress_ends-nokeypress_starts;
fsl_nokeypress_data = [round(nokeypress_starts*2)/2; round(nokeypress_durations*2)/2; ones(1,length(nokeypress_starts))]';

% save the keypress and no keypress data to FSL .txt files. 
data_directory = sprintf('EventFiles_R%d', Rnumber); [~,~] = mkdir(data_directory); cd(data_directory); 
fsl_directory = ('FSL'); [~,~] = mkdir(fsl_directory); cd(fsl_directory); 

dlmwrite(sprintf('R%d_Run_0%d_keypress_log.txt', Rnumber, Run),fsl_keypress_data,'delimiter','\t','precision','%.2f');
dlmwrite(sprintf('R%d_Run_0%d_nokeypress_log.txt', Rnumber, Run),fsl_nokeypress_data,'delimiter','\t','precision','%.2f');

%% --------------------------------------------------------------------- %%