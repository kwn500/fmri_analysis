%% --- FORMAT OUTPUT OF FACE-DECODER NEURAL NETWORK TO FSL EVENT FILE -- %% 

% Loads in the output data from the YOLO neural net identifying the number
% of faces in each frame of the visual stimulus (sampled at 2 frames per
% second). we create FSL-formatted event files, with a face TR specified as
% a TR containing an average of less than 2 faces, and at least 300
% face-specific pixels. 

% 25/7/2019 KWN

clear; clc;

%% ------------------------- CREATE FSL FILES -------------------------- %% 

runs = {'A','B','C','D','E','F','G'}; % specify runs to analyse. 

for thisRun = 1:length(runs); % for each run in turn:
    run = runs{thisRun};
    
    % load in the run-specific neural network output file (excluding
    % headers and frame-specific file names). we also exclude data from the
    % first and last frame sampled.
    C = csvread(sprintf('/home/k/kwn500/Desktop/face_neuralnet/%s_output.csv',...
        run),2,1,'B3..D602');
    
    % for each TR, we calculate the mean number of faces identified and the
    % mean number of face-specific pixels. 
    count = 1;
    rowcount = 1;
    for thisTR = 1:length(C);
        if thisTR >= 151;
            continue
        end
        face_data(rowcount) = mean([C(count,1), C(count+1,1), C(count+2,1), C(count+3,1)]);
        pixel_data(rowcount) = mean([C(count,2), C(count+1,2), C(count+2,2), C(count+3,2)]);
        count = count + 4;
        rowcount = rowcount + 1;
    end
    
    % for each TR, we identify a 'face' TR as one containing two or less
    % faces (mean) and more than an average of 300 face-specific pixels. 
    for thisRow = 1:length(face_data);
        if face_data(thisRow) <= 2 && pixel_data(thisRow) >= 300;
            event_code(thisRow)= 1;
        else
            event_code(thisRow) = 0;
        end
    end

    % create a matrix corresponding to an FSL-formatted event files, with
    % TR, event code, duration, and a mean column (1s). 
    format_data = [[0:2:299]', event_code', ones(length(event_code),1) * 2, ones(length(event_code),1)];
    
    % extract the face-specific data and remove the event code column. 
    face_data_ind = format_data(:,2) == 1;
    face_data = format_data(face_data_ind,:); face_data(:,2) = [];
    
    % repeat for the no-face data. 
    noface_data_ind = format_data(:,2) == 0;
    noface_data = format_data(noface_data_ind,:); noface_data(:,2) = [];
    
    % count the total number of face and no-face events. 
    face_count(thisRun) = length(face_data);
    noface_count(thisRun) = length(noface_data);
    
    % write this data to two run-specific .txt files. 
    dlmwrite(sprintf('/home/k/kwn500/Desktop/face_neuralnet/%s_processed_pixelmean-face.txt',run),...
        face_data,'delimiter', '\t', 'precision', '%.2f')
    dlmwrite(sprintf('/home/k/kwn500/Desktop/face_neuralnet/%s_processed_pixelmean-noface.txt',run),...
        noface_data,'delimiter', '\t', 'precision', '%.2f')
end

% save the face/no-face count data. 
save('/home/k/kwn500/Desktop/face_neuralnet/face_count_data.mat', 'face_count', 'noface_count');

%% --------------------------------------------------------------------- %%