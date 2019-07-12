%% ----- CALCULATE DIFFERENCE BETWEEN MULTIPLE SVM PARAMETER MAPS ------ %%

% Reads in two previously-created support vector parameter maps (an index
% of each voxels' 'preference' for a particular visual feature) and
% calculates the difference between these two maps, to create maps showing
% the differences in processing of two visual features in a singular ROI. 

% 03/08/2018 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

clear; clc; close all; % general housekeeping.

exp_condition = 'original'; % specify overall experiment we are analysing.
condition = 'block'; % specify the particular type of data to analyse.

% specify the list of participant numbers corresponding to the participant
% data we wish to analyse.
if isequal(exp_condition, 'original')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

elseif isequal(exp_condition, 'colour')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3517', 'R3773', 'R4059', 'R4065', 'R4244', 'R4831', 'R4928'};
end

% specify the list of ROIs we wish to analyse. 
ROIs = {'V1', 'V2', 'V3', 'V3A', 'V3B', 'V4', 'LO1', 'LO2'};

% specify the features we wish to calculate the difference between through
% multiple iterations. 
features = {'Orientation', 'Contrast', 'Shape'};

%% ---------- CALCULATE DIFFERENCE BETWEEN SVM PARAMETER MAPS ---------- %%

for ROInumber = 1:length(ROIs) % for each ROI in turn:
    
    ROI = ROIs{ROInumber}; % extract the name of the specific ROI we wish to analyse. 
    
    for participantnumber = 1:length(participants) % for each participant in turn: 
        
        participant = participants{participantnumber}; % extract the specific participant number we wish to analyse. 
        
        % if we are currently processing the MT+ ROI, and the specific
        % participant with a different high-resolution structural for the
        % motion-localiser defined MT+ ROI, specify the MT+-specific
        % session, otherwise, for all other instances, specify the
        % participant-specific default mrvista directory. 
        if isequal(participant, 'R3111') && isequal(ROI, 'MT+')
            mapdirectory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT/Gray/GLMs/supportvector/',exp_condition, participant, condition);
        else
            mapdirectory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/supportvector/',exp_condition, participant, condition);
        end
        
        % load in the three feature maps of interest from this
        % participant-specific mrvista session. 
        orientationmap = load(strcat(mapdirectory, sprintf('%s%s-%s.mat', participant, features{1}, ROI)));
        contrastmap = load(strcat(mapdirectory, sprintf('%s%s-%s.mat', participant, features{2}, ROI)));
        shapemap = load(strcat(mapdirectory, sprintf('%s%s-%s.mat', participant, features{3}, ROI)));
        
        % extract matricies of the condition-specific parameter map
        % information. 
        orientationdata = cell2mat(orientationmap.map);
        contrastdata = cell2mat(contrastmap.map);
        shapedata = cell2mat(shapemap.map);
        
        % calculate the difference in support vectors between the features 
        % of interest. 
        orientation_contrast = orientationdata - contrastdata;
        orientation_shape = orientationdata - shapedata;
        contrast_shape = contrastdata - shapedata;

        % extract the default map units and coherence thresholds from one
        % of the original maps we loaded in (this is arbitrary, we just
        % need them to save the data in a format which doesn't throw errors
        % when we load these new maps into mrvista). 
        mapUnits = orientationmap.mapUnits;
        co = orientationmap.co;
        
        % for each of the difference maps we have calculated, specify a map
        % name and save this new difference map data as a cell-array (** need
        % to loop this **). 
        mapName = 'OrientationContrast';
        map = {orientation_contrast};
        save(strcat(mapdirectory,'/difference/',(sprintf('%s%s-%s.mat', participant, mapName,ROI))), 'map', 'co', 'mapName', 'mapUnits');
        
        mapName = 'OrientationShape';
        map = {orientation_shape};
        save(strcat(mapdirectory, '/difference/',(sprintf('%s%s-%s.mat', participant, mapName,ROI))), 'map', 'co', 'mapName', 'mapUnits');
        
        mapName = 'ContrastShape';
        map = {contrast_shape};
        save(strcat(mapdirectory,'/difference/',(sprintf('%s%s-%s.mat', participant, mapName,ROI))), 'map', 'co', 'mapName', 'mapUnits');   
    end % repeat for the next participant. 
end % continue to the next ROI. 

%% --------------------------------------------------------------------- %%