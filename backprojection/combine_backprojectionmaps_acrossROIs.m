%% COMBINE CONDITION-SPECIFIC SUPPORT VECTOR PARAMETER MAPS ACROSS ROIs  %%

% Loads in the previously-computed ROI- and condition- specific support
% vector parameter maps (e.g. a map indicating each V1 voxels' 'orientation
% preference'). We combine this information across voxels in specified ROIs
% (e.g. V1, V2 & V3) and create a 'combined' parameter map- containing the
% same information as the individual maps, combined in a single map for
% more efficient data visualisation.

% 08/08/2018 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

clear; clc; close all; % general housekeeping.

exp_condition = 'original'; % specify overall experiment we are analysing.
condition = 'block'; % specify the particular type of data to analyse.

% specify the type of parameter map we wish to plot, either 'difference'-
% where the maps are stored in a different directory, or 'default'. 
maptype = 'difference'; 

% specify the list of participant numbers corresponding to the participant
% data we wish to analyse.
if isequal(exp_condition, 'original')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

elseif isequal(exp_condition, 'colour')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3517', 'R3773', 'R4059', 'R4065', 'R4244', 'R4831', 'R4928'};
end

% specify the features we wish to process (we currently create seperate
% 'combined' parameter maps for each feature.
if isequal(condition, 'block') || isequal(condition, 'feature')
    
    if isequal(maptype, 'difference')
        % if this is a 'difference' analysis condition, specify the
        % corresponding difference parameter map names, otherwise, use the
        % default names.
        features = {'OrientationContrast', 'OrientationShape', 'ContrastShape'};
    else
        features = {'Orientation', 'Contrast', 'Shape'};
    end
elseif isequal(condition, 'colour')
    if isequal(maptype, 'difference')
        features = {'RedGreen-BlueYellow', 'RedGreen-Luminance', 'BlueYellow-Luminance'};
    else 
        features = {'Red-Green', 'Blue-Yellow', 'Luminance'};
    end
elseif isequal(condition, 'colourxfeature')
    if isequal(maptype, 'difference')
        features = {'BYOrientationContrast', 'BYOrientationShape', 'BYContrastShape'};
    else
        features = {'BlueYellow-Orientation', 'BlueYellow-Contrast', 'BlueYellow-Shape'};
    end
end

% specify the ROIs we wish to combine within a singular parameter map.
ROIscombine = {'V1', 'V2', 'V3'};

%% ---------------------- COMBINE PARAMETER MAPS ----------------------- %%

for featurenum = 1:length(features) % for each feature in turn:
    
    feature = features{featurenum}; % extract the feature-name we wish to analyse.
    
    % specify a feature- and ROIs-to-combine- specific filename for writing
    % out the .mat parameter map file and naming the map for use within
    % mrvista.
    filename = sprintf('%s-%s', feature, char(ROIscombine)');
    
    for participantnum = 1:length(participants) % for each participant in turn:
        
        participant = participants{participantnum}; % extract the participant-specific number.
        
        % specify the participant-specific parameter map storage directory.
        if isequal(maptype, 'difference')
            parameter_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/supportvector/difference/', exp_condition, participant, condition);
        else
            parameter_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/supportvector/', exp_condition, participant, condition);
        end
        
        for ROInum = 1:length(ROIscombine) % for each ROI in turn:
            
            ROI = ROIscombine{ROInum}; % extract the ROI-specific name.
            
            % load the condition- and ROI- specific parameter map.
            param_map = load([parameter_directory, sprintf('%s%s-%s', participant, feature, ROI)]);
            
            % save this specific parameter map to a ROI-specific location
            % in an overall storage array.
            ROImaps(ROInum) = param_map.map;
            
            % find the positions of all non-zero entries within this
            % condition- and ROI-specific parameter map (i.e. extract the
            % position of the data corresponding to the voxels within this
            % ROI).
            voxelrows = find(all((param_map.map{1} ~=0 ),1));
            
            % store this information in an ROI-specific location in a
            % storage array.
            ROIvoxels{ROInum} = voxelrows;
        end % repeat process for the next ROI.
        
        % create an array of zeros matching the size of the whole-brain
        % parameter map matrix for this specific participant.
        full_map = (1:length(ROImaps{1})).*0;
        
        for ROIvoxelnum = 1:length(ROIvoxels) % for each ROI in turn:
            
            % replace the NaNs in the matrix we have just created with the
            % data from the specific-ROI we are processing. This is
            % basically just changing NaNs to the ROI-specific data in the
            % relevant voxel locations- combining multiple parameter maps
            % into a singular matrix.
            full_map(ROIvoxels{ROIvoxelnum}) = ROImaps{ROIvoxelnum}(ROIvoxels{ROIvoxelnum});
        end % repeat process for the next ROI.
        
        % specify the parameter map aspects in the format required by
        % mrvista, specifying the map name as the filename specified avove.
        map = {full_map};
        co = param_map.co;
        mapUnits = param_map.mapUnits;
        mapName = filename;
        
        % save this new combined parameter map with the filename specified
        % above in a 'combined' directory within the support vector
        % parameter map folder.
        output_directory = strcat(parameter_directory, 'combined/');
        [~,~] = mkdir(output_directory);
        
        save(strcat(output_directory, filename), 'map', 'co', 'mapUnits', 'mapName');
    end % repeat process for the next participant.
end % repeat the process for the next feature.

%% --------------------------------------------------------------------- %%