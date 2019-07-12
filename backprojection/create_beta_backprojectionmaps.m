%% ----- CREATE MRVISTA PARAMETER MAPS OF CONDITION-SPECIFIC BETAS ----- %%

% Very similar to the support vector parameter script, but uses the betas
% extracted from the GLM for each condition, for every voxel (orientation,
% contrast, shape & passive). 

% 22/08/2018 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %% 

clear; clc; close all; % general housekeeping.

% add function directory to path.
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions');

% specify the participant numbers we wish to analyse (at the moment, this
% must match the size of the group beta data matrix).
participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455',...
    'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

% specify the ROI we wish to analyse and it's corresponding abbreviated
% name for naming file outputs.
ROIanalysis = '01_Combined_V1_KWN'; outputROI = 'V1';

% load the ROI-specific group extracted beta data (number of repetitions x
% number of conditions x number of voxels), with a cell entry for each
% participant. 
load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/betas/%s/group_extractedBetas_%s.mat',...
    outputROI,outputROI));

exp_condition = 'original'; % specify the overall experiment to analyse. 
condition = 'block'; % specify the analysis condition. 

% specify the names of the columns of extracted betas. 
beta_conditions = {'Orientation', 'Contrast', 'Shape', 'Passive'};

% specify the default parameter map we wish to load in (we just use this as
% a template for the total number of voxels for each participant).
default_parametermap = 'OvC';

%% -------------- PROCESS BETAS AND CREATE PARAMETER MAPS -------------- %%

for ppnum = 1:length(participants) % for each participant in turn:
    
    % extract the participant-specific number.
    participant = participants{ppnum};
    
    % if we are analysing the one participant who has a seperate
    % mrvista session for the MT+ ROI:
    if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
        
        % index into the participant-specific MT+ mrvista directory and
        % open a gray window.
        cd(sprintf('/scratch/sg3/P1323/%s/fMRI/%s/mrvista/%s/MT',...
            exp_condition, participant, condition));
        VOLUME{1} = mrVista('3view');
        
    else % otherwise, for any other ROI and participant:
        
        % index into the participant-specific mrvista directory and
        % open a gray window.
        cd(sprintf('/scratch/sg3/P1323/%s/fMRI/%s/mrvista/%s',...
            exp_condition, participant, condition));
        VOLUME{1} = mrVista('3view');
    end
    
    % load in the ROI we wish to analyse.
    VOLUME{1} = loadROI(VOLUME{1}, ROIanalysis);
    
    % refresh the view to the first scan of the GLMs data type.
    VOLUME{1} = selectDataType(VOLUME{1},4);
    VOLUME{1} = setCurScan(VOLUME{1},1);
    VOLUME{1} = refreshScreen(VOLUME{1});
    
    % if we are analysing the one participant who has a seperate
    % mrvista session for the MT+ ROI:
    if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
        
        % also load in a parameter map (one within the MT+ specific
        % directory for this participant)- this is arbitrary, we
        % just use this as a template to update with our support
        % vector data.
        load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT/Gray/GLMs/original/%s.mat',...
            exp_condition, participant, condition, default_parametermap));
        
    else % otherwise, for any other ROI and participant:
        
        % also load in a parameter map- this is arbitrary, we just use
        % this as a template to update with our support vector data.
        load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/original/%s.mat',...
            exp_condition, participant, condition, default_parametermap));
    end
    
    VOLUME{1} = refreshScreen(VOLUME{1});
    
    % extract the co-ordinates of the ROI we are currently
    % analysing (i.e. which voxels (of those across the whole
    % brain) are part of the ROI we are currently interested in).
    roipos = viewGet(VOLUME{1}, 'roiindices');
    
    % create a matrix the size of the number of voxels across the
    % participants whole brain.
    allcoords = 1:length(co{1});
    allcoords = allcoords.*0;
    
    % remove the datapoints corresponding to the positions of
    % voxels contained within the current ROI.
    allcoords(roipos) = [];
    
    % extract the whole-brain data from this current arbitrary
    % parameter map.
    mapcoords = map{1};
    
    % extract the 3D matrix of multivariate beta data corresponding to the
    % participant we are currently analysing. 
    beta_data = extractedBetasAllPs{ppnum};

    for i = 1:size(beta_data,2) % for each of the feature conditions in turn:
        
        % extract the feature-specific beta values. 
        betas = beta_data(:,i,:);
        reshape_betas = squeeze(betas);
        
        % average across the repetitions of this feature condition across
        % multiple runs. 
        meanbetas = mean(reshape_betas);
         
        % extract the beta values corresponding to the passive condition
        % and also average across these condition repetitions. 
        % passive_betas = beta_data(:,4,:);                 % uncomment for -passive 'normalised' betas.
        % passive_betas = squeeze(passive_betas);           % uncomment for -passive 'normalised' betas.
        % meanpassive = mean(passive_betas);                % uncomment for -passive 'normalised' betas.
        % meanpassive = mean(squeeze(beta_data(:,4,:)));    % uncomment for -passive 'normalised' betas.

        % extract the beta values corresponding to all conditions and take
        % the average across all conditions and all repetitions to give a
        % mean value for each voxel. 
        % meanallconds = squeeze(mean(mean(beta_data)));    % uncomment for -meanallconds 'normalised' betas. 
        
        % subtract the mean passive beta for each voxel from the mean
        % feature-specific beta value for each voxel. 
        % meanbetas = meanbetas-meanpassive;                % uncomment for -passive 'normalised' betas.

        % subtract the mean all conditions beta for each voxel from the
        % mean feature-specific beta value for each voxel. 
        % meanbetas = meanbetas-meanallconds';              % uncomment for -meanallconds 'normalised' betas. 
        
        % specify a feature-specific name for this new beta parameter map.
        mapName = sprintf('%sbetas-%s',beta_conditions{i},outputROI);
        
        % replace the data in this map corresponding to the positions
        % of voxels within the ROI we are analysing with the condition-
        % specific support vectors.
        mapcoords(roipos) = meanbetas;
        
        % make any other voxels data 'missing'.
        mapcoords(allcoords)= 0;
        
        % store this updated parameter map as a cell array (needed to
        % match the default format).
        map = {mapcoords};
        mapName = 'OrientationBetas';
        
        % if we are analysing the one participant who has a seperate
        % mrvista session for the MT+ ROI:
        if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
            
            % specify and create a support-vector specific parameter
            % map output directory.
            output_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT/Gray/GLMs/betas/',...
                exp_condition, participant, condition);
            [~,~] = mkdir(output_directory);
            
            % save this new parameter map (within the specific MT+
            % directory), with it's corresponding name, and arbitrary
            % units, and with the coherence thresholds from
            % the original analysis (we never threshold by this, so again,
            % it doesn't really matter).
            save(strcat(output_directory,sprintf('%s%s', participant, mapName)),...
                'map', 'co', 'mapName', 'mapUnits');
            
        else % otherwise, for any other ROI and participant:
            
            % specify and create a support-vector specific parameter
            % map output directory.
            output_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/betas/',...
                exp_condition, participant, condition);
            [~,~] = mkdir(output_directory);
            
            % save this new parameter map, with it's corresponding
            % name, and arbitrary units, and with the coherence
            % thresholds from the original analysis (we never
            % threshold by this, so again, it doesn't really matter).
            save(strcat(output_directory,sprintf('%s%s', participant, mapName)),...
                'map', 'co', 'mapName', 'mapUnits');
        end
    end % repeat for the next feature condition. 
   clear roipos allcoords mapcoords beta_data meanbetas map 
close all % close any open mrvista windows. 
end % repeat for the next participant.

%% --------------------------------------------------------------------- %%