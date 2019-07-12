%% -------------------------- RUN MRVISTA GLM -------------------------- %%

% Specifies number of scans and parfiles, runs GLM and relevant
% contrast-of-parameter estimates.

% KWN 31/07/2018

clear; close;  % general housekeeping.
ynicInit spm8  % add spm to path.
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions'); 

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

Participant = 'R4496'; % specify participant number.

% specify overall project version (can also specify localiser here).
exp_condition = 'original'; 

condition = 'event'; % specify analysis type.

% specify whether to run glm & or contrasts (1-both, 0-contrasts only).
% currently, with the naturalistic 3-feature data, we run contrasts only, 
% and with the naturalistic 3x3 feature data, we do not run contrasts.
runglm = 1; 

%% -------------------------- RUN GLM ---------------------------------- %%

% index into the participant- and condition- specific mrvista directory.
if isequal(exp_condition, 'original') || isequal(exp_condition, 'colour')
    cd(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/', exp_condition, Participant, condition));
elseif isequal(exp_condition, 'localiser')
    cd(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/localiser', exp_condition, Participant));
elseif isequal(condition, '3feature')
    cd(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/', Participant));
elseif isequal(condition, '3x3feature')
    cd(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/multiple', Participant));
end

VOLUME{1} = mrVista('3view'); % open a mrvista gray window.

if runglm == 1 % if we have specified to run a glm:
    
    % get the number of total number of scans the participant completed.
    if isequal(exp_condition, 'original') || isequal(exp_condition, 'colour')
        numScans = viewGet(VOLUME{1}, 'numScans');
        whichScans = 1:numScans; % specify the specific scans to analyse.
        
        % specify the condition-specific parfiles to use.
        parfiles = dir(sprintf('Stimuli/Parfiles/mrvista_%s*.par', condition));
        whichParfs = {parfiles(:).name};
        
    elseif isequal(exp_condition, 'localiser')
        whichScans = 1;
        whichParfs = {'inner_outer_localiser_KWN.par'}; % specify the localiser parfile.
        
    elseif isequal(condition, '3x3feature')
         parfiles = dir(sprintf('Stimuli/Parfiles/*vista*.par'));
         whichParfs = {parfiles(:).name};
    end
    
    VOLUME{1} = selectDataType(VOLUME{1},3); % set the motion compensation data type.
    VOLUME{1} = setCurScan(VOLUME{1},1); % set the first scan view.
    
    % assign each parfile to a specific scan.
    VOLUME{1} = er_assignParfilesToScans(VOLUME{1}, whichScans, whichParfs);
    
    % group all scans to analyse within the motion compensation data type.
    VOLUME{1} = er_groupScans(VOLUME{1}, whichScans, [], 'MotionComp');
    
    % specify the default GLM parameters.
    params = er_defaultParams;
    
    if isequal(exp_condition, 'original') || isequal(exp_condition, 'colour')
        if isequal(condition, 'event') % if this is an event GLM:
            
            % specify the length of each 'block' as a single TR.
            params.eventsPerBlock = 1;
            
        else % otherwise, if this is any other analysis condition:
            
            % specify the length of each block as 5 TRs.
            params.eventsPerBlock = 5;
        end
        
        params.annotation = sprintf('%sGLM',condition);
        
    elseif isequal(exp_condition, 'localiser')
        params.eventsPerBlock = 4;
        params.annotation = 'localiser';
        
    elseif isequal(condition, '3x3feature')
        params.eventsPerBlock = 15;
        params.annotation = sprintf('3x3naturalistic');
    end
    
    % specify some further GLM parameters.
    params.glmHRF = 3; % use spm difference-of-gammas hdr.
    params.framePeriod = 3;
    params.assignParfiles = 0;
    params.lowPassFilter = 0;
    params.setHRFParams = 0;
    saveToDataType = 'GLMs';

    % run the GLM with the specified parameters.
    [VOLUME{1}, newScan] = applyGlm(VOLUME{1}, 'MotionComp', whichScans, params, saveToDataType);
end

% display the GLM data.
VOLUME{1} = selectDataType(VOLUME{1},4);
VOLUME{1} = setCurScan(VOLUME{1},1);
VOLUME{1}=refreshScreen(VOLUME{1});

% specify the condition-specific contrast names and relevant condition
% numbers.
if isequal(condition,'colour')
    contrast.names = {'RGvsIB', 'BYvsIB', 'LUMvsIB', 'RGvsLUM', 'BYvsLUM', 'RGvsBY'};
    contrast.numbers = {[1,5], [2,5], [3,5], [1,3], [2,3], [1,2]};
elseif isequal(condition,'feature') || isequal(condition, 'block')
    contrast.names = {'OvP', 'CvP', 'SvP', 'OvC', 'OvS', 'CvS'};
    contrast.numbers = {[1,4], [2,4], [3,4], [1,2], [1,3], [2,3]};
elseif isequal(condition, 'colour-passive')
    contrast.names = {'RGAvRGP', 'BYAvBYP', 'LUMAvLUMP'};
    contrast.numbers = {[1,2], [3,4], [5,6]};
elseif isequal(condition, '3feature')
    contrast.names = {'FvsO', 'FvsC', 'FvsS', 'FvsP', 'OvsC', 'OvsS', 'OvsP', 'CvsS', 'CvsP', 'SvsP'};
    contrast.numbers = {[1,2], [1,3], [1,4], [1,5], [2,3], [2,4], [2,5], [3,4], [3,5], [4,5]};
end

% if this is not an event- analysis run (as we are currently not specifying
% any contrasts for event analyses).
if ~isequal(exp_condition, 'localiser')
    if ~isequal(condition,'event')
        % for each of these specified contrasts, compute the contrast-of-parameter
        % estimate.
        for cope = 1:length(contrast.names)
            computeContrastMap2(VOLUME{1}, contrast.numbers{cope}(1), contrast.numbers{cope}(2), contrast.names{cope});
        end
    end
elseif isequal(exp_condition, 'localiser')
    computeContrastMap2(VOLUME{1}, 1, 0, 'InnervsFix');
    computeContrastMap2(VOLUME{1}, 2, 0, 'OutervsFix');
    computeContrastMap2(VOLUME{1}, 1, 2, 'InnervsOuter');
    computeContrastMap2(VOLUME{1}, 2, 1, 'OutervsInner');
end

close all; % close all mrvista windows.

%% --------------------------------------------------------------------- %%