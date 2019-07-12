%% ------------- BISECT ROIS TO 'INNER' & 'OUTER' REGIONS -------------- %% 

% loads in inner-vs-outer localiser parameter maps and takes all voxels
% within a specific ROI corresponding to positive values within the
% parameter map as 'inner' ROIs, and negative values as 'outer' ROIs, and
% saves these bisected ROIs in the individual participants ROI directory. 

% KWN 31/07/2018

clear; close; % general housekeeping. 
ynicInit spm8 % add spm to path.

%% ----------------------- PARAMETERS TO EDIT -------------------------- %%

% specify participants to analyse. 
pp_list = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4496'}; 

exp_condition = 'colour'; % specify overall project version. 

% specify ROIs to restrict. 
ROI_names{1} = '01_Combined_V1_KWN.mat';
ROI_names{2} = '02_Combined_V2_KWN.mat';
ROI_names{3} = '03_Combined_V3_KWN.mat';
ROI_names{4} = '06_Combined_V3A_KWN.mat';
ROI_names{5} = '07_Combined_V3B_KWN.mat';
ROI_names{6} = '08_Combined_V4_KWN.mat';
ROI_names{7} = '09_Combined_LO1_KWN.mat';
ROI_names{8} = '10_Combined_LO2_KWN.mat';

%% ------------------------ RESTRICT INNER ROIS ------------------------ %%

for i = 1:length(pp_list) % for each participant in turn:
    
    Participant = pp_list{i}; % specify the participant number. 
    
    % index into the participant's individiual mrvista localiser directory.
    cd(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/localiser', exp_condition, Participant));
    
    VOLUME{1} = mrVista('3view'); % open a mrVista gray 3-view window.
    
    for x = 1:length(ROI_names) % for each ROI in turn:

        % Here, to keep only the positive (inner) activation of the
        % contrast, we set the upper threshold at 2 -log(p), so that any
        % inner activation above this threshold remains. We set the lower
        % threshold to the maximum negative value, so that no outer
        % activation remains at this threshold. 
        
        % load the participant-specific inner-vs-outer parameter map. 
        VOLUME{1} = loadParameterMap(VOLUME{1}, ...
            sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/localiser/Gray/GLMs/InnervsOuter.mat',...
            exp_condition, Participant));
        
        VOLUME{1}= refreshScreen(VOLUME{1}); % refresh the window.
        
        % set a bicolour hot- cold parameter map colour bar, and manually 
        % clip these values between [-20 20];
        VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'coolhotCmap'); 
        VOLUME{1}=refreshScreen(VOLUME{1}, 1);
        VOLUME{1} = setClipMode(VOLUME{1}, 'map', [-20 20]);
        
        % get the maximum and minimumm values for this parameter map. 
        sliderMin = get(VOLUME{1}.ui.mapWinMin.sliderHandle, 'Value');
        sliderMax = get(VOLUME{1}.ui.mapWinMax.sliderHandle, 'Value');

        % set the minimum (the positive threshold) to 2 -log(p) value. 
        set(VOLUME{1}.ui.mapWinMin.sliderHandle, 'Min', 2);
        
        % set the maximum (the negative threshold) to the default negative value. 
        set(VOLUME{1}.ui.mapWinMax.sliderHandle, 'Max', sliderMin);
        
        % set these values in the mrVista window. 
        setSlider(VOLUME{1}, VOLUME{1}.ui.mapWinMin);
        setSlider(VOLUME{1}, VOLUME{1}.ui.mapWinMax);
        
        % specify the individual participant's ROI directory. 
        ROI_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs/', exp_condition, Participant);
        
        % specify the path to an individual ROI. 
        ROIs = {strcat(ROI_directory, ROI_names{x})};
        
        % load this ROI into the mrVista gray view window.
        VOLUME{1}=loadROI(VOLUME{1},ROIs,[],[],1,0); VOLUME{1}=refreshScreen(VOLUME{1},0);
        
        % restrict this ROI to the underlying inner activation. 
        VOLUME{1}= restrictROIfromMenu(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
        
        % edit this restricted ROI name and colour. 
        VOLUME{1}.ROIs.name = strcat(ROI_names{x}(1:end-4),'_Inner'); VOLUME{1}.ROIs.color = 'g';
        
        % save this restricted ROI in the individual participants ROI directory. 
        saveROI(VOLUME{1},VOLUME{1}.ROIs(VOLUME{1}.selectedROI),0);
        
        % remove this restricted ROI from the gray 3 view window. 
        VOLUME{1}=deleteROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
        
%% ------------------------ RESTRICT OUTER ROIS ------------------------ %%

        % repeat this process for the outer ROIs. Here, we set the maximum
        % value to the default positive value (to remove any inner
        % activation), and the minimum value to -2, to keep any outer
        % activation above this threshold. The same process as above is
        % then repeated. 
        
        VOLUME{1} = loadParameterMap(VOLUME{1},...
            sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/localiser/Gray/GLMs/InnervsOuter.mat',...
            exp_condition, Participant));
        VOLUME{1}= refreshScreen(VOLUME{1});

        VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'coolhotCmap'); 
        VOLUME{1}=refreshScreen(VOLUME{1}, 1);
        VOLUME{1} = setClipMode(VOLUME{1}, 'map', [-20 20]);
        
        set(VOLUME{1}.ui.mapWinMin.sliderHandle, 'Min', sliderMax);
        set(VOLUME{1}.ui.mapWinMax.sliderHandle, 'Max', -2);
        setSlider(VOLUME{1}, VOLUME{1}.ui.mapWinMin);
        setSlider(VOLUME{1}, VOLUME{1}.ui.mapWinMax);
        
        ROI_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs/', exp_condition, Participant);
        
        ROIs = {strcat(ROI_directory,ROI_names{x})};
        VOLUME{1}=loadROI(VOLUME{1},ROIs,[],[],1,0); VOLUME{1}=refreshScreen(VOLUME{1},0);
        
        VOLUME{1}= restrictROIfromMenu(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
        VOLUME{1}.ROIs.name = strcat(ROI_names{x}(1:end-4),'_Outer');
        VOLUME{1}.ROIs.color = 'r';
        
        saveROI(VOLUME{1},VOLUME{1}.ROIs(VOLUME{1}.selectedROI),0);
        
        VOLUME{1}=deleteROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
    end % repeat for next ROI.
    close all
end % repeat for next participant. 

%% --------------------------------------------------------------------- %%