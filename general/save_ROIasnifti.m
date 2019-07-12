%% ----------------- SAVE MR-VISTA ROIs AS NIFTI FILES ----------------- %%

% loads in gray-view retinotopically-defined ROIs in mr-vista and saves
% each ROI to a .nii.gz file for use in other analysis programmes.

% 10/7/2019 KWN

clear; clc; close all

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

% specify participants to analyse.
participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455',...
    'R3517', 'R3773', 'R3932', 'R4065', 'R4496'};

% specify ROIs to process.
ROIs{1} = '01_Combined_V1_KWN.mat';
ROIs{2} = '02_Combined_V2_KWN.mat';
ROIs{3} = '03_Combined_V3_KWN.mat';
ROIs{4} = '06_Combined_V3A_KWN.mat';
ROIs{5} = '07_Combined_V3B_KWN.mat';
ROIs{6} = '08_Combined_V4_KWN.mat';
ROIs{7} = '09_Combined_LO1_KWN.mat';
ROIs{8} = '10_Combined_LO2_KWN.mat';
ROIs{9} = '11_Combined_A1_KWN.mat';
ROIs{10} = '12_Combined_MT+_KWN.mat';

%% ------------------------- SAVE ROI AS NIFTI ------------------------- %%

for ppnum = 1:length(participants) % for each participant in turn:
    pp = participants{ppnum};
    
    % index into participant-specific mrvista directory and load session.
    cd(sprintf('/scratch/sg3/P1323/original/fMRI/%s/mrvista/block', pp));
    VOLUME{1} = mrVista('3view');
    
    % specify an output directory for nifti ROIs.
    outputdir = sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/anatomy/niftiROIs/original/', pp);
    [~,~] = mkdir(outputdir);
    
    for ROInum = 1:length(ROIs) % for each ROI in turn:
        ROI = ROIs{ROInum};
        
        % specify the participant-specific ROI directory.
        if isequal(ROI, '12_Combined_MT+_KWN.mat') && isequal(pp, 'R3111')
            close all % close the current mrvista window.
            
            % index into and open the condition-specific MT directory.
            cd('/scratch/groups/Projects/P1323/original/fMRI/R3111/mrvista/block/MT');
            VOLUME{1} = mrVista('3view');
            
            % specify the MT-specific ROI storage directory.
            ROI_dir_pp = {strcat('/scratch/groups/Projects/P1323/original/fMRI/R3111/anatomy/MT/ROIs/', ROI)};
        else
            ROI_dir_pp = {strcat(sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/anatomy/ROIs/', pp), ROI)};
        end
        
        % load the ROI.
        VOLUME{1}=loadROI(VOLUME{1},ROI_dir_pp,[],[],1,0); VOLUME{1}=refreshScreen(VOLUME{1},0);
        
        % save the ROI to a nifti file. 
        roiSaveForItkGray(VOLUME{1},strcat(outputdir, ROI, '.nii.gz'));
    end
    close all % close the mr-vista session.
end

%% --------------------------------------------------------------------- %%