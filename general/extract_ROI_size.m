%% ----------------- EXTRACT NUMBER OF VOXELS IN ROIS ------------------ %%

% runs for all current participants (across experiments) and records the
% size (in voxels) of all ROIs (V1-MT+), for original (dorsal/ventral),
% inner/outer and combined ROIs. 

% KWN 27/07/18

clear; clc; % general housekeeping.
ynicInit spm8 % add spm to path.
addpath('/scratch/groups/Projects/P1323/code'); % add function directory to path.

%% -------------------- SPECIFY ANALYSIS PARAMETERS  ------------------- %%

base = '/scratch/groups/Projects/'; % specify overall project directory.

% specify participant r-numbers.
pp_name{1} = 'R2268';  
pp_name{2} = 'R2548';  
pp_name{3} = 'R2590';  
pp_name{4} = 'R2904';  
pp_name{5} = 'R3111'; 
pp_name{6} = 'R3455'; 
pp_name{7} = 'R3517'; 
pp_name{8} = 'R3773';  
pp_name{9} = 'R3932'; 
pp_name{10} = 'R4059';  
pp_name{11} = 'R4065';
pp_name{12} = 'R4127'; 
pp_name{13} = 'R4244';
pp_name{14} = 'R4496';
pp_name{15} = 'R4829'; 
pp_name{16} = 'R4831';
pp_name{17} = 'R4833'; 
pp_name{18} = 'R4890';
pp_name{19} = 'R4928'; 
pp_name{20} = 'R5006';

% specify each participants roi-storage directory.
ROI_dir{1} = [base,'P1323/original/fMRI/R2268/anatomy/ROIs']; 
ROI_dir{2} = [base,'P1323/original/fMRI/R2548/anatomy/ROIs']; 
ROI_dir{3} = [base, 'P1323/original/fMRI/R2590/anatomy/ROIs']; 
ROI_dir{4} = [base,'P1323/original/fMRI/R2904/anatomy/ROIs'];
ROI_dir{5} = [base,'P1323/original/fMRI/R3111/anatomy/ROIs'];  
ROI_dir{6} = [base,'P1323/original/fMRI/R3455/anatomy/ROIs']; 
ROI_dir{7} = [base,'P1323/original/fMRI/R3517/anatomy/ROIs'];  
ROI_dir{8} = [base,'P1323/original/fMRI/R3773/anatomy/ROIs']; 
ROI_dir{9} = [base,'P1323/original/fMRI/R3932/anatomy/ROIs'];  
ROI_dir{10} = [base,'P1323/colour/fMRI/R4059/anatomy/ROIs'];  
ROI_dir{11} = [base,'P1323/original/fMRI/R4065/anatomy/ROIs']; 
ROI_dir{12} = [base,'P1323/original/fMRI/R4127/anatomy/ROIs']; 
ROI_dir{13} = [base,'P1323/colour/fMRI/R4244/anatomy/ROIs'];  
ROI_dir{14} = [base,'P1323/original/fMRI/R4496/anatomy/ROIs']; 
ROI_dir{15} = [base,'P1361/fMRI/R4829/anatomy/mrvista/ROIs'];  
ROI_dir{16} = [base,'P1323/colour/fMRI/R4831/anatomy/ROIs']; 
ROI_dir{17} = [base,'P1361/fMRI/R4833/anatomy/mrvista/ROIs'];
ROI_dir{18} = [base,'P1361/fMRI/R4890/anatomy/mrvista/ROIs'];  
ROI_dir{19} = [base,'P1323/colour/fMRI/R4928/anatomy/ROIs'];  
ROI_dir{20} = [base,'P1361/fMRI/R5006/anatomy/mrvista/ROIs']; 

% specify each participants mr-vista session directory. 
vista_dir{1} = [base, 'P1323/original/fMRI/R2268/mrvista/block'];
vista_dir{2} = [base, 'P1323/original/fMRI/R2548/mrvista/block'];
vista_dir{3} = [base, 'P1323/original/fMRI/R2590/mrvista/block'];
vista_dir{4} = [base, 'P1323/original/fMRI/R2904/mrvista/block'];
vista_dir{5} = [base, 'P1323/original/fMRI/R3111/mrvista/block'];
vista_dir{6} = [base, 'P1323/original/fMRI/R3455/mrvista/block'];
vista_dir{7} = [base, 'P1323/original/fMRI/R3517/mrvista/block'];
vista_dir{8} = [base, 'P1323/original/fMRI/R3773/mrvista/block'];
vista_dir{9} = [base, 'P1323/original/fMRI/R3932/mrvista/block'];
vista_dir{10} = [base, 'P1323/colour/fMRI/R4059/mrvista/colour'];
vista_dir{11} = [base, 'P1323/original/fMRI/R4065/mrvista/block'];
vista_dir{12} = [base, 'P1323/original/fMRI/R4127/mrvista/block'];
vista_dir{13} = [base, 'P1323/colour/fMRI/R4244/mrvista/colour'];
vista_dir{14} = [base, 'P1323/original/fMRI/R4496/mrvista/block'];
vista_dir{15} = [base, 'P1361/fMRI/R4829/retinotopy/mrvista'];
vista_dir{16} = [base, 'P1323/colour/fMRI/R4831/mrvista/colour'];
vista_dir{17} = [base, 'P1361/fMRI/R4833/retinotopy/mrvista'];
vista_dir{18} = [base, 'P1361/fMRI/R4890/retinotopy/mrvista'];
vista_dir{19} = [base, 'P1323/colour/fMRI/R4928/mrvista/colour'];
vista_dir{20} = [base, 'P1361/fMRI/R5006/retinotopy/mrvista'];

% specify name of ROIs to analyse. 
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

%% ------------------------- RECORD ROI SIZE --------------------------- %%

for i = 1:length(pp_name) % for each participant in turn:
    Participant = pp_name{i}; % specify the participant number.
    
    % index into the participant-specific mrvista session directory.
    cd(vista_dir{i});
    
    % open a gray window, and set the view to the first scan of the GLM
    % data type.
    VOLUME{1} = mrVista('3view');
    VOLUME{1} = selectDataType(VOLUME{1},4);
    VOLUME{1} = setCurScan(VOLUME{1},1);
    VOLUME{1}= refreshScreen(VOLUME{1});
    
    for x = 1:length(ROIs) % for each ROI in turn:
        
        noROI = 1; % set a default value for the ROI present/absent toggle.
        
        % for one particular participant, their MT+ ROI is associated with a
        % different high resolution structural. if this is the participant
        % we are currently analysing, and we are looking at MT+ data, set
        % the mrvista directory and ROI locations accordingly.
        if isequal(ROIs{x}, '12_Combined_MT+_KWN.mat') && isequal(Participant, 'R3111')
            close all % close the current mrvista window. 
            
            % index into and open the condition-specific MT directory. 
            cd('/scratch/groups/Projects/P1323/original/fMRI/R3111/mrvista/block/MT');
            VOLUME{1} = mrVista('3view');
            VOLUME{1} = selectDataType(VOLUME{1},4);
            VOLUME{1} = setCurScan(VOLUME{1},1);
            VOLUME{1}= refreshScreen(VOLUME{1});
            
            % specify the MT-specific ROI storage directory. 
            ROI_dir_pp = {strcat('/scratch/groups/Projects/P1323/original/fMRI/R3111/anatomy/MT/ROIs/',...
                ROIs{x})};
            
        else % otherwise, for all other participants and ROIs. 
            
            % specify a path to the specific ROI of interest.
            ROI_dir_pp = {strcat(ROI_dir{i}, '/', ROIs{x})};
        end
        
        % load this ROI into the gray window.
        VOLUME{1}=loadROI(VOLUME{1},ROI_dir_pp,[],[],1,0); VOLUME{1}=refreshScreen(VOLUME{1},0);
        
        % try to visualise the univariate analysis results for this ROI,
        % if this fails(because there are no voxels within this ROI
        % (common with the higher-order visual areas in the inner/outer
        % splits), catch the error and change the no ROI toggle to
        % indicate the absence of this particular ROI.
        try
            tc_plotScans(VOLUME{1},1);
        catch
            noROI = 0;
        end
        
        % if this ROI did not exist, set the ROI size variable to reflect
        % this.
        if noROI == 0
            ind_ROI_size = noROI;
            
        else % otherwise, if the ROI did exist:
            
            % within this plot window, visualize the GLM results.
            tc_visualizeGlm;
            
            % dump this data to the workspace.
            tc_dumpDataToWorkspace_univariate_KWN;
            
            % extract the ROI size from this data.
            ind_ROI_size = size(uvoxel_data.roi.coords,2);
        end
        
        % store the individual ROI size (in voxels) in an ROI- and
        % participant- specific location in an overall storage matrix.
        ROI_sizes(i,x) = ind_ROI_size;
        
    end % continue to the next ROI.
    close all % close all figures and the mrvista window.
    clear ROI_dir_pp % refresh some variables for the next loop iteration. 
    
end % continue to the next participant.

% save the across-individual ROI storage variable to a .mat file along with
% the participants analysed. 
save('/scratch/groups/Projects/P1323/original/fMRI/general_output/ROI_sizes_innerouter.mat', 'ROI_sizes');

%% --------------------------------------------------------------------- %%