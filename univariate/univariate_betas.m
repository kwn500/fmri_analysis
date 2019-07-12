%% ------------------- fMRI UNIVARIATE BETA ANALYSIS ------------------- %%

% For specified participants and ROIs, we extract the univariate
% condition-specific beta values and percentage variance explained from
% the GLM output and look for significant differences between conditions
% using one-way repeated measures ANOVAs.

% KWN 27/7/18

clear; close; % general housekeeping.
ynicInit spm8; % add spm to path.
addpath('/scratch/sg3/P1323/code/fmri_analysis/functions');

%% -------------------- SPECIFY PARTICIPANTS & ROIS -------------------- %%

exp_condition = 'colour'; % specify name of overall experiment (original/colour).
analysis_condition = 'colour'; % specify analysis condition of interest.

% specify the interaction analysis to run (1:4) for the naturalistic 3x3
% feature analysis.
interaction_toggle = 1;

% specify list of participant R numbers.
if isequal(exp_condition, 'original')
    participants{1} = 'R2268';   valid_runs{1} = 1:6;
    participants{2} = 'R2548';   valid_runs{2} = 1:6;
    participants{3} = 'R2590';   valid_runs{3} = 1:6;
    participants{4} = 'R2904';   valid_runs{4} = 1:6;
    participants{5} = 'R3111';   valid_runs{5} = 1:5;
    participants{6} = 'R3455';   valid_runs{6} = 1:6;
    participants{7} = 'R3517';   valid_runs{7} = 1:7;
    participants{8} = 'R3773';   valid_runs{8} = 1:6;
    participants{9} = 'R3932';   valid_runs{9} = 1:6;
    participants{10} = 'R4065';  valid_runs{10} = 1:6;
    %     participants{11} = 'R4127';  valid_runs{11} = 1:6;
    participants{11} = 'R4496';  valid_runs{11} = 1:6;
    
elseif isequal(exp_condition, 'colour')
    participants{1} = 'R2268';  valid_runs{1} = 1:8;
    participants{2} = 'R2548';  valid_runs{2} = 1:8;
    participants{3} = 'R2590';  valid_runs{3} = 1:7;
    participants{4} = 'R2904';  valid_runs{4} = 1:8;
    participants{5} = 'R3111';  valid_runs{5} = 1:7;
    participants{6} = 'R3517';  valid_runs{6} = 1:8;
    participants{7} = 'R3773';  valid_runs{7} = 1:8;
    participants{8} = 'R4059';  valid_runs{8} = 1:8;
    participants{9} = 'R4065';  valid_runs{9} = 1:8;
    participants{10} = 'R4244'; valid_runs{10} = 1:8;
    participants{11} = 'R4831'; valid_runs{11} = 1:8;
    participants{12} = 'R4928'; valid_runs{12} = 1:8;
    
elseif isequal(analysis_condition, '3x3feature')
    participants{1} = 'R2548';   valid_runs{1} = 1:8;
    participants{2} = 'R2904';   valid_runs{2} = 1:8;
    participants{3} = 'R3111';   valid_runs{3} = 1:6;
    participants{4} = 'R3517';   valid_runs{4} = 1:8;
    participants{5} = 'R3773';   valid_runs{5} = 1:7;
    participants{6} = 'R4059';   valid_runs{6} = 1:7;
    participants{7} = 'R4065';   valid_runs{7} = 1:8;
    participants{8} = 'R4127';   valid_runs{8} = 1:8;
    participants{9} = 'R4244';   valid_runs{9} = 1:7;
    participants{10} = 'R4829';  valid_runs{10} = 1:8;
    participants{11} = 'R4831';  valid_runs{11} = 1:7;
    participants{12} = 'R4833';  valid_runs{12} = 1:8;
    participants{13} = 'R4890'; valid_runs{13} = 1:8;
    participants{14} = 'R4928'; valid_runs{14} = 1:6;
    participants{15} = 'R5006'; valid_runs{15} = 1:8;
end

% we specify the anatomy locations for our naturalistic analysis.
if isequal(exp_condition, 'naturalistic')
    anat{1} = '/scratch/groups/Projects/P1361/fMRI/R2548/anatomy/T1.nii.gz';
    anat{2} = '/scratch/groups/Projects/P1361/fMRI/R2904/anatomy/nu.nii.gz';
    anat{3} = '/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/MT/nu_fast_restore.nii.gz';
    anat{4} = '/scratch/groups/Projects/P1361/fMRI/R3517/anatomy/nu.nii.gz';
    anat{5} = '/scratch/groups/Projects/P1361/fMRI/R3773/anatomy/T1.nii.gz';
    anat{6} = '/scratch/groups/Projects/P1361/fMRI/R4059/anatomy/T1.nii.gz';
    anat{7} = '/scratch/groups/Projects/P1361/fMRI/R4065/anatomy/t1.nii.gz';
    anat{8} = '/scratch/groups/Projects/P1361/fMRI/R4127/anatomy/T1.nii.gz';
    anat{9} = '/scratch/groups/Projects/P1361/fMRI/R4244/anatomy/t1.nii.gz';
    anat{10} = '/scratch/groups/Projects/P1361/fMRI/R4829/anatomy/t1.nii.gz';
    anat{11} = '/scratch/groups/Projects/P1361/fMRI/R4831/anatomy/t1.nii.gz';
    anat{12} = '/scratch/groups/Projects/P1361/fMRI/R4833/anatomy/t1.nii.gz';
    anat{13} = '/scratch/groups/Projects/P1361/fMRI/R4890/anatomy/t1.nii.gz';
    anat{14} = '/scratch/groups/Projects/P1361/fMRI/R4928/anatomy/t1.nii.gz';
    anat{15} = '/scratch/groups/Projects/P1361/fMRI/R5006/anatomy/t1.nii.gz';
end
% specify list of ROIs we wish to extract univariate data from, as well as
% an abbreviated version of the ROI name for naming output files.
ROI_names{1} = '01_Combined_V1_KWN.mat'; ROI_output{1} = 'V1';
ROI_names{2} = '02_Combined_V2_KWN.mat'; ROI_output{2} = 'V2';
ROI_names{3} = '03_Combined_V3_KWN.mat'; ROI_output{3} = 'V3';
ROI_names{4} = '06_Combined_V3A_KWN.mat'; ROI_output{4} = 'V3A';
ROI_names{5} = '07_Combined_V3B_KWN.mat'; ROI_output{5} = 'V3B';
ROI_names{6} = '08_Combined_V4_KWN.mat'; ROI_output{6} = 'V4';
ROI_names{7} = '09_Combined_LO1_KWN.mat'; ROI_output{7} = 'LO1';
ROI_names{8} = '10_Combined_LO2_KWN.mat'; ROI_output{8} = 'LO2';
ROI_names{9} = '11_Combined_A1_KWN.mat'; ROI_output{9} = 'A1';
ROI_names{10} = '12_Combined_MT+_KWN.mat'; ROI_output{10} = 'MT+';
ROI_names{11} = '14_Combined_IPS0_KWN.mat'; ROI_output{11} = 'IPS0';

% specify the output directory for extracted betas and variance explained.
if isequal(exp_condition, 'original') || isequal(exp_condition, 'colour')
    output_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/',...
        exp_condition, analysis_condition);
elseif isequal(exp_condition, 'naturalistic')
    output_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/');
end

% for the interaction analysis in the colour experiment, we
% specify the specific conditions we want to compare (it isn't particularly
% useful to compare all 13 conditions to each other- we take a more
% targeted approach), and their corresponding names, otherwise, we specify
% some default values.
if isequal(analysis_condition, 'colourxfeature')
    condition_numbers = {[1 2 3], [5 6 7], [9 10 11], [1 5 9], [2 6 10], [3 7 11]};
    condition_names = {'red-green_feature', 'blue-yellow_feature', 'luminance_feature',...
        'orientation_colour', 'contrast_colour', 'shape_colour'};
elseif isequal(analysis_condition, '3x3feature')
    
    % specify the analysis-specific interaction condition numbers and name.
    if interaction_toggle == 1
        condition_numbers = {[2 3 4]};
        condition_names = {'orientation3x3'};
    elseif interaction_toggle == 2
        condition_numbers = {[5 6 7]};
        condition_names = {'colour3x3'};
    elseif interaction_toggle == 3
        condition_numbers = {[8 9 10]};
        condition_names = {'shape3x3'};
    elseif interaction_toggle == 4
        condition_numbers = {[1 11]};
        condition_names = {'passive_face'};
    end
else
    condition_numbers = {1};
    condition_names = {'NA'};
end

% for each particular comparison in turn:
for condition_iteration = 1:length(condition_names)
    
    % extract the relevant comparison name and numbers corresponding to the
    % columns of data we wish to analyse (the condition-specific betas).
    condition_number = condition_numbers{condition_iteration};
    condition_name = condition_names{condition_iteration};
    
    if isequal(analysis_condition, 'colourxfeature')
        % specify and create a comparison-specific output directory.
        output_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/',...
            exp_condition, analysis_condition);
        output_dir = strcat(output_dir, sprintf('%s/', condition_name));
    elseif isequal(analysis_condition, '3x3feature')
        output_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3x3feature/');
        output_dir = strcat(output_dir, sprintf('%s/', condition_name));
    end
    
    for ROInum = 1:length(ROI_names) % for each ROI in turn:
        
        % specify the ROI-specific name and abbreviated output version.
        ROI = ROI_names{ROInum};
        outputROI = ROI_output{ROInum};
        
        for ppnum = 1:length(participants) % for each participant in turn:
            
            % specify the participant-specific anatomy if naturalistic
            % analysis.
            if isequal(exp_condition, 'naturalistic')
                anatomy = anat{ppnum};
            end
            
            % record the participant-specific R number.
            participant = participants{ppnum};
            
            ROIstatus = 1; % set the ROI status variable to 'present'.
            
            % for one particular participant, their MT+ ROI is associated with a
            % different high resolution structural. if this is the participant
            % we are currently analysing, and we are looking at MT+ data, set
            % the mrvista directory and ROI locations accordingly.
            if isequal(participant, 'R3111') && isequal(ROI, '12_Combined_MT+_KWN.mat')
                
                if ~(isequal(exp_condition, 'naturalistic'))
                    mrvista_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT/',...
                        exp_condition, participant, analysis_condition);
                    ROI_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/MT/ROIs/',...
                        exp_condition, participant);
                elseif isequal(exp_condition, 'naturalistic')
                    mrvista_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/R3111/mrvista/multiple/MT/');
                    ROI_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/MT/ROIs/');
                end
            elseif isequal(participant, 'R3517') && isequal(exp_condition, 'naturalistic')
                mrvista_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans/multiple/';
                ROI_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/anatomy/ROIs/';
                
            else % otherwise, in all other cases:
                
                % specify the directory containing the mrvista session of
                % interest.
                mrvista_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/',...
                    exp_condition, participant, analysis_condition);
                ROI_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs/',...
                    exp_condition, participant);
                
            end
            cd(mrvista_dir); % index into the mrvista session directory.
            
            % for one participant, mrvista does not by default select their
            % corresponding high-resolution structural scan, for this
            % participant, we specify the location of this file whilst opening
            % a mrvista gray window.
            if ~isequal(exp_condition, 'naturalistic')
                if isequal(participant,'R2548')
                    
                    if isequal(exp_condition, 'original')
                        anat = '/scratch/groups/Projects/P1323/original/fMRI/R2548/anatomy/T1.nii.gz';
                    elseif isequal(exp_condition, 'colour')
                        anat = '/scratch/groups/Projects/P1323/colour/fMRI/R2548/anatomy/nu.nii.gz';
                    end
                    VOLUME{1} = mrVista('3view', anat);
                    
                else % otherwise, we just open a mrvista gray 3-view window.
                    VOLUME{1} = mrVista('3view');
                end
            else
                VOLUME{1} = mrVista('3view');
            end
            
            % specify the path to the specific ROI of interest.
            ROIs = {strcat(ROI_dir, '/', ROI)};
            
            % load this ROI into the mrvista gray window.
            VOLUME{1}=loadROI(VOLUME{1},ROIs,[],[],1,0); VOLUME{1}=refreshScreen(VOLUME{1},0);
            
            % set the view to the first scan of the GLM data type.
            VOLUME{1} = selectDataType(VOLUME{1},4);
            VOLUME{1} = setCurScan(VOLUME{1},1);
            
            % try to visualise the GLM results and dump the associated data to
            % the workspace. if the ROI contains no voxels (which may happen
            % with inner/outer ROI splits in higher visual areas), catch the
            % error and toggle the ROI status variable to reflect this.
            try
                tc_plotScans(VOLUME{1},1);
                tc_visualizeGlm;
                tc_dumpDataToWorkspace_univariate_KWN;
            catch
                ROIstatus = 0;
            end
            
            if ROIstatus == 1 % if the ROI of interest did contain voxels:
                
                % extract the condition-specific univariate betas (for
                % orientation, contrast, shape and passive respectively), and
                % the percent variance explained in terms of the GLM within
                % this ROI (one value- not condition specific).
                if isequal(analysis_condition, 'colourxfeature') || isequal(analysis_condition, '3x3feature')
                    betas = uvoxel_data.glm.betas(condition_number);
                elseif isequal(analysis_condition, '3feature')
                    betas = uvoxel_data.glm.betas(1:5);
                else
                    betas = uvoxel_data.glm.betas(1:4);
                end
                varexp = uvoxel_data.glm.varianceExplained;
                
                % otherwise, if the ROI of interest did not contain voxels:
            elseif ROIstatus == 0
                
                % set the extracted beta values and percentage variance
                % explained to default missing values.
                betas = [NaN NaN NaN NaN];
                varexp = NaN;
            end
            
            % specify a participant-specific output directory and create this
            % folder.
            pp_output_dir = strcat(output_dir, participant, '/'); [~,~] = mkdir(pp_output_dir);
            
            % save the participant and ROI-specific extracted univariate betas
            % and percentage variance explained.
            save(strcat(pp_output_dir, outputROI, '_', analysis_condition, '_univariate_betas.mat'), 'betas');
            save(strcat(pp_output_dir, outputROI, '_', analysis_condition, '_univariate_varexp.mat'),'varexp');
            
            % store the individual participant and ROI-specific information in
            % a participant-specific location in overall storage variables.
            betas_pp(ppnum,:) = betas;
            varexp_pp(ppnum,:) = varexp;
            
            % refresh these extracted variables for the next participant.
            clear betas varexp
            
            close all % close all figures and mrvista windows.
        end % continue to the next participant.
        
        % specify and create a group output directory.
        group_output_dir = strcat(output_dir, 'group/'); [~,~] = mkdir(group_output_dir);
        
        % save the ROI-specific univariate betas and percentage variance
        % explained across all participants to .mat files.
        save(strcat(group_output_dir, outputROI, '_', analysis_condition, '_univariate_betas.mat'),...
            'participants', 'betas_pp');
        save(strcat(group_output_dir, outputROI, '_', analysis_condition, '_univariate_varexp.mat'),...
            'participants', 'varexp_pp');
        
        % specify and create the rmanova and figure-specific output directories.
        stats_output_dir = strcat(group_output_dir, 'stats/'); [~,~] = mkdir(stats_output_dir);
        figure_output_dir = strcat(stats_output_dir, 'figures/'); [~,~] = mkdir(figure_output_dir);
        
        % run one-way repeated measures anovas to assess the differences
        % between orientation, contrast and shape attentional modulation, and
        % plot the data (with corresponding significance) as bar charts. the
        % original and normalised toggles reflect whether we subtract the
        % passive condition activation from the orientation, contrast and shape
        % columns- this doesn't make any difference to the significance of
        % results.
        if ~isequal(exp_condition, 'naturalistic')
            if isequal(analysis_condition, 'colourxfeature')
                rmanova_plot_interaction(betas_pp, outputROI, analysis_condition, 'original',...
                    stats_output_dir, figure_output_dir, condition_name);
            else
                rmanova_plot(betas_pp, outputROI, analysis_condition, 'original', stats_output_dir, figure_output_dir, condition_name);
                rmanova_plot(betas_pp, outputROI, analysis_condition, 'normalised', stats_output_dir, figure_output_dir, condition_name);
            end
        end
        
        % otherwise, we analyse the naturalistic 3feature data in the same
        % way, but call it here for ease.
        if isequal(exp_condition, 'naturalistic')
            
            if isequal(analysis_condition, '3feature')
                extracted_betas = dir('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/group/*univariate_betas.mat');
            elseif isequal(analysis_condition, '3x3feature')
                extracted_betas = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3x3feature/%s/group/*univariate_betas.mat',...
                    condition_name));
            end
            betafiles = {extracted_betas.name};
            
            for i = 1:length(betafiles)
                
                % for each ROI, we load in the univariate beta data and
                % perform a one-way repeated measures anova assessing the
                % difference in BOLD signal modulation across attention
                % conditions.
                
                if isequal(analysis_condition, '3feature')
                    load(strcat('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/group/', betafiles{i}));
                elseif isequal(analysis_condition, '3x3feature')
                    load(strcat(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3x3feature/%s/group/',...
                        condition_name), betafiles{i}));
                end
                ROI = betafiles{i}(1:end-21);
                betadata = betas_pp;
                
                if isequal(analysis_condition, '3feature')
                    [p{i},f{i}, dfm{i}, dfr{i}, means{i}, stderr{i}, mauchlyp{i}] = rmanova_plot_naturalistic(betadata,...
                        betadata, 'original', stats_output_dir, figure_output_dir, ROI);
                    
                    [~, ~, ~, adj_p]=fdr_bh(p);
                    
                    % store and save the important analysis information to
                    % a .mat file. 
                    univariate_block.p = p;
                    univariate_block.mauchlyp = mauchlyp;
                    univariate_block.adj_p = adj_p;
                    univariate_block.f = f;
                    univariate_block.dfm = dfm;
                    univariate_block.dfr = dfr;
                    univariate_block.means = means;
                    univariate_block.stderrs = stderr;
                    
                    save('/scratch/groups/Projects/P1361/fMRI/vista_output/output/univariate_naturalistic3feature-mauchly.mat',...
                        'univariate_block');
                    
                elseif isequal(analysis_condition, '3x3feature')
                    [p{i},f{i}, dfm{i}, dfr{i}, means{i}, stderr{i}, mauchlyp{i}] = rmanova_plot_naturalistic3x3adj_p(betadata,...
                        stats_output_dir, figure_output_dir, ROI, toggle, adj_p);
                    adj_p = adj_p_all(i);
                    
                    % store and save the important analysis information to
                    % a .mat file. 
                    univariate_block.p = p;
                    univariate_block.reorder_p = reorder_p;
                    univariate_block.mauchlyp = mauchlyp;
                    univariate_block.adj_p = adj_p_all;
                    univariate_block.f = f;
                    univariate_block.dfm = dfm;
                    univariate_block.dfr = dfr;
                    univariate_block.means = means;
                    univariate_block.stderrs = stderr;
                    
                    save(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/output/univariate_naturalistic3x3feature-mauchly_%s.mat',condition_name),...
                        'univariate_block');
                end
            end
        end
        
        clear betas_pp varexp_pp % refresh these group variables for the next iteration.
        close all % close all figures and mrvista windows.
    end % continue to the next ROI.
end % continue to the next comparison.

%% --------------------------------------------------------------------- %%