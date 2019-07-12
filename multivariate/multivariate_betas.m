%% -------- EXTRACT MULTIVARIATE (PER VOXEL) TIMESERIES DATA ----------- %%

% For each participant and each ROI specified, we extract the timeseries
% data for each voxel across all runs, along with the average percentage
% variance explained within the GLM for each voxel. We then normalise these
% timeseries and extract the relevant condition-specific beta values for
% each voxel (by deconvolving the HDR with the timeseries). We then select
% the top 100 voxels (in terms of percentage variance explained for further
% analysis).

% KWN 30/07/2018

clear; close all; % general housekeeping.
ynicInit spm8 % add spm to the path.
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions');

%% ----------------------- PARAMETERS TO EDIT -------------------------- %%

exp_condition = 'colour'; % specify overall project version.
condition = 'feature'; % specify analysis type.

% specify the type of ROIs we wish to analyse (inner_outer vs combined).
ROI_condition = 'combined';

% specify number of voxels to select (highest %age variance explained).
if isequal(ROI_condition, 'combined')
    voxels = 100;
elseif isequal(ROI_condition, 'inner_outer')
    voxels = 70;
end

% specify participant numbers and number of valid runs for each individual
% (this varies whether we are analysing data from the original, or the
% newer colour project).
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
elseif isequal(exp_condition, 'naturalistic')
    participants{1} = 'R2548';   valid_runs{1} = 1:8;
    participants{2} = 'R2904';   valid_runs{2} = 1:8;
    participants{3} = 'R3111';   valid_runs{3} = 1:6;
    participants{4} = 'R3517';   valid_runs{4} = 1:8;
    participants{5} = 'R3773';   valid_runs{5} = 1:7;
    % participants{6} = 'R4059';   valid_runs{6} = 1:7;
    participants{6} = 'R4065';   valid_runs{6} = 1:8;
    participants{7} = 'R4127';   valid_runs{7} = 1:8;
    participants{8} = 'R4244';   valid_runs{8} = 1:7;
    participants{9} = 'R4829';  valid_runs{9} = 1:8;
    % participants{11} = 'R4831';  valid_runs{11} = 1:7;
    participants{10} = 'R4833';  valid_runs{10} = 1:8;
    participants{11} = 'R4890'; valid_runs{11} = 1:8;
    participants{12} = 'R4928'; valid_runs{12} = 1:6;
    participants{13} = 'R5006'; valid_runs{13} = 1:8;
    
    % for the naturalistic analysis, we also specify participant
    % anatomy files.
    anat{1} = '/scratch/groups/Projects/P1361/fMRI/R2548/anatomy/T1.nii.gz';
    anat{2} = '/scratch/groups/Projects/P1361/fMRI/R2904/anatomy/nu.nii.gz';
    anat{3} = '/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/vAnatomy.nii.gz';
    anat{4} = '/scratch/groups/Projects/P1361/fMRI/R3517/anatomy/nu.nii.gz';
    anat{5} = '/scratch/groups/Projects/P1361/fMRI/R3773/anatomy/T1.nii.gz';
    % anat{6} = '/scratch/groups/Projects/P1361/fMRI/R4059/anatomy/T1.nii.gz';
    anat{6} = '/scratch/groups/Projects/P1361/fMRI/R4065/anatomy/t1.nii.gz';
    anat{7} = '/scratch/groups/Projects/P1361/fMRI/R4127/anatomy/T1.nii.gz';
    anat{8} = '/scratch/groups/Projects/P1361/fMRI/R4244/anatomy/t1.nii.gz';
    anat{9} = '/scratch/groups/Projects/P1361/fMRI/R4829/anatomy/t1.nii.gz';
    % anat{11} = '/scratch/groups/Projects/P1361/fMRI/R4831/anatomy/t1.nii.gz';
    anat{10} = '/scratch/groups/Projects/P1361/fMRI/R4833/anatomy/t1.nii.gz';
    anat{11} = '/scratch/groups/Projects/P1361/fMRI/R4890/anatomy/t1.nii.gz';
    anat{12} = '/scratch/groups/Projects/P1361/fMRI/R4928/anatomy/t1.nii.gz';
    anat{13} = '/scratch/groups/Projects/P1361/fMRI/R5006/anatomy/t1.nii.gz';
end

% specify the ROIs we wish to analyse and their corresponding shorter
% 'output names' to use in later output file naming.
if isequal(ROI_condition, 'combined')
    ROI_names{1} = '01_Combined_V1_KWN.mat';   ROI_output{1} = 'V1';
    ROI_names{2} = '02_Combined_V2_KWN.mat';   ROI_output{2} = 'V2';
    ROI_names{3} = '03_Combined_V3_KWN.mat';   ROI_output{3} = 'V3';
    ROI_names{4} = '06_Combined_V3A_KWN.mat';  ROI_output{4} = 'V3A';
    ROI_names{5} = '07_Combined_V3B_KWN.mat';  ROI_output{5} = 'V3B';
    ROI_names{6} = '08_Combined_V4_KWN.mat';   ROI_output{6} = 'V4';
    ROI_names{7} = '09_Combined_LO1_KWN.mat';  ROI_output{7} = 'LO1';
    ROI_names{8} = '10_Combined_LO2_KWN.mat';  ROI_output{8} = 'LO2';
    ROI_names{9} = '11_Combined_A1_KWN.mat';   ROI_output{9} = 'A1';
    ROI_names{10} = '12_Combined_MT+_KWN.mat'; ROI_output{10} = 'MT+';
    ROI_names{11} = '14_Combined_IPS0_KWN.mat';   ROI_output{11} = 'IPS0';
    
elseif isequal(ROI_condition, 'inner_outer')
    ROI_names{1} = '01_Combined_V1_KWN_Inner.mat'; ROI_output{1} = 'V1I';
    ROI_names{2} = '01_Combined_V1_KWN_Outer.mat'; ROI_output{2} = 'V1O';
    ROI_names{3} = '02_Combined_V2_KWN_Inner.mat'; ROI_output{3} = 'V2I';
    ROI_names{4} = '02_Combined_V2_KWN_Outer.mat'; ROI_output{4} = 'V2O';
    ROI_names{5} = '03_Combined_V3_KWN_Inner.mat'; ROI_output{5} = 'V3I';
    ROI_names{6} = '03_Combined_V3_KWN_Outer.mat'; ROI_output{6} = 'V3O';
    ROI_names{7} = '08_Combined_V4_KWN_Inner.mat'; ROI_output{7} = 'V4I';
    ROI_names{8} = '08_Combined_V4_KWN_Outer.mat'; ROI_output{8} = 'V4O';
end

% specify the condition-specific parfile naming structure, some experiments
% name the parfiles with 'scan_' and some 'run_' so we specify this here,
% to be incremented within the script and we also specify the
% condition-specific parfile ending.
if isequal(condition, 'block')
    parfile_ext = {'*scan_0','_interblock&fix-5.par'};
elseif isequal(condition, 'feature') || isequal(condition, 'colour')
    parfile_ext = {'*run_0', '_edited5.par'};
elseif isequal(condition, 'colourxfeature')
    parfile_ext = {'*run_0'};
end

% specify the length of TRs used in the scanning session.
if isequal(exp_condition, 'naturalistic')
    tr = 2;
else
    tr = 3;
end

% get a standard HRF from spm, which we will use to convolve with our
% event sequence and extract betas.
[spmhrf,p] = spm_hrf(tr);

%% ------------------- SPECIFY OUTPUT DIRECTORIES ---------------------- %%

% specify the overarching condition-specific data storage directory.
if isequal(exp_condition, 'naturalistic')
    root_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/%s/',condition);
else
    root_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/', exp_condition, condition);
end

if isequal(ROI_condition, 'combined')
    % specify the analysis-condition output directories.
    timeseries_dir = strcat(root_dir, 'timeseries/'); [~,~] = mkdir(timeseries_dir);
    varexp_dir = strcat(root_dir, 'varexp/'); [~,~] = mkdir(varexp_dir);
    normalised_dir = strcat(root_dir, 'normalised_timeseries/'); [~,~] = mkdir(normalised_dir);
    beta_dir = strcat(root_dir, 'betas/'); [~,~] = mkdir(beta_dir);
    selected_dir = strcat(root_dir, 'betasxvarexp/'); [~,~] = mkdir(selected_dir);
    
elseif isequal(ROI_condition, 'inner_outer')
    
    % specify the analysis-condition output directories.
    timeseries_dir = strcat(root_dir, 'timeseries/inner_outer/'); [~,~] = mkdir(timeseries_dir);
    varexp_dir = strcat(root_dir, 'varexp/inner_outer/'); [~,~] = mkdir(varexp_dir);
    normalised_dir = strcat(root_dir, 'normalised_timeseries/inner_outer/'); [~,~] = mkdir(normalised_dir);
    beta_dir = strcat(root_dir, 'betas/inner_outer/'); [~,~] = mkdir(beta_dir);
    selected_dir = strcat(root_dir, 'betasxvarexp/inner_outer/'); [~,~] = mkdir(selected_dir);
end

if isequal(exp_condition, 'naturalistic')
    base_dir = '/scratch/groups/Projects/P1361/fMRI/';
else
    base_dir = (sprintf('/scratch/groups/Projects/P1323/%s/fMRI/', exp_condition));
end

%% --------------- EXTRACT MULTIVARIATE TIMESERIES DATA ---------------- %%

for ROInumber = 1:length(ROI_names) % for each ROI in turn:
    
    % specify the corresponding full and abbreviated ROI names.
    ROI = ROI_names{ROInumber}; outputROI = ROI_output{ROInumber};
    
    for thisparticipant = 1:length(participants) % for each participant in turn:
        
        participant = participants{thisparticipant}; % specify the participant number.
        
        % for one particular participant, their MT+ ROI is associated with a
        % different high resolution structural. if this is the participant
        % we are currently analysing, and we are looking at MT+ data, set
        % the mrvista directory and ROI locations accordingly.
        if isequal(participant, 'R3111') && isequal(ROI, '12_Combined_MT+_KWN.mat')
            if isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
                mrvista_dir = '/scratch/groups/Projects/P1361/fMRI/R3111/mrvista/MT/';
                ROI_dir = '/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/MT/ROIs/';
                
                % ensure we update the MT+ ROI location for a single
                % participant.
                anat{3} = '/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/MT/nu_fast_restore.nii.gz';
            elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3x3feature')
                mrvista_dir = '/scratch/groups/Projects/P1361/fMRI/R3111/mrvista/multiple/MT/';
                ROI_dir = '/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/MT/ROIs/';
            else
                mrvista_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT',...
                    exp_condition, participant, condition);
                ROI_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/MT/ROIs',...
                    exp_condition, participant);
            end
        elseif isequal(participant, 'R3517') && isequal(exp_condition, 'naturalistic')
            if isequal(condition, '3feature')
                mrvista_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans/';
                ROI_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/anatomy/ROIs/';
            elseif isequal(condition, '3x3feature')
                mrvista_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans/multiple/';
                ROI_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/anatomy/ROIs/';
            end
            
        else % otherwise, in all other cases:
            
            % ensure we reset the ROI location for a single participant.
            if isequal(exp_condition, 'naturalistic')
                if isequal(condition, '3feature')
                    anat{3} = '/scratch/groups/Projects/P1361/fMRI/R3111/anatomy/vAnatomy.nii.gz';
                    
                    % specify the directory containing the mrvista session of interest.
                    mrvista_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/', participant);
                    
                    % specify the ROI storage directory.
                    ROI_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/anatomy/ROIs', participant);
                elseif isequal(condition, '3x3feature')
                    % specify the directory containing the mrvista session of interest.
                    mrvista_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/multiple/', participant);
                    
                    % specify the ROI storage directory.
                    ROI_dir = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/anatomy/ROIs/', participant);
                end
            else
                mrvista_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s',...
                    exp_condition, participant, condition);
                ROI_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs',...
                    exp_condition, participant);
            end
        end
        
        cd(mrvista_dir); % index into the mrvista session directory.
        
        % for one participant, mrvista does not by default select their
        % corresponding high-resolution structural scan, so, for this
        % participant, we specify the location of this file whilst opening
        % a mrvista gray window.
        if isequal(exp_condintion, 'naturalistic')
            VOLUME{1} = mrVista('3view', anat{thisparticipant});
        else
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
        end
        
        % ensure we are working with the motion compensation datatype.
        VOLUME{1} = selectDataType(VOLUME{1},3); VOLUME{1} = setCurScan(VOLUME{1},1);
        
        % load the current ROI of interest into the mrVista volume view, we
        % do this differently for inner/outer ROIs, which are stored in a
        % seperate directory.
        if ~isempty(strfind(ROI, 'Inner')) || ~isempty(strfind(ROI, 'Outer'))
            VOLUME{1} = loadROI(VOLUME{1}, strcat('/inner_outer', '/', ROI));
        else
            VOLUME{1} = loadROI(VOLUME{1}, ROI);
        end
        
        if isequal(exp_condition, 'original')
            nTRs = 131; % specify the number of TRs per scan (after dummy volumes are removed).
        elseif isequal(exp_condition, 'colour')
            nTRs = 99;
        elseif isequal(exp_condition, 'naturalistic')
            nTRs = 150;
        end
        
        % specify the locations of the first TR of the first scan (1) and
        % the last TR of the first scan (nTRs).
        startRow = 1; endRow = nTRs;
        
        if isequal(participant, 'R3111') && isequal(exp_condition, 'original')
            valid_runs{thisparticipant} = 1:5;
        end
        
        ROIcoords = viewGet(VOLUME{1}, 'roicoords');
        
        for thisScan = valid_runs{thisparticipant} % for each scan in turn:
            
            % extract the raw voxel data for the specific ROI and scan.
            % 'false' allows detrending of the data and conversion to % signal.
            [rawVoxData] = extractAllTSeries2(VOLUME{1}, thisScan, ROIcoords, 'false');
            
            % collect this voxel data in a matrix which collates across scans for each participant.
            rawData(startRow:endRow,:) = rawVoxData;
            
            % increment the number of TRs to match the next scan (e.g. scan 1 = 1:99, scan 2 = 100:199).
            startRow = size(rawData,1) + 1; endRow = size(rawData,1) + nTRs;
        end
        
        % save this data across scans in an overall data matrix (extractedData).
        extractedData = rawData;
        
        % remove any voxels outside of the gray matter in which we have GLM
        % data- this is only the case for A1 and MT+ voxels.
        extractedData = extractedData(:,all(~isnan(extractedData)));
        
        mv_plotScans(VOLUME{1},1); % open the visualise multivariate results window.
        
        mv_visualizeGlm; % within this plot window, visualize the GLM results.
        
        tc_dumpDataToWorkspace_univariate_KWN; % dump this GLM voxel data to the workspace.
        
        varexp = uvoxel_data.glm.varianceExplained; % extract the voxel variance explained data.
        
        % specify the ROI-specific multivariate timeseries output directory.
        ROI_timeseries_dir = strcat(timeseries_dir, sprintf('%s', outputROI),'/');
        [~,~] = mkdir(ROI_timeseries_dir);
        
        % save the participant- and ROI-specific multivariate timeseries
        % data within the corresponding output directory.
        save(strcat(ROI_timeseries_dir, sprintf('%s_tSeries_%s.mat', participant, outputROI)),'extractedData');
        
        % add the participant-specific extracted timeseries data to the
        % group storage array.
        groupTSeries{thisparticipant} = extractedData;
        
        % specify the ROI-specific multivariate percentage variance
        % explained output directory.
        ROI_varexp_dir = strcat(varexp_dir, sprintf('%s', outputROI),'/');
        [~,~] = mkdir(ROI_varexp_dir);
        
        % save the participant- and ROI-specific multivariate percentage
        % variance explained data within the corresponding output
        % directory.
        save(strcat(ROI_varexp_dir,sprintf('%s_varexp_%s.mat', participant, outputROI)),'varexp');
        
        % add the participant-specific extracted percentage variance
        % explained data to the group storage array.
        groupvarexp{thisparticipant} = varexp;
        
        %% -------------- NORMALISE MULTIVARIATE TIMESERIES DATA --------------- %%
        
        % specify the locations of the first TR of the first scan (1) and
        % the last TR of the first scan (nTRs).
        startRow = 1; endRow = nTRs;
        
        for thisScan = valid_runs{thisparticipant} % for each scan in turn:
            
            data = extractedData(startRow:endRow,:); % extract a scans-worth of data.
            
            % create a polynomial data column, with second column 1:99, and a column of 1s.
            A = [(1:nTRs).^2;1:nTRs;ones(1,nTRs)].';
            
            % divide the polynomial by the timeseries data for the specific scan.
            Best_Fit = A*(A\data);
            
            data = data-Best_Fit; % remove linear trends (ramp) from the data.
            
            % calculate the stdev of each voxel across all TRs, and tile to
            % repeat value for each TR within a run.
            sdResponse = repmat(std(data),size(data,1),1);
            
            % divide the amplitude for each voxel by the standard deviation
            % across all timepoints (normalise).
            meanResponse = (data)./ sdResponse;
            
            % store this normalised timeseries data in a run-specific location
            % in a storage array.
            normResponse(startRow:endRow,:) = meanResponse;
            
            startRow = startRow + nTRs; endRow = (startRow + nTRs) - 1; % increment number of TRs for next scan.
        end % repeat for next scan.
        
        % specify the ROI-specific multivariate normalised timeseries
        % output directory.
        ROI_normalised_dir = strcat(normalised_dir, outputROI, '/'); [~,~] = mkdir(ROI_normalised_dir);
        
        % save the participant- and ROI-specific multivariate timeseries
        % data within the corresponding output directory.
        save(strcat(ROI_normalised_dir, sprintf('%s_normalised_tSeries_%s.mat',...
            participant, outputROI)), 'normResponse');
        
        % add the participant-specific normalised timeseries data to the
        % group storage array.
        normalisedTSeries{thisparticipant} = normResponse;
        
        clear normResponse
        
        %%  ---------------------- GET EVENT TIMINGS --------------------------- %%
        
        % specify the participant-specific parfile storage directory.
        if isequal(exp_condition, 'naturalistic')
            if isequal(participant, 'R3517')
                if isequal(condition, '3feature')
                    parfile_dir = [base_dir, sprintf('%s/mrvista/session1-allscans/Stimuli/Parfiles/', participant)];
                elseif isequal(condition, '3x3feature')
                    parfile_dir = '/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans/multiple/Stimuli/Parfiles/';
                end
            else
                if isequal(condition, '3feature')
                    parfile_dir = [base_dir, sprintf('%s/mrvista/Stimuli/Parfiles/', participant)];
                elseif isequal(condition, '3x3feature')
                    parfile_dir = [base_dir, sprintf('%s/mrvista/multiple/Stimuli/Parfiles/', participant)];
                end
            end
        else
            parfile_dir = [base_dir, sprintf('%s/mrvista/%s/Stimuli/Parfiles/', participant, condition)];
        end
        
        % specify the locations of the first TR of the first scan (1) and
        % the last TR of the first scan (nTRs).
        startRow = 1; endRow = nTRs;
        
        if isequal(participant, 'R3111') && isequal(exp_condition, 'original')
            valid_runs{thisparticipant} = [1 2 3 4 7];
        end
        
        for thisScan = valid_runs{thisparticipant} % for each scan in turn:
            
            % specify the partcipant- and scan-specific parfile name.
            if isequal(condition, 'event')
                parfile = dir(strcat(parfile_dir, strcat('*scan_0', sprintf('%d', thisScan'), '_3tr.par')));
            elseif isequal(condition, 'colourxfeature')
                parfile = dir(strcat(parfile_dir, strcat(parfile_ext{1}, sprintf('%d.par', thisScan'))));
            elseif isequal(exp_condition, 'naturalistic')
                parfile = dir(strcat(parfile_dir, sprintf('%d*.par', thisScan)));
            else
                parfile = dir(strcat(parfile_dir, strcat(parfile_ext{1}, sprintf('%d', thisScan'), parfile_ext{2})));
            end
            
            % load in this participant- and scan-specific parfile, storing
            % the time column and event column in seperate arrays.
            if isequal(exp_condition, 'naturalistic')
                [timing, event] = textread(strcat(parfile_dir, parfile.name),'%f%d%*[^\n]');
            else
                [timing, event] = textread(strcat(parfile_dir, parfile.name),'%f %d');
            end
            
            % combine the timing and event information into a singular
            % matrix to simultaneously remove interblock information from
            % timing information and event codes.
            combined_data = [timing, event];
            
            % re-seperate these timing and event columns of data for
            % further processing
            timing = combined_data(:,1); event = combined_data(:,2);
            
            timingTRs = timing / tr;   % convert timings from seconds to TRs.
            timingTRs = timingTRs + 1; % add one to the TR timing of events.
            
            eventlist = zeros(nTRs,1); % create a 1 column array of zeros x number of TRs in single scan.
            indexTr = find(timingTRs); % store the positions of timing points in the parfile.
            
            for thisEvent = 1:length(timingTRs) % for each event in turn:
                
                % find the TR corresponding to the event number.
                trToGet = timingTRs(thisEvent);
                
                % add the event number to the corresponding location in
                % this matrix of zeros.
                eventlist(trToGet) = event(indexTr(thisEvent));
            end
            
            % save the padded event-sequence in a scan-specific location
            % in a storage array before recoding.
            eventlistByType(startRow:endRow) = eventlist;
            
            % here we recode events to estimate the haemodynamic response-
            % change all events to 1 (regardless of type).
            eventlist(eventlist==2)=1; % BY/contrast.
            eventlist(eventlist==3)=1; % LUM/shape
            eventlist(eventlist==4)=1; % passive (feature_block only).
            eventlist(eventlist==5)=1; % interblock (colour condition only).
            eventlist(eventlist==6)=1; % mismatched events (event analysis only).
            eventlist(eventlist==7)=1;
            eventlist(eventlist==8)=1;
            eventlist(eventlist==9)=1;
            eventlist(eventlist==10)=1;
            eventlist(eventlist==11)=1;
            eventlist(eventlist==12)=1;
            eventlist(eventlist==13)=1; % for all of the colour-interaction analyses.
            
            % save this modified padded event sequence in a scan-specific
            % location in an overall storage array (eventually we will
            % have one an event sequence for all scans for a participant).
            eventlistAllScans(startRow:endRow) = eventlist;
            
            % increment the number of TRs to match the next scan.
            startRow = size(eventlistAllScans,2) + 1; endRow = size(eventlistAllScans,2) + nTRs;
        end % repeat for next scan.
        
        % transpose this data to a 1 column x number of TRs event sequence.
        allEventTypes = eventlistByType';
        
        fulleventseq{thisparticipant} = allEventTypes; % save this data for each participant.
        
        % refresh variables for use in the next iteration.
        clear data A Best_Fit_ sdResponse meanResponse normResponse rawData startRow endRow
        clear rawVoxData eventlist allEvents indexTr eventlistAllScans eventlistByType
        
        close all % close the mrvista window.
        
    end % repeat process for next participant.
    
    % save the ROI-specifiic multivariate group timeseries, percentage
    % variance explained and normalised timeseries data in the data-specific
    % directories.
    save(strcat(ROI_timeseries_dir, sprintf('group_tSeries_%s.mat', outputROI)), 'groupTSeries');
    save(strcat(ROI_varexp_dir, sprintf('group_varexp_%s.mat', outputROI)), 'groupvarexp');
    save(strcat(ROI_normalised_dir, sprintf('group_normalised_tSeries_%s.mat', outputROI)), 'normalisedTSeries');
    
    %%  -------------------- CALCULATE BETA VALUES ------------------------- %%
    
    for thisparticipant = 1:length(participants) % for each participant in turn:
        
        % get the corresponding participant number.
        participant = participants{thisparticipant};
        
        % write to the command line the participant we are currently processing.
        toPrint = sprintf('Getting event betas for %s', participant); disp(toPrint)
        
        % get the participant-specific unmodified event sequence across all scans.
        eventlistByType = fulleventseq{thisparticipant};
        
        % get the participant-specific normalised timeseries data.
        participantTSeries = normalisedTSeries{thisparticipant};
        
        % get the total number of TRs across all scans.
        timingsAllScans = 1:length(eventlistByType);
        
        % create a design matrix (TRs, events) across all scans.
        fulldesign = [timingsAllScans',eventlistByType];
        
        totalTRs = size(fulldesign,1); % get the total number of TRs across all scans.
        
        allEVs = find(eventlistByType); % get the position of events (TRs when events occured).
        
        % create a row vector 1:number of events across scans.
        allTimes = 1:length(timingsAllScans(allEVs));
        
        % save a column array of the number of TRs across all scans.
        timePoints = timingsAllScans';
        
        % create a zero-array TRs x number of events across all scans.
        event = zeros(totalTRs,length(allTimes));
        
        % for each event in turn, code the event position with a 1.
        for thisEvent = 1:length(allTimes)
            event(allEVs(thisEvent),thisEvent) = 1;
        end
        
        hrfthisROI = spmhrf; % get a copy of the SPM haemodynamic response.
        
        % convolve this haemodynamic response function with the number of
        % TRs across all scans.
        cm = GetConvolutionMatrix(hrfthisROI',totalTRs);
        
        % convolve our individual events with our ROI-specific HDR.
        for thisEv = 1:size(event,2) % for each event in turn:
            
            % multiply this specific event in the design matrix with the
            % convolution matrix, and store in a new convolved design matrix.
            convEvents(:,thisEv) = event(:,thisEv)'*cm;
        end
        
        % add a column matrix of 1s (1:total TRs) to design matrix.
        meanCol = ones(size(timePoints)); designMatrix = [meanCol,convEvents];
        
        % add a 'ramp' column matrix to account for low frequency drift.
        ramp = linspace(0,1,length(timePoints)); designMatrix = [designMatrix,ramp'];
        
        X = designMatrix; % create a copy of our design matrix.
        
        startEV = 1; % specify the starting number of events (always 1).
        
        % if this is feature analysis condition (or the equivalent from
        % the original experiment):
        if isequal(condition, 'block')
            
            % calculate the number of events per type
            % (orientation, contrast, shape, passive).
            nEV = 4*length(valid_runs{thisparticipant});
        elseif isequal(condition, 'feature')
            
            nEV = 3*length(valid_runs{thisparticipant});
            
        elseif isequal(condition, 'colour')
            
            % calculate the number of events per type
            % (red-green, blue-yellow, luminance).
            nEV = 4*length(valid_runs{thisparticipant});
            
        elseif isequal(condition, 'colourxfeature')
            nEV = 1*length(valid_runs{thisparticipant});
        elseif isequal(condition, '3feature')
            nEV = 3*length(valid_runs{thisparticipant});
        elseif isequal(condition, '3x3feature')
            nEV =length(valid_runs{thisparticipant});
        end
        
        endEV = nEV; % specify the total number of events of each type.
        
        if isequal(condition, '3feature')
            endEV = 0;
        end
        
        for thisVoxel = 1:size(participantTSeries,2) % for each voxel in turn:
            
            % extract the data for each TR for this specific voxel.
            Y = double(participantTSeries(:,thisVoxel));
            
            % take the pseudoinverse of the design matrix multipled with
            % the voxel-specific data across all scans.
            betasThisVoxel = inv((X'*X))*X'*Y;
            
            betasEVs = betasThisVoxel(2:end-1); % get a beta for each event for that voxel.
            
            % we now have the betas for an individual voxel for each event,
            % so we'll sort them by event type:
            % extract the original (unpadded) event sequence.
            goodEvent = eventlistByType(eventlistByType>0);
            
            % if this is feature analysis condition (or the equivalent from
            % the original experiment).
            if isequal(condition, 'feature') || isequal(condition, 'block')
                
                % find the position of orientation events and extract the
                % corresponding betas values for this voxel.
                idxO = find(goodEvent==1); betasO = betasEVs(idxO);
                idxC = find(goodEvent==2); betasC = betasEVs(idxC); % contrast
                idxS = find(goodEvent==3); betasS = betasEVs(idxS); % shape
                idxP = find(goodEvent==4); betasP = betasEVs(idxP); % passive
                idxIB = find(goodEvent == 5); betasIB = betasEVs(idxIB); % interblock
                
                betasIB = betasIB(1:nEV);
                
                % concatenate the betas of interest.
                allBetas = [betasO,betasC,betasS,betasP, betasIB];
                
            elseif isequal(condition, 'colour') % if this is colour analysis condition:
                
                % perform the same process as above, but for red-green,
                % blue-yellow, luminance and interblock events.
                idxRG = find(goodEvent==1); betasRG = betasEVs(idxRG);
                idxBY = find(goodEvent==2); betasBY = betasEVs(idxBY);
                idxLUM = find(goodEvent==3); betasLUM = betasEVs(idxLUM);
                idxIB = find(goodEvent==5); betasIB = betasEVs(idxIB);
                
                % extract only the number of betas for the interblock
                % condition matching the number of betas extracted for the
                % three colour conditions, to make writing out the betas
                % matrix easier.
                betasIB = betasIB(1:nEV);
                
                % concatenate the betas of interest.
                allBetas = [betasRG,betasBY,betasLUM,betasIB];
                
            elseif isequal(condition, 'event')
                % find the position of orientation events and extract the
                % corresponding betas values for this voxel.
                idxO = find(goodEvent==1); betasO = betasEVs(idxO);
                idxC = find(goodEvent==2); betasC = betasEVs(idxC); % repeat for contrast events.
                idxS = find(goodEvent==3); betasS = betasEVs(idxS); % shape events.
                idxNC = find(goodEvent==4); betasNC = betasEVs(idxNC); % no change events.
                idxMM = find(goodEvent==6); betasMM = betasEVs(idxMM); % and mismatch events.
                
                % find the length of the largest array of betas.
                maxlength = max([length(betasO), length(betasC), length(betasS), length(betasNC), length(betasMM)]);
                
                % set the number of events equal to the length of the condition with the most events.
                nEV = maxlength;
                
                % if this is the first voxel processed for this participant,
                % set the end event equal to the maximum number of events.
                if thisVoxel==1; endEV = nEV; end
                
                % create an array of nans equal to the length of the longest
                % beta array.
                padBetas = zeros(1,maxlength)'*NaN;
                
                % pad each of the extracted beta arrays with nans so that they
                % are all equal to the greatest number of events.
                padBetasO = padBetas; padBetasO(1:length(betasO)) = betasO;
                padBetasC = padBetas; padBetasC(1:length(betasC)) = betasC;
                padBetasS = padBetas; padBetasS(1:length(betasS)) = betasS;
                padBetasNC = padBetas; padBetasNC(1:length(betasNC)) = betasNC;
                padBetasMM = padBetas; padBetasMM(1:length(betasMM)) = betasMM;
                
                % concatenate the betas of interest.
                allBetas = [padBetasO,padBetasC,padBetasS,padBetasNC,padBetasMM];
            elseif isequal(condition, 'colourxfeature')
                
                idxRGO = find(goodEvent==1); betasRGO = betasEVs(idxRGO);
                idxRGC = find(goodEvent==2); betasRGC = betasEVs(idxRGC);
                idxRGS = find(goodEvent==3); betasRGS = betasEVs(idxRGS);
                idxRGP = find(goodEvent==4); betasRGP = betasEVs(idxRGP);
                idxBYO = find(goodEvent==5); betasBYO = betasEVs(idxBYO);
                idxBYC = find(goodEvent==6); betasBYC = betasEVs(idxBYC);
                idxBYS = find(goodEvent==7); betasBYS = betasEVs(idxBYS);
                idxBYP = find(goodEvent==8); betasBYP = betasEVs(idxBYP);
                idxLUMO = find(goodEvent==9); betasLUMO = betasEVs(idxLUMO);
                idxLUMC = find(goodEvent==10); betasLUMC = betasEVs(idxLUMC);
                idxLUMS = find(goodEvent==11); betasLUMS = betasEVs(idxLUMS);
                idxLUMP = find(goodEvent==12); betasLUMP = betasEVs(idxLUMP);
                idxIB = find(goodEvent == 13); betasIB = betasEVs(idxIB);
                
                betasIB = betasIB(1:nEV);
                
                % concatenate the betas of interest.
                allBetas = [betasRGO, betasRGC, betasRGS, betasRGP, betasBYO,...
                    betasBYC, betasBYS, betasBYP, betasLUMO, betasLUMC,betasLUMS, betasLUMP,betasIB];
                
            elseif isequal(condition, '3feature')
                idxF = find(goodEvent==1); betasF = betasEVs(idxF);
                idxO = find(goodEvent==2); betasO = betasEVs(idxO);
                idxC = find(goodEvent==3); betasC = betasEVs(idxC);
                idxS = find(goodEvent==4); betasS = betasEVs(idxS);
                idxP = find(goodEvent == 5); betasP = betasEVs(idxP);
                
                maxlength = 24;
                
                betasO(length(betasO)+1:24) = NaN;
                betasC(length(betasC)+1:24) = NaN;
                betasS(length(betasS)+1:24) = NaN;
                betasF(length(betasF)+1:24) = NaN;
                betasP(length(betasP)+1:24) = NaN;
                
                allBetas = [betasO, betasC, betasS, betasF, betasP];
                
                if thisVoxel == 1
                    endEV = maxlength;
                end
                
                nEV = maxlength;
                
            elseif isequal(condition, '3x3feature')
                idxP = find(goodEvent==1); betasP = betasEVs(idxP);
                idxV = find(goodEvent==2); betasV = betasEVs(idxV);
                idxH = find(goodEvent==3); betasH = betasEVs(idxH);
                idxD = find(goodEvent==4); betasD = betasEVs(idxD);
                idxR = find(goodEvent==5); betasR = betasEVs(idxR);
                idxG = find(goodEvent==6); betasG = betasEVs(idxG);
                idxB = find(goodEvent==7); betasB = betasEVs(idxB);
                idxC = find(goodEvent==8); betasC = betasEVs(idxC);
                idxS = find(goodEvent==9); betasS = betasEVs(idxS);
                idxT = find(goodEvent==10); betasT = betasEVs(idxT);
                idxF = find(goodEvent==11); betasF = betasEVs(idxF);
                
                maxlength = length(valid_runs{thisparticipant});
                
                betasP(end+1:maxlength) =NaN;
                betasV(end+1:maxlength) =NaN;
                betasH(end+1:maxlength) =NaN;
                betasD(end+1:maxlength) =NaN;
                betasR(end+1:maxlength) =NaN;
                betasG(end+1:maxlength) =NaN;
                betasB(end+1:maxlength) =NaN;
                betasC(end+1:maxlength) =NaN;
                betasS(end+1:maxlength) =NaN;
                betasT(end+1:maxlength) =NaN;
                betasF(end+1:maxlength) =NaN;
                
                
                % concatenate the betas of interest.
                allBetas = [betasP, betasV, betasH, betasD, betasR, betasG, betasB, betasC, betasS, betasT,betasF];
                
            end
            
            % stack these beta values for each event- one column per event.
            betasAllVoxels(startEV:endEV,:) = allBetas;
            
            % increment event values for next loop iteration.
            startEV = length(betasAllVoxels) + 1; endEV = (startEV + nEV) - 1;
            
            % add the beta values for this voxel to a 3D matrix: repetion x condition x voxels.
            betaList(:,:,thisVoxel) = allBetas;
            
            clear allBetas Y goodEvent betasAllVoxels
        end
        
        ROI_beta_dir = strcat(beta_dir, outputROI, '/'); [~,~] = mkdir(ROI_beta_dir);
        
        save(strcat(ROI_beta_dir, sprintf('%s_extractedBetas_%s.mat',...
            participant, outputROI)),'betaList'); % save this extracted beta data.
        
        % add these participant-specific extracted betas to a group storage array.
        extractedBetasAllPs{thisparticipant} = betaList;
        
        % refresh some important variables for the next participant.
        clear convEvents betaList allBetas timingsAllScans eventlistByType
    end % repeat process for the next participant.
    
    % save extracted betas for all participants.
    save(strcat(ROI_beta_dir, sprintf('group_extractedBetas_%s.mat', outputROI)), 'extractedBetasAllPs');
    
    close all % close the mrvista windows.
end % continue to the next ROI.

%% -------- SELECT VOXELS BY TOP PERCENTAGE VARIANCE EXPLAINED --------- %%

for ROInum = 1:length(ROI_names) % for each ROI in turn:
    
    outputROI = ROI_output{ROInum}; % extract the ROI-specific output name.
    
    for thisparticipant = 1:length(participants) % for each participant in turn:
        
        participant = participants{thisparticipant}; % extract the participant number.
        
        % load the participant and ROI-specific beta and variance explained files.
        load(strcat(beta_dir, sprintf('%s/', outputROI), sprintf('%s_extractedBetas_%s.mat', participant, outputROI)));
        load(strcat(varexp_dir, sprintf('%s/', outputROI), sprintf('%s_varexp_%s.mat', participant, outputROI)));
        
        % sort the variance explained in descending order, and also sort
        % the indicies of each data point (i.e. the voxel number the data
        % has come from) in the same way.
        [sorted_varexp_vals, sorted_varexp_ind] = sort(varexp,'descend');
        
        % select the number of voxels specified to retain with the highest
        % percentage variances, and again, extract the indicies of each of
        % these data points.
        sorted_varexp_values = sorted_varexp_vals(1:voxels);
        sorted_varexp_indicies = sorted_varexp_ind(1:voxels);
        
        % select the beta values corresponding to the voxels with the top n
        % percentage variance explained.
        final_betas = betaList(:,:, sorted_varexp_indicies);
        
        % store these participant-specific extracted betas in a group array.
        group_betas{thisparticipant} = final_betas;
        
        % create the ROI-specific output storage directory.
        ROI_selected_dir = strcat(selected_dir, outputROI, '/', sprintf('%d', voxels), '/'); [~,~] = mkdir(ROI_selected_dir);
        
        % save the participant- and ROI- specific selected voxels data.
        save(strcat(ROI_selected_dir, sprintf('%s_extractedBetas_%s.mat', participant, outputROI)), 'final_betas');
        
    end % continue to the next participant.
    
    % save the group ROI-specific selected voxel data.
    save(strcat(ROI_selected_dir, sprintf('group_extractedBetas_%s.mat', outputROI)), 'group_betas');
    close all; % close any remaining open mrvista windows.
end % continue to the next ROI.

%% --------------------------------------------------------------------- %%