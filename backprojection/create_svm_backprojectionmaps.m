%%  CREATE MRVISTA PARAMETER MAPS OF CONDITION-SPECIFIC SUPPORT VECTORS  %%

% For specified ROIs, we create condition-specific parameter maps
% reflecting the mean support vector weights for each voxel (for example,
% we create a map indicating each voxels 'orientation preference'). We
% back-project these voxel-specific support vector weights onto the
% corresponding voxels in a default mrvista parameter map.

% 07/08/2018 KWN

clear; close all; clc;

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

exp_condition = 'colour'; % specify overall experiment we are analysing.
condition = 'feature'; % specify the particular type of data to analyse.

% specify a list of ROIs to analyse and their corresponding ROI number.
ROIs{1} = 'V1';   ROIs_number{1} = '01';
ROIs{2} = 'V2';   ROIs_number{2} = '02';
ROIs{3} = 'V3';   ROIs_number{3} = '03';
ROIs{4} = 'V4';   ROIs_number{4} = '08';

% specify a list of the participants data held within the group betas
% dataset (i.e. the participants with valid data we tested in this
% experiment).
if isequal(exp_condition, 'original')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455',...
        'R3517', 'R3773', 'R3932', 'R4065', 'R4127','R4496'};
elseif isequal(exp_condition, 'colour')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3517',...
        'R3773', 'R4059', 'R4065', 'R4244', 'R4831', 'R4928'};
elseif isequal(exp_condition, 'naturalistic')
    participants = {'R2548', 'R2904', 'R3111', 'R3517', 'R4059','R4127',...
        'R4244', 'R4831','R4833', 'R4890', 'R4928', 'R5006'};
end

% specify the default parameter map we wish to load in (we just use this as
% a template for the total number of voxels for each participant).
if isequal(condition, 'block') || isequal(condition, 'feature')
    default_parametermap = 'OvP';
    
    % specify the columns of data corresponding to the conditions we wish to
    % analyse (currently we only specify one column of data as we are just
    % comparing activation to the passive condition), and specify the
    % corresponding names of these conditions.
    condition_numbers= [1, 2, 3];
    condition_names = {'Orientation', 'Contrast', 'Shape'};
    baseline_condition = 4; % specify number of baseline condition (passive).
    
elseif isequal(condition, 'colour-passive')
    default_parametermap = 'RGAvsRGP';
    
    condition_numbers = 1; % 1, 3 or 5
    condition_names = {'BYA'}; % BYA, RGA or LUMA
    baseline_condition = 2; % 2, 4 or 6
    
elseif isequal(condition, 'colourxfeature')
    default_parametermap = 'RGOvsRGP';
    
    condition_numbers = [1 2 3]; % [1 2 3], [5 6 7], [8 9 10]
    condition_names = {'BlueYellow-Orientation', 'BlueYellow-Contrast', 'BlueYellow-Shape'};
    
    baseline_condition = 4; % 4, 8, 11
    
elseif isequal(exp_condition, 'naturalistic')
    default_parametermap = 'OvsP';
    
    condition_numbers= [1,2,3,4];
    condition_names = {'Orientation', 'Contrast', 'Shape', 'Face'};
    baseline_condition = 5;
end

%% --------- EXTRACT SUPPORT VECTORS AND CREATE PARAMETER MAPS --------- %%

for ROIanalyse = 1:length(ROIs) % for each ROI in turn:
    
    % specify the specific ROI name and corresponding number.
    ROI = ROIs{ROIanalyse};
    ROI_number = ROIs_number{ROIanalyse};
    
    % specify the ROI-specific output beta storage directory, and load the
    % corresponding group betas .mat file.
    if ~isequal(exp_condition, 'naturalistic')
        data_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/betas/%s/',...
            exp_condition, condition, ROI);
    else
        data_directory = sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/betas/%s/', ROI);
    end
    load(strcat(data_directory, sprintf('group_extractedBetas_%s.mat', ROI)));
    
    for participant = 1:length(extractedBetasAllPs) % for each participant in turn:
        
        % specify the corresponding participant number.
        pp_number = participants{participant};
        
        % extract the participant-specific beta value data (number of
        % condition repetitions x number of conditions x number of voxels
        % within the ROI of interest).
        pp_data = extractedBetasAllPs{participant};
        
        for conditionanalyse = 1:length(condition_numbers) % for each specific-condition analysis in turn:
            
            if isequal(condition, 'colour-passive') || isequal(exp_condition, 'naturalistic')
                data1 = squeeze(pp_data(:,conditionanalyse,:));
                data2 = squeeze(pp_data(:,baseline_condition,:));
                
                % remove any nan values from the colour and colour-passive
                % data columns.
                data1(~any(~isnan(data1), 2),:)=[];
                data2(~any(~isnan(data2), 2),:)=[];
                
                % vertically stack and z-score the data.
                conditiondata = [data1;data2];
                conditionzscore = zscore(conditiondata,[],2);
                
                % specify condition labels matching the organisation of
                % columns and calculate weights.
                labels = [ones(size(data1,1),1); ones(size(data2,1),1)*2];
                conditionsizes = [size(data1,1), size(data2,1)];
                [maxsize, ind] = max(conditionsizes);
                
                weights = 1:2;
                weights(ind) = [];
                
                output.weights{ind} = 1;
                output.weights{weights} = maxsize/size(conddata{weights},1);
                
                % run the weighted svm classifier.
                modelSV=libsvmtrain(labels,conditionzscore,...
                    sprintf('-t 2 -w1 %.2f -w2 %.2f -q', output.weights{1}, output.weights{2}));
            else
                
                % extract the two columns of data corresponding to the
                % conditions we wish to analyse (for now we are looking at an
                % attention condition (e.g. orientation), versus passive).
                condition_data = pp_data(:,[condition_numbers(conditionanalyse) baseline_condition],:);
                
                % record the size of this matrix- i.e. number of condition
                % repetitions, number of conditions to analyse and number of
                % voxels.
                [reps,conditions,voxels]=size(condition_data);
                
                % create a 1D column of labels matching the organisation of the
                % data column (i.e. orientation coded as '1', passive as '-1',
                % etc).
                labels = [ones(reps,1);ones(reps,1)*-1];
                
                % vertically stack the condition-specific data to produce a 2D
                % matrix (number of repetitions across all conditions x number of
                % voxels).
                traindata=reshape(condition_data,[reps*conditions,voxels]);
                
                % normalise this data by z-scoring.
                traindata=zscore(traindata,[],2);
                
                % run the svm classification on this normalised data with the
                % corresponding labels. we use a radial basis function kernel, and
                % with the number of cross-validations corresponding to the number
                % of condition repetitions for the individual participant.
                svm_parameters = sprintf('-t 2, -v %d -q', reps);
                classificationaccuracy=libsvmtrain(labels,traindata,svm_parameters);
                
                % we then rerun this classification without the crossvalidation
                % to extract the support vectors associated with this
                % classification.
                modelSV=libsvmtrain(labels,traindata,'-t 2 -q');
            end
            
            % calculate the weighted mean for each data point- not every condition
            % repetition data point contributes 'enough' to be a support vector.
            for thisSV=1:(sum(modelSV.nSV))
                wMean(thisSV,:)=full(modelSV.sv_coef(thisSV)*double(modelSV.SVs(thisSV,:)));
            end
            
            % get the number of support vectors contributing to the first condition.
            svInClass1=modelSV.nSV(1);
            
            % get an index of each voxels 'preference'- we calculate the
            % mean across the support vectors contributing to the first condition.
            meanC1=full(nanmean(wMean(1:svInClass1,:)));
            
            % if we are analysing the one participant who has a seperate
            % mrvista session for the MT+ ROI:
            if isequal(pp_number, 'R3111') && isequal(ROI, 'MT+')
                
                % index into the participant-specific MT+ mrvista directory and
                % open a gray window.
                if isequal(exp_condition, 'naturalistic')
                    cd(sprintf('/scratch/groups/Projects/P1361/fMRI/R3111/mrvista/MT'));
                else
                    cd(sprintf('/scratch/sg3/P1323/%s/fMRI/%s/mrvista/%s/MT',...
                        exp_condition, pp_number, condition));
                end
                VOLUME{1} = mrVista('3view');
                
            elseif isequal(pp_number, 'R3517') && isequal(exp_condition, 'naturalistic')
                cd('/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans');
                
            else % otherwise, for any other ROI and participant:
                
                % index into the participant-specific mrvista directory and
                % open a gray window.
                if isequal(exp_condition, 'naturalistic')
                    cd(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista', pp_number));
                else
                    cd(sprintf('/scratch/sg3/P1323/%s/fMRI/%s/mrvista/%s',...
                        exp_condition, pp_number, condition));
                end
                VOLUME{1} = mrVista('3view');
            end
            
            % load in the ROI we wish to analyse.
            VOLUME{1} = loadROI(VOLUME{1}, sprintf('%s_Combined_%s_KWN.mat', ROI_number, ROI));
            
            VOLUME{1} = selectDataType(VOLUME{1},4);
            VOLUME{1} = setCurScan(VOLUME{1},1);
            VOLUME{1} = refreshScreen(VOLUME{1});
            
            % if we are analysing the one participant who has a seperate
            % mrvista session for the MT+ ROI:
            if isequal(pp_number, 'R3111') && isequal(ROI, 'MT+')
                
                % also load in a parameter map (one within the MT+ specific
                % directory for this participant)- this is arbitrary, we
                % just use this as a template to update with our support
                % vector data.
                if isequal(exp_condition, 'naturalistic')
                    load(sprintf('/scratch/groups/Projects/P1361/fMRI/R3111/mrvista/MT/Gray/GLMs/%s.mat', default_parametermap));
                else
                    load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT/Gray/GLMs/original/%s.mat',...
                        exp_condition, pp_number, condition, default_parametermap));
                end
                
            elseif isequal(pp_number, 'R3517') && isequal(exp_condition, 'naturalistic')
                load(sprintf('/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans/Gray/GLMs/%s.mat',default_parametermap));
                
            else % otherwise, for any other ROI and participant:
                
                % also load in a parameter map- this is arbitrary, we just use
                % this as a template to update with our support vector data.
                if isequal(exp_condition, 'naturalistic')
                    load(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/Gray/GLMs/%s.mat', pp_number, default_parametermap));
                else
                    load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/original/%s.mat',...
                        exp_condition, pp_number, condition, default_parametermap));
                end
            end
            
            VOLUME{1} = refreshScreen(VOLUME{1});
            
            % specify the condition- and ROI-specific name of the
            % parameter map we wish to create
            mapName = sprintf('%s-%s',condition_names{conditionanalyse},ROI);
            
            % extract the co-ordinates of the ROI we are currently
            % analysing (i.e. which voxels (of those across the whole
            % brain) are part of the ROI we are currently interested in).
            roipos = viewGet(VOLUME{1}, 'roiindices');
            
            % create a matrix the size of the number of voxels across the
            % participants whole brain.
            allcoords = 1:length(co{1});
            
            % remove the datapoints corresponding to the positions of
            % voxels contained within the current ROI.
            allcoords(roipos) = [];
            
            % extract the whole-brain data from this current arbitrary
            % parameter map.
            mapcoords = map{1};
            
            % replace the data in this map corresponding to the positions
            % of voxels within the ROI we are analysing with the condition-
            % specific support vectors.
            mapcoords(roipos) = meanC1;
            
            % make any other voxels data 'missing'.
            mapcoords(allcoords)= 0;
            
            % store this updated parameter map as a cell array (needed to
            % match the default format).
            map = {mapcoords};
            
            % if we are analysing the one participant who has a seperate
            % mrvista session for the MT+ ROI:
            if isequal(pp_number, 'R3111') && isequal(ROI, 'MT+')
                
                % specify and create a support-vector specific parameter
                % map output directory.
                if isequal(exp_condition, 'naturalistic')
                    output_directory = sprintf('/scratch/groups/Projects/P1361/fMRI/R3111/mrvista/MT/Gray/GLMs/supportvector/');
                else
                    output_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT/Gray/GLMs/supportvector/',...
                        exp_condition, pp_number, condition);
                end
                [~,~] = mkdir(output_directory);
                
                % save this new parameter map (within the specific MT+
                % directory), with it's corresponding name, and arbitrary
                % units, and with the coherence thresholds from
                % the original analysis (we never threshold by this, so again,
                % it doesn't really matter).
                save(strcat(output_directory,sprintf('%s', pp_number, mapName)), 'map', 'co', 'mapName', 'mapUnits');
                
            elseif isequal(pp_number, 'R3517') && isequal(exp_condition, 'naturalistic')
                
                output_directory = ('/scratch/groups/Projects/P1361/fMRI/R3517/mrvista/session1-allscans/Gray/GLMs/supportvector/');
                [~,~] = mkdir(output_directory);
                
                save(strcat(output_directory,sprintf('%s', pp_number, mapName)), 'map', 'co', 'mapName', 'mapUnits');
            else % otherwise, for any other ROI and participant:
                
                % specify and create a support-vector specific parameter
                % map output directory.
                if isequal(analysis_condition, 'naturalistic')
                    output_directory = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/Gray/GLMs/supportvector/', pp_number);
                else
                output_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Gray/GLMs/supportvector/',...
                    exp_condition, pp_number, condition);
                end
                [~,~] = mkdir(output_directory);
                
                % save this new parameter map, with it's corresponding
                % name, and arbitrary units, and with the coherence
                % thresholds from the original analysis (we never
                % threshold by this, so again, it doesn't really matter).
                save(strcat(output_directory,sprintf('%s', pp_number, mapName)), 'map', 'co', 'mapName', 'mapUnits');
            end
            
            % refresh some variables for the next loop iteration.
            clear clear modelSV wMean meanC1 svInClass1 roipos allcoords mapcoords map
            close all
        end % repeat for the next condition analysis.
        
        close all % close all the mrvista windows for this participant.
    end % repeat this process for the next participant.
end % repeat this process for the next ROI.

%% --------------------------------------------------------------------- %%