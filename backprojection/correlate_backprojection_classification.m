%% ------ CORRELATE EUCLIDEAN DISTANCE BETWEEN BACKPROJECTION MAPS ----- %%
%% -------------------- WITH CLASSIFICATION SCORES --------------------- %%

% loads in feature-specific backprojection maps and calculates euclidean
% (RMSE) distance between them, and we correlate this distance with
% classification scores to examine the relationship between spatial pattern
% of attentional focus and classification. 

% 10/7/19 KWN

clear; clc;

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

% specify ROI to analyse and corresponding number. 
ROI = 'V4';
ROInum = 4;

exp_cond = 'original'; % specify overall experiment to analyse. 
cond = 'block'; % specify condition to analyse. 

% specify data directory,
data_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/backprojection/data/%s/interpdata/',...
    exp_cond, cond, ROI);

% specify conditions to analyse. 
if isequal(cond, 'block') || isequal(cond, 'feature')
    conditions = {'orientation', 'contrast', 'shape'};
elseif isequal(cond, 'colour')
    conditions = {'red-green','blue-yellow','luminance'};
end

for thiscond = 1:length(conditions) % for each condition in turn:
    
    % extract the participant filenames. 
    files = dir(strcat(data_dir,conditions{thiscond},'/*.mat'));
    filenames = {files.name};
    
    % for each participant, load in and store the interpolated data. 
    for thispp = 1:length(filenames)
        load(strcat(data_dir,conditions{thiscond},'/',filenames{thispp}));
        voxeldata(thiscond,thispp,:,:) = interpdata;
    end
end

%% ----- CALCULATE RMSE DISTANCE AND CORRELATE WITH CLASSIFICATION ----- %% 

if isequal(exp_cond, 'original')
    load('/scratch/groups/Projects/P1323/original/fMRI/vista_output/output/multivariate_block_backprojection_twoway.mat');
elseif isequal(exp_cond, 'colour')
    if isequal(cond, 'colour')
        load('/scratch/groups/Projects/P1323/colour/fMRI/vista_output/output/colour/multivariate_colour_backprojection_twoway.mat');
    elseif isequal(cond, 'feature')
        load('/scratch/groups/Projects/P1323/colour/fMRI/vista_output/output/feature/multivariate_feature_backprojection_twoway.mat');
    end
end

% specify numeric comparisons between conditions to analyse. 
cond_comparisons = {[1,2],[1,3],[2,3]};

for thiscomp = 1:size(cond_comparisons,2) % for each comparison in turn:
    for thispp = 1:size(voxeldata,2) % for each participant in turn:
        
        % extract the interpolation data corresponding to these two
        % conditions and calculate the RMSE distance between them.
        data1 = squeeze(voxeldata(cond_comparisons{thiscomp}(1),thispp,:,:));
        data2 = squeeze(voxeldata(cond_comparisons{thiscomp}(2),thispp,:,:));
        
        rmse(thiscomp,thispp) = sum(sum((data1-data2).^2));
    end
    
    % this is currently to remove the participants who don't have a good
    % rendering of the full visual field in V4 (need to improve on this and
    % make condition-specific).
    if isequal(ROI,'V4')
        rmse(thiscomp,3) = nanmean(rmse(thiscomp,:));
        rmse(thiscomp,11) = nanmean(rmse(thiscomp,:));
        rmse(thiscomp,12) = nanmean(rmse(thiscomp,:));
    end
    
    % extract the condition-specific two-way classification data.
    classdata(thiscomp,:) = multivariate_block_twoway.data{ROInum}(:,thiscomp);
    
    % perform shapiro-wilk normality tests on the rmse and classification
    % data.
    [rmse_sw.h(thiscomp),rmse_sw.p(thiscomp),rmse_sw.w(thiscomp)] = swtest(rmse(thiscomp,:));
    [classdata_sw.h(thiscomp),classdata_sw.p(thiscomp),classdata_sw.w(thiscomp)] = swtest(classdata(thiscomp,:));
    
    % correlate rmse and classification values. 
    [corrval.r(thiscomp),corrval.p(thiscomp)]  = corr(rmse(thiscomp,:)',classdata(thiscomp,:)','type', 'Spearman');
end

% write out important analysis output to a .mat file. 
backprojection_output.sw.rmse = rmse_sw;
backprojection_output.sw.classacc = classdata_sw;
backprojection_output.rmse = rmse;
backprojection_output.classacc = classdata;
backprojection_output.corr = corrval;

save(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/output/%s_backprojection_classification_corr.mat',...
    exp_cond, ROI), 'backprojection_output');

%% --------------------------------------------------------------------- %%