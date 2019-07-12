%% ---- PLOT DENDROGRAMS OF CLUSTERING ACROSS 3X3 FEATURE CONDITIONS --- %%

% loads in naturalistic 3x3 univariate betas and produces a dendrogram
% illustrating the hierarchical clustering amongst conditions for specified
% visual ROIs.

% 12/7/19 KWN

clear; clc; close all;

%% ------------------- SPECIFY ANALYSIS PARAMETERS --------------------- %%

% specify ROIs to analyse.
ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2'};

% specify overall condition file names.
conditions = {'orientation', 'colour', 'shape', 'passive_face'};

% specify individual features to analyse.
feature_conditions = {'Vertical', 'Horizontal', 'Diagonal', 'Red', 'Green',...
    'Blue', 'Circular', 'Square', 'Triangular', 'Passive', 'Face'};

% specify figure output directory.
outputdir = '/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3x3feature/cluster_analysis/';

%% ------------------------ CREATE DENDROGRAMS ------------------------- %%

for thisROI = 1:length(ROIs) % for each ROI in turn:
    ROI = ROIs{thisROI};
    
    for thiscond = 1:length(conditions) % for each overarching condition in turn:
        condition = conditions{thiscond};
        
        % load in the univariate ROI-specific data (across all participants).
        data{thiscond} = load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3x3feature/%s3x3/group/%s_%s3x3_univariate_betas.mat',...
            condition, ROI, condition));
    end
    
    % combine this data to create a participants x condition matrix.
    catdata = [data{1}.betas_pp, data{2}.betas_pp, data{3}.betas_pp, data{4}.betas_pp];
    
    % produce a dendogram detailing the linkage between these conditions.
    tree = linkage(catdata','average');
    fig = figure();
    [~, ~, outperm] = dendrogram(tree,0);
    
    % label the x-axis with the correct order of conditions for the
    % dendrogram's most efficient presentation.
    for i = 1:length(outperm)
        newconditions{i} = feature_conditions{outperm(i)};
    end
    
    % specify figure formatting.
    xticklabels(newconditions); xtickangle(45);
    title(ROI);
    ylim([0 2.5]);
    
    % save the figure as a .png file.
    saveas(fig, strcat(outputdir, sprintf('%s_3x3feature_dendrogram_univariate.png', ROI)));
    close;
end

%% --------------------------------------------------------------------- %%