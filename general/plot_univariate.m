%% -------------------- PLOT FMRI UNIVARIATE RESULTS ------------------- %%

% produces a variety of plots to summarise our univariate findings, beyond
% what are output from the univariate analysis scripts.

% KWN 9/7/2019

clear; clc;
addpath('/Users/kirstie/Documents/code/fmri_analysis/functions/');

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

% specify ROIs to analyse. 
ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2', 'A1'};

% specify experimental condition to analyse, including interaction if
% examining colourxfeature experiment findings, otherwise, leave this blank. 
% also specify whether we wish to perform parametric or non-parametric
% analyses.
exp_condition = 'colour';
condition = 'colourxfeature';
interaction = 'luminance_feature';
analysis_condition = 'non-parametric';

% specify the location of experiment-specific data.
if isequal(exp_condition, 'original')
    data.directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/group/betas/',...
        exp_condition, condition);
    output.directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/processed',...
        exp_condition);
elseif isequal(exp_condition, 'colour') && (isequal(condition, 'feature') || isequal(condition,'colour'))
    data.directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/group/',...
        exp_condition, condition);
    output.directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/processed',...
        exp_condition, condition);
elseif isequal(condition, 'colourxfeature')
    data.directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/%s/group/',...
        exp_condition, condition, interaction);
    output.directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/%s/processed',...
        exp_condition, condition, interaction);
elseif isequal(condition, 'naturalistic')
    data.directory = sprintf('/Users/kirstie/Documents/analysis/P1361naturalistic/univariate/raw/betas/');
    output.directory = sprintf('/Users/kirstie/Documents/analysis/P1361naturalistic/univariate/processed/');
end

%% ---------------------- ATTENTION VERSUS PASSIVE --------------------- %%

for thisROI = 1:length(ROIs) % for each ROI in turn:
    
    % load in the ROI-specific univariate beta values. 
    data.values = load(strcat(data.directory, sprintf('%s_%s_univariate_betas.mat', ROIs{thisROI}, condition)));
    
    % calculate the mean attention value across all attention conditions 
    % and extract the univariate passive beta for each participant.
    if isequal(condition, 'naturalistic')
        data.attention_all(:,thisROI) =  mean(data.values.betas_pp(:,1:4),2);
        data.passive_all(:,thisROI) = data.values.betas_pp(:,5);
    else
        data.attention_all(:,thisROI) =  mean(data.values.betas_pp(:,1:3),2);
        data.passive_all(:,thisROI) = data.values.betas_pp(:,4);
    end
    
    % run a wilcoxon signed-rank test to compare the attention and passive
    % conditions. 
    if isequal(analysis_condition, 'parametric')
        [data.h(thisROI),data.p(thisROI),data.ci(thisROI,:),data.stats] = ttest(data.attention, data.passive);
        data.tstat(thisROI) = data.stats.tstat;
    elseif isequal(analysis_condition, 'non-parametric')
        [data.p(thisROI),data.h(thisROI),data.stats] = signrank(data.attention, data.passive, 'method', 'approximate');
        data.zstat(thisROI) = data.stats.zval;
    end
    
    % extract the mean of the attention and passive univariate betas.
    data.attendmean(thisROI) = mean(data.attention);
    data.passivemean(thisROI) = mean(data.passive);  
    
    % calculate the difference between the attention and passive condition
    % BOLD signal modulation and calculate the mean and standard deviation.
    data.plotdata.full(thisROI,:) = data.attention-data.passive;
    data.plotdata.mean(thisROI) = mean(data.plotdata.full(thisROI,:));
    data.plotdata.stderr(thisROI) = (std(data.plotdata.full(thisROI,:)) / sqrt(length(data.plotdata.full(thisROI,:))));
end

% perform benjamini-hochberg false disovery rate correction.
[~, ~, ~, adj_p]=fdr_bh(data.p, 0.05, 'dep');

% create a table of data for output.
if isequal(analysis_condition, 'parametric')
    data.output = [data.attendmean; data.passivemean; data.tstat; data.p]';
    data.rowheaders = {'ROIs', 'attendmean', 'passivemean', 'tstat', 'pvalue'};
elseif isequal(analysis_condition, 'non-parametric')
    data.output = [data.attendmean; data.passivemean; data.zstat; data.p; adj_p]';
    data.rowheaders = {'ROIs', 'attendmean', 'passivemean', 'zstat', 'pvalue', 'adj_pvalue'};
end

% add column and row-headers to the data table. 
data.output = [ROIs', num2cell(data.output)];
data.output = [data.rowheaders; data.output];

% write this output data to a .csv file.
if isequal(analysis_condition, 'parametric')
    cell2csv(strcat(output.directory,'/attendvspassivettest.csv'), data.output);
elseif isequal(analysis_condition, 'non-parametric')
    cell2csv(strcat(output.directory,'/attendvspassive-wsrtest.csv'), data.output);
end

% plot data as a bar figure and save this to a .pdf file.
data.plotdata.colour = [.5 .5 .5];
data.pvalthreshold= [0.05, 0.01, 0.001];
figure(1)
superbar(data.plotdata.mean, 'E', data.plotdata.stderr, 'P', adj_p,...
    'BarFaceColor', data.plotdata.colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2,'PStarThreshold', data.pvalthreshold,...
    'PStarShowNS', false, 'PStarShowGT', false);
set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs); 
xlabel('ROIs');
ylabel('Attention-Passive Difference'); ylim([-.1, .5]);
saveas(gcf, strcat(output.directory, '/attendvspassive-wsrtest.pdf'));

%% ------------------------- INNER VERSUS OUTER ------------------------ %%

% specify inner and outer ROIs to analyse.
ROIs = {'V1I', 'V1O', 'V2I', 'V2O', 'V3I', 'V3O', 'V4I', 'V4O'};

% specify data storage and output directories.
fig.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/univariate/block/group/';
data.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/univariate/block/group/stats/';
output.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/univariate/processed/';

for thisROI = 1:length(ROIs) % for each ROI in turn:
    
    % load in the the ROI-specific univariate betas.
    fig.data = load(strcat(fig.directory, sprintf('%s_block_univariate_betas.mat', ROIs{thisROI})));
    
    % calculate the mean and standard deviation across the three attention
    % conditions.
    fig.plotdata.full{thisROI} = fig.data.betas_pp(:,1:3);
    fig.plotdata.mean(thisROI,:) = mean(fig.plotdata.full{thisROI});
    fig.plotdata.stderr(thisROI,:) = (std(fig.plotdata.full{thisROI})) / sqrt(length(fig.plotdata.full{thisROI}));
    
    % read in the output from the univariate repeated measures anovas and
    % extract the p-values corresponding to the differences between the
    % attention conditions for the inner and outer ROIs.
    fig.pdata = readtable(strcat(data.directory, sprintf('%s_block_original.csv', ROIs{thisROI})));
    fig.plotdata.mainp(thisROI) = table2array(fig.pdata(1,1));
    fig.plotdata.posthoc(thisROI,1,:) = table2array(fig.pdata(1,5));
    fig.plotdata.posthoc(thisROI,2,:) = table2array(fig.pdata(1,6));
    fig.plotdata.posthoc(thisROI,3,:) = table2array(fig.pdata(1,7));
end

% create a p-value matrix containing only significant p-values for plotting.
fig.plotdata.sigind = find(fig.plotdata.mainp > .05);
fig.plotdata.posthoc(fig.plotdata.sigind,:) = NaN;
fig.plotdata.posthoc(fig.plotdata.posthoc>.05)= NaN;

fig.plotdata.pmatrix = nan(numel(fig.plotdata.mean), numel(fig.plotdata.mean));

for thisROI = 1:length(ROIs) % for each ROI in turn:
       
    % combine p-value data across all ROIs into a singular matrix.
    fig.plotdata.pmatrix(thisROI,thisROI+length(ROIs)) = fig.plotdata.posthoc(thisROI,1);
    fig.plotdata.pmatrix(thisROI,thisROI+length(ROIs)*2) = fig.plotdata.posthoc(thisROI,2);
    fig.plotdata.pmatrix(thisROI+length(ROIs), thisROI+length(ROIs)*2) = fig.plotdata.posthoc(thisROI,3);
end

% format the p-value matrix for plotting.
fig.plotdata.pmatrixedit = fig.plotdata.pmatrix'; 
lidx = tril(true(size(fig.plotdata.pmatrix)), -1); 
fig.plotdata.pmatrix(lidx) = fig.plotdata.pmatrixedit(lidx);

% plot the inner and outer ROI data as a bar chart with associated
% significance asterisks and save as .pdf file.
fig.plotdata.colour = [.5 .5 .5];
figure(2);
superbar(fig.plotdata.mean, 'E', fig.plotdata.stderr, 'P', fig.plotdata.pmatrix,...
    'BarFaceColor', fig.plotdata.colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2);
set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs); 
xlabel('ROIs');
ylabel('Univariate Betas'); ylim([-0.5, 2]);
legend('Orient', 'Contrast', 'Shape')
saveas(gcf, strcat(output.directory, '/univariate_block_innerouter.pdf'));

%% ----------------------- UNIVARIATE ATTENTION ------------------------ %% 

% specify ROIs to analyse. 
ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2', 'A1'};

% here, we repeat the same process as above, but for the whole ROI, rather
% than the inner and outer subsections.

for thisROI = 1:length(ROIs) % for each ROI in turn:
    
    fig.data = load(strcat(fig.directory, sprintf('%s_%s_univariate_betas.mat', ROIs{thisROI}, condition)));
    
    if isequal(condition, 'naturalistic')
        fig.plotdata.full{thisROI} = fig.data.betas_pp(:,1:4);
    else
        fig.plotdata.full{thisROI} = fig.data.betas_pp(:,1:3);
    end
    
    fig.plotdata.mean(thisROI,:) = mean(fig.plotdata.full{thisROI});
    fig.plotdata.stderr(thisROI,:) = (std(fig.plotdata.full{thisROI})) / sqrt(length(fig.plotdata.full{thisROI}));
    
    fig.pdata = readtable(strcat(data.directory, sprintf('%s_%s_original.csv', ROIs{thisROI}, condition)));
    fig.plotdata.mainp(thisROI) = table2array(fig.pdata(1,1));
    fig.plotdata.posthoc(thisROI,1,:) = table2array(fig.pdata(1,5));
    fig.plotdata.posthoc(thisROI,2,:) = table2array(fig.pdata(1,6));
    fig.plotdata.posthoc(thisROI,3,:) = table2array(fig.pdata(1,7));
    
    if isequal(condition, 'naturalistic')
        fig.plotdata.posthoc(thisROI,4,:) = table2array(fig.pdata(1,8));
    end
end

fig.pvalthreshold = [0.05, 0.01, 0.001];
[~, ~, ~, adj_p]=fdr_bh(fig.plotdata.mainp, 0.05, 'dep');

fig.plotdata.sigind = find(adj_p > fig.pvalthreshold(1));
fig.plotdata.posthoc(fig.plotdata.sigind,:) = NaN;
fig.plotdata.posthoc(fig.plotdata.posthoc>.05)= NaN;

fig.plotdata.pmatrix = nan(numel(fig.plotdata.mean), numel(fig.plotdata.mean));

for thisROI = 1:length(ROIs)
       
    fig.plotdata.pmatrix(thisROI,thisROI+length(ROIs)) = fig.plotdata.posthoc(thisROI,1);
    fig.plotdata.pmatrix(thisROI,thisROI+length(ROIs)*2) = fig.plotdata.posthoc(thisROI,2);
    fig.plotdata.pmatrix(thisROI+length(ROIs), thisROI+length(ROIs)*2) = fig.plotdata.posthoc(thisROI,3);
end

fig.plotdata.pmatrixedit = fig.plotdata.pmatrix'; lidx = tril(true(size(fig.plotdata.pmatrix)), -1); fig.plotdata.pmatrix(lidx) = fig.plotdata.pmatrixedit(lidx);

fig.plotdata.colour = [.5 .5 .5];
figure(3);
superbar(fig.plotdata.mean, 'E', fig.plotdata.stderr, 'P', fig.plotdata.pmatrix,...
    'BarFaceColor', fig.plotdata.colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2);
set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs); 
xlabel('ROIs');
ylabel('Univariate Betas'); ylim([-0.5, 1.8]);
if isequal(condition, 'colour') || ~(strfind(interaction,'colour'))
    legend('Red-Green', 'Blue-Yellow', 'Luminance')
else
legend('Orient', 'Contrast', 'Shape')
end
saveas(gcf, strcat(output.directory, sprintf('univariate_%s.pdf', condition)));

%% --------------------- UNIVARIATE EVENT-RELATED ---------------------- %% 

% here, we again repeat the same process as above, but for the whole ROI
% and for the event-related univariate betas.

ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2', 'A1'};

for thisROI = 1:length(ROIs)

    fig.data = load(strcat(fig.directory, sprintf('%s_event_univariate_betas.mat', ROIs{thisROI})));
    
    fig.plotdata.full{thisROI} = fig.data.betas_pp(:,1:3);
    fig.plotdata.mean(thisROI,:) = mean(fig.plotdata.full{thisROI});
    fig.plotdata.stderr(thisROI,:) = (std(fig.plotdata.full{thisROI})) / sqrt(length(fig.plotdata.full{thisROI}));
    
    fig.pdata = readtable(strcat(data.directory, sprintf('%s_event_original.csv', ROIs{thisROI})));
    fig.plotdata.mainp(thisROI) = table2array(fig.pdata(1,1));
    fig.plotdata.posthoc(thisROI,1,:) = table2array(fig.pdata(1,5));
    fig.plotdata.posthoc(thisROI,2,:) = table2array(fig.pdata(1,6));
    fig.plotdata.posthoc(thisROI,3,:) = table2array(fig.pdata(1,7));
end

fig.pvalthreshold = [0.05/6, 0.01/6, 0.001/6];
fig.plotdata.sigind = find(fig.plotdata.mainp > fig.pvalthreshold(1));
fig.plotdata.posthoc(fig.plotdata.sigind,:) = NaN;
fig.plotdata.posthoc(fig.plotdata.posthoc>.05)= NaN;

fig.plotdata.pmatrix = nan(numel(fig.plotdata.mean), numel(fig.plotdata.mean));

for thisROI = 1:length(ROIs)
       
    fig.plotdata.pmatrix(thisROI,thisROI+length(ROIs)) = fig.plotdata.posthoc(thisROI,1);
    fig.plotdata.pmatrix(thisROI,thisROI+length(ROIs)*2) = fig.plotdata.posthoc(thisROI,2);
    fig.plotdata.pmatrix(thisROI+length(ROIs), thisROI+length(ROIs)*2) = fig.plotdata.posthoc(thisROI,3);
end

fig.plotdata.pmatrixedit = fig.plotdata.pmatrix'; lidx = tril(true(size(fig.plotdata.pmatrix)), -1); fig.plotdata.pmatrix(lidx) = fig.plotdata.pmatrixedit(lidx);

fig.plotdata.colour = [.5 .5 .5];
figure(4);
superbar(fig.plotdata.mean, 'E', fig.plotdata.stderr, 'P', fig.plotdata.pmatrix,...
    'BarFaceColor', fig.plotdata.colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2);
set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs); 
xlabel('ROIs');
ylabel('Univariate Betas'); ylim([-0.5, 1.5]);
legend('Orient', 'Contrast', 'Shape')
saveas(gcf, strcat(output.directory, '/univariate_event.pdf'));

%% --------------------------------------------------------------------- %%