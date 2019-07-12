%% ------------------- PLOT FMRI MULTIVARIATE RESULTS ------------------ %%

% produces a variety of plots to summarise our multivariate classification
% findings, beyond what are output from the multivariate analysis scripts.

% KWN 10/7/2019

clear; clc;
addpath('/Users/kirstie/Documents/code/fmri_analysis/functions/');

%% ------------------------ THREE-WAY DECODING ------------------------- %%

% specify ROIs to analyse.
ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2'};

% specify data and output directories.
data.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/libsvmclassification/threeway/';
output.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/processed/';

% load in the 3-way classification data. 
original = load(strcat(data.directory, '3wayclassificationdata.mat'));

% extract the ROI-specific mean classification accuracies, test statistic
% and associated significance values.
data.accuracy = original.ROIclassification_accuracy([1,4,5,6,7],:)';
data.stat = original.tstat(:,[1,4,5,6,7]);
data.p = original.tsig(:,[1,4,5,6,7]); data.p(data.p>.05)=NaN;

% calculate the mean and standard error of classification scores. 
data.plotdata.mean = mean(data.accuracy);
data.plotdata.stderr = (std(data.accuracy) / sqrt(length(data.accuracy)));

% plot the data across ROIs on a single bar chart and save as .pdf file.
data.plotdata.colour = [.5 .5 .5];
data.pthreshold = [0.05/length(ROIs) 0.01/length(ROIs) 0.001/length(ROIs)];

figure(1)
superbar(data.plotdata.mean, 'E', data.plotdata.stderr, 'P', data.p,...
    'BarFaceColor', data.plotdata.colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2, 'PStarThreshold', data.pthreshold);
set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs); 
xlabel('ROIs');
ylabel('Classification Accuracy (%)'); ylim([25, 70]);
hline = refline([0 33.33]);
hline.Color = 'r'; hline.LineWidth = 1; hline.LineStyle = '--';

saveas(gcf, strcat(output.directory, '/multivariate_threeway.pdf'));

%% ------------------------- TWO-WAY DECODING -------------------------- %%

% specify data and output directories.
data.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/libsvmclassification/twoway/';
output.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/processed/';

% load in two-way classifcation data. 
original = load(strcat(data.directory, 'twowaydata.mat'));

% extract the ROI-specific mean classification accuracies, test statistic
% and associated significance values.
data.accuracy = cell2mat(original.alldata(:,[1,2,3,10,11,12,13,14,15,16,17,18,19,20,21]));
data.stat = cell2mat(original.tstatdata(:,[1,2,3,10,11,12,13,14,15,16,17,18,19,20,21]));
data.p = cell2mat(original.pvaldata(:,[1,2,3,10,11,12,13,14,15,16,17,18,19,20,21]));
data.p(data.p>0.05/(length(ROIs)*3))=NaN;

% calculate the mean and standard error of classification scores. 
data.plotdata.mean = mean(data.accuracy);
data.plotdata.stderr = (std(data.accuracy) / sqrt(length(data.accuracy)));

data.plotdata.colour = [.5 .5 .5];
data.pthreshold = [0.05/(length(ROIs)*3) 0.01/(length(ROIs)*3) 0.001/(length(ROIs)*3)];

% plot the data across ROIs on individual subplots and save as .pdf file.
ROIstart = 1;
for thisROI = 1:length(ROIs) 
    subplot(2,3,thisROI)
    
    superbar(data.plotdata.mean(:,[ROIstart:ROIstart+2]), 'E', data.plotdata.stderr(:,[ROIstart:ROIstart+2]), 'P', data.p(:,[ROIstart:ROIstart+2]),...
        'BarFaceColor', data.plotdata.colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
        'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
        'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2, 'PStarThreshold', data.pthreshold, 'PStarShowGT', false);
    set(gca, 'xtick', 1:3, 'xticklabel', {'OC', 'OS', 'CS'});
    xlabel('Pairwise Comparisons');
    ylabel('Classification Accuracy (%)'); ylim([25, 100]);
    hline = refline([0 50]);
    hline.Color = 'r'; hline.LineWidth = 1; hline.LineStyle = '--';
    
    ROIstart = ROIstart + 3;
end

saveas(gcf, strcat(output.directory, '/multivariate_twoway.pdf'));

%% --------------------------------------------------------------------- %%