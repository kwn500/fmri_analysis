%%  ANALYSE DIFFERENCE IN BOLD MODULATION BETWEEN INNER AND OUTER ROIS   %%

% loads in the univariate betas for each participant for the inner and
% outer ROI counterparts, and calculates a mean for each across the three
% attention conditions (orientation, contrast and shape). then performs a
% paired-sample t-test analysing difference in signal modulation between
% the ROI locations. 

% KWN 9/7/2019

clear; clc;

%%

ROIs = {'V1', 'V2', 'V3', 'V4'}; % specify ROIs to analyse. 

for i = 1:length(ROIs) % for each ROI in turn:

    % load in the ROI-specific inner and outer univariate betas across all
    % participants.
    inner = load(sprintf('/Users/kirstie/Documents/analysis/originalRFP/univariate/raw/inner_outer/%sI_block_univariate_betas.mat', ROIs{i}));
    outer = load(sprintf('/Users/kirstie/Documents/analysis/originalRFP/univariate/raw/inner_outer/%sO_block_univariate_betas.mat', ROIs{i}));
    
    % calculate the mean for each participants across the three attention
    % conditions. 
    innermean = mean(inner.betas_pp(:,1:3),2);
    outermean = mean(outer.betas_pp(:,1:3),2);
    
    % perform a paired-sample t-test on the inner versus outer data and 
    % store the ROI-specific values.  
    [data.H(i),data.P(i),data.CI{i},data.STATS{i}] = ttest(innermean, outermean);
end

% save the statistical output to a .mat file. 
save('/Users/kirstie/Documents/analysis/originalRFP/univariate/inner_outer_ttest.mat', 'data');

%% --------------------------------------------------------------------- %%