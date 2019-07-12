%%  PLOT MULTIDIMENSIONAL SCALING CORRELATIONS & SPATIAL REPRESENTATION  %%

% reads in condition-specific univariate beta values sorted by percentage
% variance explained, performs bootstrapping to ensure all conditions have
% the same number of data points and correlates these univariate
% beta values and plots the similarity between conditions in space.

% 12/7/19 KWN

clear; clc;
addpath(genpath('/scratch/sg3/P1323/code/fmri_analysis/functions/'));

%% --------------------- SPECIFY ANALYSIS PARAMETERS ------------------- %%

cmap = fireice; % load in colourmap.

% specify ROIs to analyse,
ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2', 'IPS0', 'A1'};

% specify conditions to analyse.
condlabels = {'passive', 'vertical', 'horizontal', 'diagonal', 'red',...
    'green', 'blue', 'circular', 'square', 'triangular', 'face'};

% initialise participant counter variable.
ppcount = 0;

%% --------------------------- BOOTSTRAP DATA -------------------------- %%

for thisROI = 1:length(ROIs) % for each ROI in turn:
    ROI = ROIs{thisROI};
    
    % extract the participant univariate beta files.
    files = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3x3feature/betasxvarexp/%s/100/R*', ROI));
    filenames = {files.name};
    
    reshapecondppdata=[];
    
    for thisFile = 1:length(filenames) % for each participant in turn:
        
        % load in the participant-specific beta values.
        load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3x3feature/betasxvarexp/%s/100/%s', ROI, filenames{thisFile}));
        
        % increment participant counter variable.
        ppcount = ppcount +1;
        
        reshapeconddata = [];
        
        for thiscond = 1:size(condlabels) % for each condition in turn:
            
            % extract the condition specific data and remove any nan
            % values.
            conddata = squeeze(final_betas(:,thiscond,:));
            conddata(all(isnan(conddata),2),:) = [];
            
            % if we do not have data for the full number of condition
            % repetitions.
            if size(conddata,1) ~= 8
                
                % calculate the number of bootstrapped samples required.
                origsize = size(conddata,1);
                bstrpsize = 8 - origsize;
                datavalues = 1:origsize;
                sampnumber = origsize;
                
                count = origsize + 1;
                
                % create bootstrapped samples from the condition-specific
                % data set.
                for bstrpsamp = 1:bstrpsize
                    datasamp = randsample(datavalues,sampnumber,'true');
                    dataselect(bstrpsamp,:) = mean(conddata(datasamp,:),1);
                    conddata(count,:) = dataselect(bstrpsamp,:);
                    count = count+1;
                end
            end
            
            % add this data to overall storage arrays.
            conddata = mean(conddata,1);
            reshapeconddata = [reshapeconddata; conddata];
        end
        processdata(thisROI,ppcount,:,:) = reshapeconddata;
        reshapecondppdata = [reshapecondppdata,squeeze(processdata(thisROI,thisFile,:,:))];
    end
    
    % combine data across all ROIs in a singular matrix.
    superROIdata(thisROI,:,:) = reshapeconddata;
    ppcount = 0; % refresh participant counter variable.
end

%%  COMPUTE AND PLOT CORRELATIONS BETWEEN CONDITION-SPECIFIC BETA VALUES %%

% calculate the correlations between conditions for each ROI.
for thisROI = 1:size(ROIs)
    [corrmat(thisROI,:,:), pval(thisROI,:,:)] = ...
        corr(squeeze(superROIdata(thisROI,:,:))', squeeze(superROIdata(thisROI,:,:))');
end

for thisROI = 1:size(ROIs)
    
    % subtract the condition-specific correlations from the control ROI
    % (A1) correlation values for normalisation.
    fig = figure(thisROI);
    data = squeeze(corrmat(thisROI,:,:))-squeeze(corrmat(7,:,:));
    
    % plot the data as a correlation matrix and save as .png file.
    imagesc(data, [-1 1]);
    title(ROIs{thisROI});
    set(gca, 'YTick', 1:11); set(gca, 'YTickLabel', condlabels); set(gca, 'YTickLabelRotation', 45);
    set(gca, 'XTick', 1:11); set(gca, 'XTickLabel', condlabels); set(gca, 'XTickLabelRotation', 45);
    colormap(cmap); colorbar;
    saveas(fig, sprintf('%s_rsa_matrix.png', ROIs{thisROI}));
end

%% --------------- PLOT SPATIAL SIMILARITY OF CONDITIONS --------------- %%

for thisROI = 1:size(ROIs) % for each ROI in turn:
    
    % perform multidimensional scaling to calculate distances between
    % conditions.
    dm = squeeze(corrmat(thisROI,:,:));
    dm = 1-(dm+1)/2;
    Y = cmdscale(dm,2);
    
    % plot the position of each category on a scatterplot for each ROI and
    % save as .png file.
    fig = figure(thisROI+7);
    plot(Y(:,1),Y(:,2),'.');
    text(Y(:,1)+.01, Y(:,2), condlabels);
    set(gca, 'YLim', [-.4, .7]);
    set(gca, 'XLim', [-.4, .7]);
    title(ROIs{thisROI});
    saveas(fig, sprintf('%s_rsa_mds.png', ROIs{thisROI}));
end

%% --------------------------------------------------------------------- %%