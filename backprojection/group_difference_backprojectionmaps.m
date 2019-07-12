%% ------------------ CALCULATE GROUP DIFFERENCE MAPS ------------------ %%

% loads in the group data (an interpolated grid of values for each
% participant), with no averaging or normalisation, and calculates the
% difference between these two conditions, then plots the mean group 
% difference map, with z-score thresholding. 

% 09/04/2018 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

clear; clc; close all; % general houskeeping. 
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions'); % add function directory to path.

outputROIs = {'V1', 'V2', 'V3', 'V4'}; % specify the ROIs we wish to analyse. 

% specify the combination of features we wish to calculate differences
% between. 
features{1} = {'orientation', 'contrast'};
features{2} = {'orientation', 'shape'};
features{3} = {'contrast', 'shape'};

cmap = fireice; % load the hot-cold colourmap. 
groupthresh = [1.7 -1.7]; % specify the upper and lower z-score thresholding desired. 

%% ------------------- LOAD IN RFP FOR OVERLAY PLOTS ------------------- %%

[distances] = RFP_generator(0,0.5,0.5); distances = double(distances); % create a reference RFP.
distances = imrotate(distances,90,'crop'); % rotate by 90 for correct display on the surface.
        
% replace any values outside of a range of values (for the line
% thickness of the RFP) to NaN, and the RFP line outline to 1.
distances(distances < 181 | distances > 191) = NaN;
distances(~isnan(distances)) = 1;
RFPind = find(distances==1); % store the positions of RFP line data.

% create a grid of x and y values matching the size of the interpolated
% difference data plot for later overlay. 
[distxsamples, distysamples] = deal(linspace(-6, 6, 500));
[distXQ, distYQ] = meshgrid(distxsamples, distysamples);

%% -------------------- SPECIFY INTERPOLATION GRID --------------------- %%

% specify the grid of radius and eccentricity the data was originally
% interpolated over, for later plotting. 
xsamples = logspace(-2, log10(6), 500);
ysamples = linspace(0, pi*2, 500);

[XQ, YQ] = meshgrid(xsamples, ysamples); % grid across the full desired sampling range.
[XQ, YQ] = pol2cart(YQ, XQ); % convert from polar to cartesian coordinates.
        
%% ----------------- CALCULATE & PLOT DIFFERENCE MAPS ------------------ %%

for x = 1:length(outputROIs) % for each ROI in turn:
    
    outputROI = outputROIs{x}; % extract the ROI-specific output file name. 
    
    % specify the ROI-specific output data directory. 
    outputdir = sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/difference/%s/', outputROI);
    
    for i = 1:length(features) % for each combination of features in turn: 
        
        % load in the two condition-specific group interpolated data matrices. 
        data1 = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/original/group/%s/group_%s%s_interpolateddata_polar.mat', outputROI, features{i}{1}, outputROI));
        data2 = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/original/group/%s/group_%s%s_interpolateddata_polar.mat', outputROI, features{i}{2}, outputROI));
       
        % calculate the difference between these two condition matrices. 
        diff = data1.featuremap - data2.featuremap;
        diff = nanmean(diff,3)./nanstd(diff,[],3); % zscore the data. 
        
        % threshold, so that we only display values beyond the specified
        % thresholding range. 
        diff(diff<groupthresh(1) & diff>groupthresh(2)) = 0;
        
        figure(1); surf(XQ,YQ,diff); shading interp  % plot the data as a 3D surface.
        colorbar;colormap(cmap); caxis([-6 6]); view(90, -90);
        
        % add the RFP overlay to the plot. 
        hold on; plot3(distXQ, distYQ,  distances-10, 'Color', [1 1 1], 'LineWidth',2);

        % save this plot with an ROI and difference-feature specific name. 
        saveas(gca, strcat(outputdir, sprintf('group_%s%s-%smanualdifference-zscorethresh%.2f.png', features{i}{1}, features{i}{2}, outputROI, groupthresh(1))));

        close all; % close any open figure windows. 
    end % continue for the next feature combination.
end % continue for the next ROI. 

%% --------------------------------------------------------------------- %%