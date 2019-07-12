%% ------ ANALYSE PARAMETER MAP DIFFERENCES BY FEATURE THRESHOLDS ------ %%

% Loads in the participants' standardised interpolated support vector
% parameter maps, for orientation, contrast and shape, along with each
% participants feature threshold (from previous psychophysics testing).
% Then splits the data into lower and higher performance (currently
% splitting the data into two halves arbitrarily and average the parameter
% maps for these two data halves to look for differences in the pattern of
% activation across these two data groups.

% 03/09/2018 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

clear; clc; % general housekeeping.

% add function directory to path.
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions');

% specify the threshold data storage directory.
group_directory =('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/threshold_correlation/');

% specify the output directory.
outputdir = ('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/performance/');

% specify the abbreviated name of the ROI we wish to analyse.
outputROIs = {'V1', 'V2', 'V3', 'V4'};

groupthresh = [1.7 -1.7]; % max & min boundaries for group figure thresholding.

%% -------------------- GENERATE RFP FOR OVERLAY ----------------------- %%

[distances] = RFP_generator(0,0.5,0.5); distances = double(distances); % create a reference RFP.
distances = imrotate(distances,90,'crop'); % rotate by 90 for correct display on the surface.

% replace any values outside of a range of values (for the line
% thickness of the RFP) to NaN, and the RFP line outline to 1.
distances(distances < 181 | distances > 191) = NaN;
distances(~isnan(distances)) = 1;
RFPind = find(distances==1); % store the positions of RFP line data.
[distxsamples, distysamples] = deal(linspace(-6, 6, 500));
[distXQ, distYQ] = meshgrid(distxsamples, distysamples);

%% ----------------- CREATE HIGHER/LOWER PARAMETER MAPS ---------------- %%

% load in the participants' feature thresholds (orientation, contrast and
% shape) from previous psychophysics testing.
load('/scratch/groups/Projects/P1323/original/fMRI/general_output/psychophysics_thresholds_original.mat');
data(11,:) = []; % remove participant data with variance explained issue. 

% store the corresponding columns of data accordingly.
pps = data(:,1);
othresh = [(1:length(pps))', data(:,2)]; othresh = sortrows(othresh,2);
cthresh = [(1:length(pps))', data(:,3)]; cthresh = sortrows(cthresh,2);
sthresh = [(1:length(pps))', data(:,4)]; sthresh = sortrows(sthresh,2);

% extract the higher and lower threshold data for each feature.
lower_othresh = othresh(1:6,:); higher_othresh = othresh(7:11,:);
lower_cthresh = cthresh(1:6,:); higher_cthresh = cthresh(7:11,:);
lower_sthresh = sthresh(1:6,:); higher_sthresh = sthresh(7:11,:);

for i =1:length(outputROIs)
    outputROI = outputROIs{i};
    
    % load in the participants' interpolated parameter map data for the three
    % visual features respectively.
    ointerp = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/original/group/%s/group_orientation%s_interpolateddata_polar.mat', outputROI, outputROI));
    cinterp = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/original/group/%s/group_contrast%s_interpolateddata_polar.mat', outputROI, outputROI));
    sinterp = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/backprojection/interpolation/original/group/%s/group_shape%s_interpolateddata_polar.mat',outputROI, outputROI));
    
    % extract the higher and lower interpolated map data for each feature.
    lower_ointerp = ointerp.featuremap(:,:,lower_othresh(:,1)); higher_ointerp = ointerp.featuremap(:,:,higher_othresh(:,1));
    lower_cinterp = cinterp.featuremap(:,:,lower_cthresh(:,1)); higher_cinterp = cinterp.featuremap(:,:,higher_cthresh(:,1));
    lower_sinterp = sinterp.featuremap(:,:,lower_sthresh(:,1)); higher_sinterp = sinterp.featuremap(:,:,higher_sthresh(:,1));
    
    % normalise (zscore) each set of feature-specific higher and lower data
    % and then threshold so we only present data above and below the
    % specified values above.
    normalisedlo = nanmean(lower_ointerp,3)./nanstd(lower_ointerp,[],3);
    normalisedlo(normalisedlo < groupthresh(1) & normalisedlo > groupthresh(2)) = 0; % threshold.
    
    normalisedlc = nanmean(lower_cinterp,3)./nanstd(lower_cinterp,[],3);
    normalisedlc(normalisedlc < groupthresh(1) & normalisedlc > groupthresh(2)) = 0; % threshold.
    
    normalisedls = nanmean(lower_sinterp,3)./nanstd(lower_sinterp,[],3);
    normalisedls(normalisedls < groupthresh(1) & normalisedls > groupthresh(2)) = 0; % threshold.
    
    normalisedho = nanmean(higher_ointerp,3)./nanstd(higher_ointerp,[],3);
    normalisedho(normalisedho < groupthresh(1) & normalisedho > groupthresh(2)) = 0; % threshold.
    
    normalisedhc = nanmean(higher_cinterp,3)./nanstd(higher_cinterp,[],3);
    normalisedhc(normalisedhc < groupthresh(1) & normalisedhc > groupthresh(2)) = 0; % threshold.
    
    normalisedhs = nanmean(higher_sinterp,3)./nanstd(higher_sinterp,[],3);
    normalisedhs(normalisedhs < groupthresh(1) & normalisedhs > groupthresh(2)) = 0; % threshold.
    
    % add this normalised data to a singular storage cell array for easier
    % figure generation.
    datacombined{1} = normalisedlo; datacombined{2} = normalisedlc; datacombined{3} = normalisedls;
    datacombined{4} = normalisedho; datacombined{5} = normalisedhc; datacombined{6} = normalisedhs;
    
    % specify the performance-specific figure titles.
    titles{1} = 'High Performance Orientation'; titles{2} = 'High Performance Contrast'; titles{3} = 'High Performance Shape';
    titles{4} = 'Low Performance Orientation';  titles{5} = 'Low Performance Contrast';  titles{6} = 'Low Performance Shape';
    
    % create a meshgrid of eccentricity and polar angle values matching the
    % interpolation script, in polar co-ordinates, converted to cartesian.
    
    xsamples = logspace(-2,log10(6),500);
    ysamples = linspace(0, pi*2, 500);
    [XQ, YQ] = meshgrid(xsamples, ysamples);
    [XQ, YQ] = pol2cart(YQ, XQ);
    
    cmap = fireice; % load in the hot/cold colour map.
    
    % for each of the 6 data matrices (higher and lower x three visual
    % features), plot the data as a 3D surface with an RFP overlain.
    fig= figure('Position', get(0, 'Screensize'));  % open a full-screen figure for better visualisation.
    [ha, pos] = tight_subplot(2,3,[0.05, 0.03],0.05,0.03); % create an arrangement of 6 tightly-tiled subplots.
    
    for i = 1:length(datacombined) % for each analysis condition in turn:
        
        axes(ha(i)); % within an analysis-specific subplot:
        plotdata = datacombined{i}; % specify the analysis-specific data.
        
        surf(XQ,YQ,plotdata);shading interp; % plot this data as a 3D surface mesh. 
        view(90,-90); colorbar; colormap(cmap); caxis([-4 4]);title(titles{i});
        
        % add the RFP overlay. 
        hold on; plot3(distXQ, distYQ,  distances-10, 'Color', [1 1 1], 'LineWidth',2);
    end % repeat for the next analysis condition. 

    % save the figure as a .png image file with an ROI-specific name.
    fig.Color = 'white'; saveas(gca, [outputdir, sprintf('group_performance_%s_thresh%.2f.png', outputROI, groupthresh(1))]);
    
    close all; % close any open figures. 
end % repeat for the next ROI. 

%% --------------------------------------------------------------------- %%