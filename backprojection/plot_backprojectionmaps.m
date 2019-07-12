%% --- STANDARDISE & AVERAGE SUPPORT VECTOR MAPS ACROSS PARTICIPANTS --- %%

% For each participant and specific-ROI, loads up the pRF session and
% extracts eccentricity and polar angle information for each voxel within
% that ROI. Also loads in the previously-created (vs passive) support
% vector parameter maps. We interpolate the eccentricity and polar angle
% data, to form a standard grid for each participant, convert this data to
% cartesian co-ordinates and plot the interpolated parameter map data, to
% demonstrate the arrangement of feature-specific voxel activation
% associated with a voxels particular eccentricity and polar angle.

% added compatibilty for colour experiment (24/09/2018).

% 29/08/18 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

clear all; close all; clc; % general housekeeping.
ynicInit spm8

addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions'); % add function directory to path.

exp_condition = 'colour';
condition = 'feature';

% specify participant numbers and their corresponding pRF directories
% (these participants vary as a function of experiment we are analysing).
baseprf_dir = '/groups/labs/wadelab/data/pRF/Himmel/TemporalContrastpRF/';

if isequal(exp_condition, 'original')
    participants{1} = 'R2268'; prf_dir{1} = [baseprf_dir, '/pRF_data_MMH_Collected/R2268'];
    participants{2} = 'R2548'; prf_dir{2} = [baseprf_dir, '/pRF/R2548/CombinedSessions'];
    participants{3} = 'R2590'; prf_dir{3} = [baseprf_dir, '/pRF_data_MMH_Collected/R2590'];
    participants{4} = 'R2904'; prf_dir{4} = [baseprf_dir, '/pRF_data_MMH_Collected/R2904'];
    participants{5} = 'R3111'; prf_dir{5} = [baseprf_dir, '/pRF/R3111/CombinedSessions'];
    participants{6} = 'R3455'; prf_dir{6} = [baseprf_dir, '/pRF/R3455/CombinedSessions'];
    participants{7} = 'R3517'; prf_dir{7} = [baseprf_dir, '/pRF_data_MMH_Collected/R3517'];
    participants{8} = 'R3773'; prf_dir{8} = [baseprf_dir, '/pRF_data_MMH_Collected/R3773'];
    participants{9} = 'R3932'; prf_dir{9} = [baseprf_dir, '/pRF_data_MMH_Collected/R3932'];
    participants{10} = 'R4065'; prf_dir{10} = [baseprf_dir, '/pRF_data_MMH_Collected/R4065'];
    %     participants{11} = 'R4127'; prf_dir{11} = [baseprf_dir, '/pRF_data_MMH_Collected/R4127'];
    participants{11} = 'R4496'; prf_dir{11} = [baseprf_dir, '/pRF_data_MMH_Collected/R4496'];
elseif isequal(exp_condition, 'colour')
    participants{1} = 'R2268'; prf_dir{1} = [baseprf_dir, '/pRF_data_MMH_Collected/R2268'];
    participants{2} = 'R2548'; prf_dir{2} = [baseprf_dir, '/pRF/R2548/CombinedSessions'];
    participants{3} = 'R2590'; prf_dir{3} = [baseprf_dir, '/pRF_data_MMH_Collected/R2590'];
    participants{4} = 'R2904'; prf_dir{4} = [baseprf_dir, '/pRF_data_MMH_Collected/R2904'];
    participants{5} = 'R3111'; prf_dir{5} = [baseprf_dir, '/pRF/R3111/CombinedSessions'];
    participants{6} = 'R3517'; prf_dir{6} = [baseprf_dir, '/pRF_data_MMH_Collected/R3517'];
    participants{7} = 'R3773'; prf_dir{7} = [baseprf_dir, '/pRF_data_MMH_Collected/R3773'];
    participants{8} = 'R4059'; prf_dir{8} = [baseprf_dir, '/pRF_data_MMH_Collected/R4059'];
    participants{9} = 'R4065'; prf_dir{9} = [baseprf_dir, '/pRF_data_MMH_Collected/R4065'];
    participants{10} = 'R4244'; prf_dir{10} = [baseprf_dir, '/pRF_data_MMH_Collected/R4244_redo'];
    participants{11} = 'R4831'; prf_dir{11} = '/scratch/groups/Projects/P1323/colour/fMRI/R4831/retinotopy/mrvista';
    participants{12} = 'R4928'; prf_dir{12} = '/scratch/groups/Projects/P1323/colour/fMRI/R4928/retinotopy/mrvista';
elseif isequal(exp_condition, 'naturalistic')
    participants{1} = 'R2548'; prf_dir{1} = [baseprf_dir, '/pRF/R2548/CombinedSessions'];
    participants{2} = 'R2904'; prf_dir{2} = [baseprf_dir, '/pRF_data_MMH_Collected/R2904'];
    participants{3} = 'R3111'; prf_dir{3} = [baseprf_dir, '/pRF/R3111/CombinedSessions'];
    participants{4} = 'R3517'; prf_dir{4} = [baseprf_dir, '/pRF_data_MMH_Collected/R3517'];
    participants{5} = 'R4059'; prf_dir{5} = [baseprf_dir, '/pRF_data_MMH_Collected/R4059'];
    participants{6} = 'R4127'; prf_dir{6} = [baseprf_dir, '/pRF_data_MMH_Collected/R4127'];
    participants{7} = 'R4244'; prf_dir{7} = [baseprf_dir, '/pRF_data_MMH_Collected/R4244_redo'];
    participants{8} = 'R4831'; prf_dir{8} = '/scratch/groups/Projects/P1323/colour/fMRI/R4831/retinotopy/mrvista';
    participants{9} = 'R4833';  prf_dir{9} = '/scratch/groups/Projects/P1361/fMRI/R4833/retinotopy/';
    participants{10} = 'R4890'; prf_dir{10} = '/scratch/groups/Projects/P1361/fMRI/R4890/retinotopy/mrvista';
    participants{11} = 'R4928'; prf_dir{11} = '/scratch/groups/Projects/P1323/colour/fMRI/R4928/retinotopy/mrvista';
    participants{12} = 'R5006'; prf_dir{12} = '/scratch/groups/Projects/P1361/fMRI/R5006/retinotopy';
end

% specify root directory.
if isequal(exp_condition, 'naturalistic')
    basedir = sprintf('/scratch/groups/Projects/P1361/fMRI/');
else
    basedir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/', exp_condition);
end

% specify ROI for analysis, and it's abbreviated format for file outputs.
ROIanalysis = '03_Combined_V3_KWN'; outputROI = 'V3';

feature = 'Orientation'; % feature to analyse.
analysis_type = 'supportvector'; % data type to analyse (betas/supportvector).
parametermapname = 'Orientation'; % parameter map to load (mrvista).
parametermaptype  = 'original'; % parameter map type (original/difference).
coordtoggle = 'polar'; % coordinate system for plotting (cartesian/polar).
overlaytoggle = 1; % toggle for RFP overlay (1--yes, 0--no).
annulustoggle = 1; % toggle for annulus analysis (1--yes, 0--no).
indthresh = [1.7 -1.7]; % max & min boundaries for ind figure thresholding (display above these values).
groupthresh = [1.7 -1.7]; % max & min boundaries for group figure thresholding.

%% --------- LOAD REFERENCE RADIAL FREQUENCY PATTERN OUTLINE ----------- %%

if overlaytoggle == 1 % if an RFP overlay is requested:
    [distances] = RFP_generator(0,0.5,0.5); distances = double(distances); % create a reference RFP.
    distances = imrotate(distances,90,'crop'); % rotate by 90 for correct display on the surface.
    
    % replace any values outside of a range of values (for the line
    % thickness of the RFP) to NaN, and the RFP line outline to 1.
    distances(distances < 181 | distances > 191) = NaN;
    distances(~isnan(distances)) = 1;
    RFPind = find(distances==1); % store the positions of RFP line data.
end

%% ------ PRODUCE PARTICIPANT-SPECIFIC STANDARDISED PARAMETER MAPS ----- %%

for participantnumber = 1:length(participants) % for each participant in turn:
    
    participant = participants{participantnumber}; % extract the participant-specific number.
    
    % if the participant does not already have extracted prf data for this ROI:
    if sum(exist([basedir, sprintf('/%s/anatomy/ROIs/pRFdata/%s%s_pRF_data.mat',...
            participant, participant, outputROI)], 'file')) == 0
        extractpRFdata(prf_dir, exp_condition, participantnumber, participant,...
            ROIanalysis, outputROI); % extract and save the prf data for all voxels in this ROI.
    end
    
    if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
        load([basedir, sprintf('/%s/anatomy/MT/ROIs/pRFdata/%s%s_pRF_data.mat',...
            participant, participant, outputROI)]); % load ROI-specific data.
    else
        load([basedir, sprintf('/%s/anatomy/ROIs/pRFdata/%s%s_pRF_data.mat',...
            participant, participant, outputROI)]); % load ROI-specific data.
    end
    
    % load in the ROI- & feature-specific parameter map.
    if isequal(parametermaptype,'original')
        if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
            load([basedir, sprintf('/%s/mrvista/%s/MT/Gray/GLMs/%s/%s%s-%s.mat',...
                participant, condition, analysis_type, participant, parametermapname, outputROI)]);
        else
            load([basedir, sprintf('/%s/mrvista/%s/Gray/GLMs/%s/%s%s-%s.mat',...
                participant, condition, analysis_type, participant, parametermapname, outputROI)]);
        end
    elseif isequal(parametermaptype,'difference')
        if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
            load([basedir, sprintf('%s/mrvista/%s/MT/Gray/GLMs/%s/difference/%s%s-%s.mat',...
                participant, condition, analysis_type, participant, parametermapname, outputROI)]);
        else
            load([basedir, sprintf('%s/mrvista/%s/Gray/GLMs/%s/difference/%s%s-%s.mat',...
                participant, condition, analysis_type, participant, parametermapname, outputROI)]);
        end
    end
    
    if isequal(exp_condition, 'naturalistic') && isequal(participant, 'R3517')
        load([basedir, sprintf('/%s/mrvista/session1-allscans/Gray/GLMs/supportvector/%s%s-%s.mat', participant, participant, parametermapname, outputROI)]);
    else
        load([basedir, sprintf('/%s/mrvista/Gray/GLMs/supportvector/%s%s-%s.mat', participant, participant, parametermapname, outputROI)]);
    end
    
    z = map{1}; z(z == 0) = []; % extract the ROI-specific parameter map values.
    
    if isequal(coordtoggle, 'cartesian')
        
        x = ROIx; y = ROIy; % extract the ROI-specific x and y coordinate values.
        
        [xsamples, ysamples] = deal(linspace(-6, 6, 500)); % specify range of interpolation values.
        [XQ, YQ] = meshgrid(xsamples, ysamples); % grid across the full desired sampling range.
        
    elseif isequal(coordtoggle, 'polar')
        
        x = ROIecc; y = ROIpol; % extract the ROI-specific eccentricity and polar angle values.
        
        xsamples = logspace(-2, log10(6), 500); % specify range of interpolation values.
        ysamples = linspace(0, pi*2, 500);
        
        [XQ, YQ] = meshgrid(xsamples, ysamples); % grid across the full desired sampling range.
        [x, y, z] = pol2cart(y,x,z); [XQ, YQ] = pol2cart(YQ, XQ); % convert from polar to cartesian coordinates.
    end
    
    % create 500x500 cartesian grid for plotting of the RFP.
    [distxsamples, distysamples] = deal(linspace(-6, 6, 500));
    [distXQ, distYQ] = meshgrid(distxsamples, distysamples);
    
    if annulustoggle == 1 % if we specified an annulus analysis:
        % create an annulus meshgrid, with a ring of data with an average
        % radius 2 degrees from the centre in polar coordinates.
        annulusxsamples = linspace(1.5, 2.5, 100);
        annulusysamples = linspace(0, pi*2, 360);
        [annulusXQ, annulusYQ] = meshgrid(annulusxsamples, annulusysamples);
        [annulusXQ, annulusYQ] = pol2cart(annulusYQ, annulusXQ);
        [annulusx, annulusy, annulusz] = pol2cart(y,x,z); % convert to cartesian coordinates.
        annulusinterpdata = griddata(annulusx, annulusy, annulusz, annulusXQ, annulusYQ, 'natural');
        
        % normalise this data for plotting.
        annulusmean = nanmean(nanmean(annulusinterpdata));
        annulusstd = nanmean(nanstd(annulusinterpdata));
        plotannulus = (annulusinterpdata-annulusmean)./annulusstd;
    end
    
    interpdata = griddata(x,y,z,XQ,YQ,'natural'); % interpolate this data across the specified grid.
    
    cmap = fireice; % load hot/cold colour map.
    
    plotmean = nanmean(nanmean(interpdata)); % calculate overall mean.
    plotstd = nanmean(nanstd(interpdata)); % calculate overall std.
    plotdata = (interpdata-plotmean)./plotstd; % z-score every element in the matrix.
    
    % set any value in the thresholding range (between the minimum and
    % maximum) to 0- we display values outside of this range.
    plotdata(plotdata < indthresh(1) & plotdata > indthresh(2)) = 0;
    
    % plot the normalised interpolated data as a 3D surface mesh.
    figure(1); surf(XQ,YQ,plotdata); shading interp
    colorbar; colormap(cmap); caxis([-6 6]); view(90, -90);
    
    if overlaytoggle == 1 % if we specified an RFP overlay, add this to the plot.
        hold on; plot3(distXQ, distYQ,  distances-10, 'Color', [1 1 1], 'LineWidth',2);
    end
    
    % create a participant-specific output directory.
    ppdir = ([basedir, sprintf('/vista_output/multivariate/%s/backprojection/interpolation/%s/%s/',...
        condition, parametermaptype, participant)]); [~,~] = mkdir(ppdir);
    
    featuremap(:,:,participantnumber) = interpdata; % add participant data to group interpolated storage matrix.
    plotfeaturemap(:,:,participantnumber) = plotdata; % add participant plotting data to group storage matrix.
    
    % save plot as 3D .fig and 2D .png image.
    savefig(figure(1), strcat(ppdir,sprintf('%s%s%s_surface_%s_thresh%.2f.fig',...
        participant, lower(feature), outputROI, coordtoggle, indthresh(1))));
    saveas(gca, strcat(ppdir, sprintf('%s%s%s_imageplane_%s_thresh%.2f.png',...
        participant, lower(feature), outputROI, coordtoggle, indthresh(1))));
    
    if annulustoggle == 1 % if we specified an annulus overlay:
        
        hold on; % add this to the plot, with a mean radius of 2 degrees.
        annulusmean = nanmean(plotannulus,2)';
        annulusmap(:,:,participantnumber) = plotannulus;
        
        theta = linspace(0,360,360); theta = deg2rad(theta);
        annulusmean = annulusmean+2; [x,y] = pol2cart(theta,annulusmean);
        
        plot3(x,y,ones(1,360)*-3,'Color', [0.5 0.5 0.5],'LineWidth', 3);
        
        % save the figures as before with this annulus overlay.
        savefig(figure(1), strcat(ppdir,sprintf('%s%s%s_surface_%s_thresh%.2f_annulus.fig',...
            participant, lower(feature), outputROI, coordtoggle, indthresh(1))));
        saveas(gca, strcat(ppdir, sprintf('%s%s%s_imageplane_%s_thresh%.2f_annulus.png',...
            participant, lower(feature), outputROI, coordtoggle, indthresh(1))));
    end
    
    close all % close any remaining open windows.
end % repeat the process for the next participant.

%% --------- PRODUCE AVERAGE GROUP STANDARDISED PARAMETER MAPS --------- %%

% create a group output storage directory.
groupdir = [basedir, sprintf('/vista_output/multivariate/%s/backprojection/interpolation/%s/group/',...
    condition, parametermaptype)]; [~,~] = mkdir(groupdir);

save(strcat(groupdir, sprintf('group_%s%s_interpolateddata_%s.mat', lower(feature),...
    outputROI, coordtoggle)),... 'featuremap'); % save group interpolated data.
    
save(strcat(groupdir, sprintf('group_%s%s_plotdata_%s.mat', lower(feature),...
    outputROI, coordtoggle)), 'plotfeaturemap'); % save group plot data.

if annulustoggle == 1 % if we specified an annulus overlay, load in this annulus data.
    save(strcat(groupdir, sprintf('group_%s%s_annulusdata.mat', lower(feature), outputROI)), 'annulusmap');
end

% z-score across every participant in the interpolated data matrix.
normaliseddata = nanmean(featuremap, 3)./nanstd(featuremap,[], 3);
normaliseddata(normaliseddata < groupthresh(1) & normaliseddata > groupthresh(2)) = 0; % threshold.

% plot the normalised interpolated data as a 3D surface mesh.
figure(1); surf(XQ,YQ,normaliseddata); shading interp
colorbar;colormap(cmap); caxis([-6 6]); view(90, -90);

if overlaytoggle == 1 % if we specified an RFP overlay, add this to the plot.
    hold on; plot3(distXQ, distYQ,  distances-10, 'Color', [1 1 1], 'LineWidth',2);
end

% save plot as 3D .fig and 2D .png image.
savefig(figure(1), strcat(groupdir,sprintf('group%s%s_surface_%s_thresh%.2f.fig',...
    lower(feature), outputROI, coordtoggle, groupthresh(1))));
saveas(gca, strcat(groupdir, sprintf('group%s%s_imageplane_%s_thresh%.2f.png',...
    lower(feature), outputROI, coordtoggle, groupthresh(1))));

if annulustoggle == 1 % if we specified an annulus overlay:
    
    hold on; % add this to the plot, with a mean radius of 2 degrees.
    annulusmean = nanmean(annulusmap,3); annulusmean = nanmean(annulusmean,2)';
    annulusmean = annulusmean+2; [x,y] = pol2cart(theta,annulusmean);
    
    plot3(x,y,ones(1,360)*-3,'Color', [0.5 0.5 0.5],'LineWidth', 3);
    
    % save the figures as before with this annulus overlay.
    savefig(figure(1), strcat(groupdir,sprintf('group%s%s_surface_%s_thresh%.2f_annulus.fig',...
        lower(feature), outputROI, coordtoggle, groupthresh(1))));
    saveas(gca, strcat(groupdir, sprintf('group%s%s_imageplane_%s_thresh%.2f_annulus.png',...
        lower(feature), outputROI, coordtoggle, groupthresh(1))));
end

close all; % close any open figure windows.

% for checking correct orientation- plot point in each of 4 quadrants.
% [x1,y1] = pol2cart(deg2rad(45),4); [x2,y2] = pol2cart(deg2rad(135),4);
% [x3,y3] = pol2cart(deg2rad(225),4); [x4,y4] = pol2cart(deg2rad(315),4); hold on
% plot(x1,y1,'r.'); plot(x2,y2,'y.'); plot(x3,y3,'g.'); plot(x4,y4,'b.');

%% --------- PLOT SINGULAR CONDITION DATA X ALL PARTICIPANTS ----------- %%

fig= figure('Position', get(0, 'Screensize'));  % open a full-screen figure for better visualisation.
[ha, pos] = tight_subplot(3,4,[0.03, 0.01],0.03,0.03); % create an arrangement of 12 tightly-tiled subplots.

for i = 1:length(participants) % for each participant-in turn:
    
    axes(ha(i)); % within a participant-specific subplot:
    data = plotfeaturemap(:,:,i); % extract the participant's normalised data.
    
    % plot the normalised interpolated data as a 3D surface viewed on a 2D plane.
    surf(XQ,YQ,data); shading interp; view(90,-90);
    colorbar; colormap(cmap); caxis([-6 6]); title(participants{i});
    xlim([-6 6]); ylim([-6 6]); set(gca, 'XTick', []); set(gca, 'YTick', []);
    
    if overlaytoggle == 1 % if we specified an RFP overlay, add this to the plot.
        hold on; plot3(distXQ, distYQ,  distances-10, 'Color', [1 1 1], 'LineWidth',2);
    end
end

% save figure as a .png image file.
fig.Color = 'white'; saveas(gca, [groupdir, sprintf('group%s%s_xparticipants_%s_thresh%.2f.png',...
    lower(feature), outputROI, coordtoggle, groupthresh(1))]);
close all; % close any open figure windows.

if annulustoggle == 1 % if we specified an annulus overlay, repeat this process adding the annulus overlay.
    fig= figure('Position', get(0, 'Screensize'));  % open a full-screen figure for better visualisation.
    [ha, pos] = tight_subplot(3,4,[0.03, 0.01],0.03,0.03); % create an arrangement of 12 tightly-tiled subplots.
    
    for i = 1:length(participants) % for each participant-in turn:
        
        axes(ha(i)); % within a participant-specific subplot:
        data = plotfeaturemap(:,:,i); % extract the participant's normalised data.
        
        % plot the normalised interpolated data as a 3D surface viewed on a 2D plane.
        surf(XQ,YQ,data); shading interp; view(90,-90);
        colorbar; colormap(cmap); caxis([-6 6]);
        xlim([-6 6]); ylim([-6 6]); set(gca, 'XTick', []); set(gca, 'YTick', []);
        title(participants{i});
        
        if overlaytoggle == 1 % if we specified an RFP overlay, add this to the plot.
            hold on; plot3(distXQ, distYQ,  distances-10, 'Color', [1 1 1], 'LineWidth',2);
        end
        
        hold on;
        annulusmean = mean(annulusmap(:,:,i),2)';
        theta = linspace(0,360,360); theta = deg2rad(theta);
        annulusmean = annulusmean+2;
        [x,y] = pol2cart(theta,annulusmean);
        
        plot3(x,y,ones(1,360)*-3,'Color', [0.5 0.5 0.5],'LineWidth', 3);
    end
    fig.Color = 'white'; saveas(gca, [groupdir, sprintf('group%s%s_xparticipants_%s_thresh%.2f_annulus.png',...
        lower(feature), outputROI, coordtoggle, groupthresh(1))]);
end

close all; % close any open figure windows.
%% --------------------------------------------------------------------- %%