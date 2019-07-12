%% ----------------- PLOT CONN DIFFERENCE CORRELATIONS ----------------- %% 

% loads in .mat files containing t-statistics and associated p-values for
% comparisons between attentional focus conditions. assesses these
% differences for significance and plots the resulting difference matrices.

% 11/7/19 KWN

clear; clc;
addpath(genpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions/'));

%% -------------------- SPECIFY ANALYSIS CONDITIONS -------------------- %%

% specify condition comparison names.
conds_full{1} = {'Passive-Orientation', 'Passive-Contrast', 'Passive-Shape'};
conds_full{2} = {'Orientation-Passive', 'Contrast-Passive', 'Shape-Passive'};
conds_full{3} = {'Orientation-Contrast', 'Orientation-Shape', 'Contrast-Shape'};
conds_full{4} = {'Contrast-Orientation', 'Shape-Orientation', 'Shape-Contrast'};

% specify condition comparisons generated in conn output.
condlabels_full{1} = {'ORIENTATION(-1).PASSIVE(1)', 'CONTRAST(-1).PASSIVE(1)', 'PASSIVE(1).SHAPE(-1)'};
condlabels_full{2} = {'ORIENTATION(1).PASSIVE(-1)', 'CONTRAST(1).PASSIVE(-1)', 'PASSIVE(-1).SHAPE(1)'};
condlabels_full{3} = {'CONTRAST(-1).ORIENTATION(1)', 'ORIENTATION(1).SHAPE(-1)', 'CONTRAST(1).SHAPE(-1)'};
condlabels_full{4} = {'CONTRAST(1).ORIENTATION(-1)', 'ORIENTATION(-1).SHAPE(1)', 'CONTRAST(-1).SHAPE(1)'};

% specify condition comparison output filename.
output_filename_full{1} = 'conn_connectivity_difference-passive-condition';
output_filename_full{2} = 'conn_connectivity_difference-conditions_inverse';
output_filename_full{3} = 'conn_connectivity_difference-conditions';
output_filename_full{4} = 'conn_connectivity_difference-conditions_inverse';

ROIs = {'V1', 'V3AB', 'hV4', 'LO1', 'LO2', 'IPS0'}; % specify ROIs.

cmap = fireice; % load hot-cold colourmap.

%% ------------- VISUALLY DISPLAY DIFFERENCE T-STATISTICS -------------- %%

for i = 1:length(conds_full) % for each condition comparison in turn:
    
    % select the condition specific analysis parameters. 
    conds = conds_full{i};
    condlabels = condlabels_full{i};
    output_filename = output_filename_full{i};
    
    % load a full-screen figure with three subplots.
    fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    ha = tight_subplot(1,3, [.05, .05], .08, .08);
    
    for thiscond = 1:length(conds) % for each condition comparison:
        
        % load the comparison specific .mat file.
        load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/conn_output/RFPachromatic/results/secondlevel/visualROIs/AllSubjects/%s/ROI.mat',...
            condlabels{thiscond}));
        
        % organise the t statistics to match the original correlation matrix structure.
        tmat = [ROI(1).F(1:6); ROI(2).F(1:6); ROI(3).F(1:6);...
            ROI(4).F(1:6); ROI(5).F(1:6); ROI(6).F(1:6)];
        tmat = tril(tmat); tmat(tmat == 0) = NaN; % set the upper half of the matrix to NaNs.
        
        % repeat the same process for the associated p value matrix.
        tmat_p = [ROI(1).p(1:6); ROI(2).p(1:6); ROI(3).p(1:6);...
            ROI(4).p(1:6); ROI(5).p(1:6); ROI(6).p(1:6)];
        tmat_p = tril(tmat_p); tmat_p(tmat_p == 0) = NaN;
        
        % threshold p values so that we only display significant t statistics.
        tmat_p_thresh = tmat_p <.05;
        tmat_thresh = tmat; tmat_thresh(tmat_p>.05) = 0;
        
        axes(ha(thiscond)); % plot the resulting t statistic difference matrix.
        imagesc(tmat_thresh, 'AlphaData', ~isnan(tmat_thresh), [-5, 5]);
        colormap(cmap); colorbar;
        set(gca, 'TickLength', [0 0]); title(conds{thiscond});
        xlim([0.5 length(ROIs)+0.5]); ylim([0.5 length(ROIs)+0.5]);
        set(gca, 'XTickLabel', ROIs); set(gca, 'YTickLabel', ROIs);
        xticks(1:length(ROIs)); xtickangle(45);
        
        clear ROI % clear stored data for next iteration.
    end
    
    % save figure as a .png file.
    saveas(gcf, sprintf('/scratch/groups/Projects/P1323/original/fMRI/conn_output/%s',output_filename), 'png'); % save the figure.
end
%% --------------------------------------------------------------------- %% 