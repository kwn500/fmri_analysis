%% ------------------------- PROCESS CONN DATA ------------------------- %% 

% extracts condition-specific correlation coefficients from conn analysis
% and plots these in the original formatting for comparison with manual
% connectivity analysis. 

% 11/7/19 KWN 

clear;clc;
addpath(genpath('/scratch/sg3/P1323/code/fmri_analysis/functions'));
cmap = fireice;

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

% specify the conditions to analyse (must match conn file naming). 
condition = 'feature';

if isequal(condition, 'feature')
    conds = {'ORIENTATION', 'CONTRAST', 'SHAPE', 'PASSIVE'};
elseif isequal(condition, 'colour-passive')
    conds = {'RG-ATTENTION', 'RG-PASSIVE', 'BY-ATTENTION', 'BY-PASSIVE',...
        'LUM-ATTENTION', 'LUM-PASSIVE'};
elseif isequal(condition, 'colourfeature')
    conds = {'RG-ORIENTATION', 'RG-CONTRAST', 'RG-SHAPE', 'RG-PASSIVE',...
        'BY-ORIENTATION', 'BY-CONTRAST', 'BY-SHAPE', 'BY-PASSIVE',...
        'LUM-ORIENTATION', 'LUM-CONTRAST', 'LUM-SHAPE', 'LUM-PASSIVE'};
elseif isequal(condition, '3feature')
    conds = {'FACE','ORIENTATION', 'COLOUR', 'SHAPE', 'PASSIVE'};
elseif isequal(condition, '3x3feature')
    conds = {'PASSIVE', 'VERTICAL', 'HORIZONTAL', 'DIAGONAL', 'RED',...
        'GREEN', 'BLUE', 'CIRCULAR', 'SQUARE', 'TRIANGULAR', 'FACE'};
end

% specify ROIs to analyse.
ROIs = {'V1', 'V3AB', 'hV4', 'LO1', 'LO2', 'IPS0'}; 

%% ------------ EXTRACT AND PLOT CORRELATION COEFFICIENTS -------------- %%

for thiscond = 1:length(conds) % for each condition:
    
    % extract the condition-specific data struct. 
    data{thiscond} = load(sprintf('/raid/kwn/P1361/fMRI/conn_output/3x3feature/3x3featureconn/results/secondlevel/visualROIs/AllSubjects/%s/ROI.mat',...
        upper(conds{thiscond})));
    
    for thisROI = 1:length(ROIs)  % for each ROI;
        
        % re-order the ROI-specific correlation matrices. 
        corrdata{thiscond,thisROI} = data{thiscond}.ROI(thisROI).y(:,1:6);
    end
end

% save this data to a .mat file.
save('/raid/kwn/P1361/fMRI/conn_output/3x3feature/correlation_data', 'data');

for thiscond = 1:length(conds) % for each condition:
    
    % stack the ROI correlation matrices in the order matching the original analysis. 
    cmat{thiscond} = [nanmean(corrdata{thiscond,1}); nanmean(corrdata{thiscond,2}); nanmean(corrdata{thiscond,3});...
        nanmean(corrdata{thiscond,4}); nanmean(corrdata{thiscond,5}); nanmean(corrdata{thiscond,6})]; 
    
    % remove the on-diagonals and set the upper half of the matrix to NaN.
    cmat{thiscond} = tril(cmat{thiscond}); cmat{thiscond}(cmat{thiscond} == 0) = NaN;
end

% stack all condition correlation matrices and calculate the global mean. 
stackedmat = cat(3, cmat{1}, cmat{2}, cmat{3}, cmat{4});
globalmean = nanmean(stackedmat, 3);

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
ha = tight_subplot(4,3, [.05, .05], .08, .08);

for thiscond = 1:length(conds) % for each condition
    
    % plot the correlation matrix. 
    axes(ha(thiscond));
    
    % remove the global mean and fisher transform the resulting values. 
    cmatfisher{thiscond} = atanh(cmat{thiscond}- globalmean); 
    imagesc(cmatfisher{thiscond}, 'AlphaData', ~isnan(cmatfisher{thiscond}), [-.1, .1]);
    colormap(cmap); colorbar;
    set(gca, 'TickLength', [0 0]); title(conds{thiscond});
    xlim([0.5 length(ROIs)+0.5]); ylim([0.5 length(ROIs)+0.5]);
    set(gca, 'XTickLabel', ROIs); set(gca, 'YTickLabel', ROIs); 
    xticks(1:length(ROIs)); xtickangle(45);
    
    % vectorise the fisher-transformed correlation coefficients.
    cmatfisher{thiscond}(isnan(cmatfisher{thiscond})) = [];
end

% save the figure as a .png file. 
saveas(gcf, '/raid/kwn/P1361/fMRI/conn_output/3x3feature/conn_connectivity', 'png'); % save the figure. 
%% --------------------------------------------------------------------- %%