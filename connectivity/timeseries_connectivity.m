%% ------------------ TIMESERIES CONNECTIVITY ANALYSIS ----------------- %%

% loads in ROI- and condition-specific timeseries data, which we
% pre-process (i.e. mean removal), and correlate across ROIs. we plot these
% fisher transformed correlation coefficients at the individual and group
% levels and calculate the distance between ROI activation as a function of
% attentional condition.

% 18/2/2019

clear; clc; close all;
addpath(genpath('/Users/kirstie/Documents/code/fmri_analysis/functions'));

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

params.expcond = 'original'; % project (original, colour, naturalistic)
params.cond = 'block'; % analysis (block, colour, feature, na)
params.condition_toggle = 1; % specify toggle for specific interaction analyses.

params.iterations = 10000; % specify desired iterations.

% specify and create output directory.
params.outputdir = (sprintf('/Users/kirstie/Documents/analysis/%s/%s/connectivity/',...
    params.expcond, params.cond));
[~,~] = mkdir(params.outputdir);

% specify ROIs to analyse.
params.ROIs = {'V1', 'V3AB',  'V4', 'LO1', 'LO2', 'IPS0'};

% specify analysis-specific condition labels.
if isequal(params.cond, 'feature') || isequal(params.cond, 'block')
    params.condlabels = {'Orientation', 'Contrast', 'Shape', 'Passive'};
elseif isequal(params.cond,'colour')
    params.condlabels = {'RedGreenAttention', 'RedGreenPassive',...
        'BlueYellowAttention', 'BlueYellowPassive', 'LuminanceAttention', 'LuminancePassive'};
    
    % if this is an interaction analysis, specify the individual interaction
    % conditions to analyse.
elseif isequal(params.cond, 'colourxfeature')
    if params.condition_toggle == 1
        params.condlabels = {'RGOrientation', 'BYOrientation', 'LUMOrientation'};
    elseif params.condition_toggle == 2
        params.condlabels = {'RGContrast', 'BYContrast', 'LUMContrast'};
    elseif params.condition_toggle == 3
        params.condlabels = {'RGShape', 'BYShape', 'LUMShape'};
    elseif params.condition_toggle == 4
        params.condlabels = {'RGOrientation', 'RGContrast', 'RGShape', 'RGPassive'};
    elseif params.condition_toggle == 5
        params.condlabels = {'BYOrientation', 'BYContrast', 'BYShape', 'BYPassive'};
    elseif params.condition_toggle == 6
        params.condlabels = {'LUMOrientation', 'LUMContrast', 'LUMShape', 'LUMPassive'};
    end
    
elseif isequal(params.expcond,'naturalistic') && isequal(params.cond, '3feature')
    params.condlabels = {'Face', 'Orientation', 'Contrast', 'Shape', 'Passive'};
elseif isequal(params.expcond, 'naturalistic') && isequal(params.cond, '3x3feature')
    if params.condition_toggle == 1
        params.condlabels = {'Vertical', 'Horizontal', 'Diagonal'};
    elseif params.condition_toggle == 2
        params.condlabels = {'Red', 'Green', 'Blue'};
    elseif params.condition_toggle == 3
        params.condlabels = {'Circular', 'Square', 'Triangular'};
    end
end

% specify analysis-specific participant lists.
if isequal(params.expcond, 'original')
    pps = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455',...
        'R3517','R3773', 'R3932', 'R4065', 'R4127', 'R4496'};
elseif isequal(params.expcond, 'colour')
    pps = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3517',...
        'R3773', 'R4059', 'R4065', 'R4244', 'R4831', 'R4928'};
elseif isequal(params.expcond, 'naturalistic') && isequal(params.cond, '3feature')
    pps = {'R2548', 'R2904', 'R3111', 'R3517', 'R4059', 'R4065',...
        'R4127', 'R4244', 'R4831', 'R4833', 'R4890', 'R4928'};
elseif isequal(params.expcond, 'naturalistic') && isequal(params.cond, '3x3feature')
    pps = {'R2548', 'R2904', 'R3111', 'R3517', 'R3773', 'R4059',...
        'R4065', 'R4127', 'R4244', 'R4831', 'R4833', 'R4890', 'R4928', 'R5006'};
end

%% ---------------- CALCULATE CORRELATION COEFFICIENTS ----------------- %%

for thispp = 1:length(pps) % for each participant in turn:
    
    clear dataList % refresh timeseries storage.
    
    for thisROI = 1:length(params.ROIs) % for each ROI in turn:
        
        % load participant & ROI processed timeseries data.
        if isequal(params.expcond, 'colour')
            load(sprintf('/Users/kirstie/Documents/analysis/%s/%s/connectivity/processed_timeseries/%s/%s_%s_timeseries.mat',...
                params.expcond, params.cond, params.ROIs{thisROI}, pps{thispp}, params.ROIs{thisROI}));
        elseif isequal(params.expcond, 'original')
            load(sprintf('/Users/kirstie/Documents/analysis/original/connectivity/processed_timeseries/%s/%s_%s_timeseries.mat',...
                params.ROIs{thisROI}, pps{thispp}, params.ROIs{thisROI}));
        elseif isequal(params.expcond, 'naturalistic')
            load(sprintf('/Users/kirstie/Documents/analysis/naturalistic/3feature/connectivity/processed_timeseries/%s/%s_%s_timeseries.mat',...
                params.ROIs{thisROI}, pps{thispp}, params.ROIs{thisROI}));
        end
        
        % calculate voxel mean for each attention condition.
        if isequal(params.cond, 'feature') || isequal(params.cond, 'block')
            tSeriesvoxels = struct2cell(tSeries);
            for thisfield = 1:size(tSeriesvoxels,1)-1
                dataList(thisfield,:,thisROI) = mean(tSeriesvoxels{thisfield},2);
            end
            
            % if this is a colour analysis condition:
        elseif isequal(params.cond, 'colour')
            
            % not all conditions have the same number of entries (different
            % numbers of condition repetitions) so we bootstrap these
            % values from the condition data for that participant for a
            % specified number of samples and calculate the mean so that
            % all conditions have the same amount of data to analyse.
            maxsize = size(tSeries.rga,1);
            minsize = size(tSeries.rgp,1);
            bstrp_size = maxsize-minsize;
            data_values = [1:minsize]';
            samp_number = floor(minsize/2);
            
            count = minsize + 1;
            
            % for each of the specified bootstrapped samples in turn:
            for i = 1:bstrp_size
                
                % sample from the original data with replacement.
                rgpsamp(i,:) = randsample(data_values,samp_number, 'true');
                bypsamp(i,:) = randsample(data_values,samp_number, 'true');
                lumpsamp(i,:) = randsample(data_values,samp_number,'true');
                
                tSeries.rgp(count,:) = mean(tSeries.rgp(rgpsamp(i,:),:));
                tSeries.byp(count,:) = mean(tSeries.byp(bypsamp(i,:),:));
                tSeries.lump(count,:) = mean(tSeries.lump(lumpsamp(i,:),:));
                
                count = count + 1; % increment bootstrap counter variable.
            end
            clear rgpsamp bypsamp lumpsamp % clear bootstrapped variables.
            
            % calculate voxel mean for each attention condition.
            tSeriesvoxels = struct2cell(tSeries);
            for thisfield = 1:size(tSeriesvoxels,1)
                dataList(thisfield,:,thisROI) = mean(tSeriesvoxels{thisfield},2);
            end
            
            % if this is an interaction analysis, extract the
            % condition-specific data for further processing.
        elseif isequal(params.cond, 'colourxfeature')
            if params.condition_toggle == 1
                tSeriesnew.rgo = tSeries.rgo;
                tSeriesnew.byo = tSeries.byo;
                tSeriesnew.lumo = tSeries.lumo;
            elseif params.condition_toggle == 2
                tSeriesnew.rgc = tSeries.rgc;
                tSeriesnew.byc = tSeries.byc;
                tSeriesnew.lumc = tSeries.lumc;
            elseif params.condition_toggle == 3
                tSeriesnew.rgs = tSeries.rgs;
                tSeriesnew.bys = tSeries.bys;
                tSeriesnew.lums = tSeries.lums;
            elseif params.condition_toggle == 4
                tSeriesnew.rgo = tSeries.rgo;
                tSeriesnew.rgc = tSeries.rgc;
                tSeriesnew.rgs = tSeries.rgs;
                tSeriesnew.rgp = tSeries.rgp;
            elseif params.condition_toggle == 5
                tSeriesnew.byo = tSeries.byo;
                tSeriesnew.byc = tSeries.byc;
                tSeriesnew.bys = tSeries.bys;
                tSeriesnew.byp = tSeries.byp;
            elseif params.condition_toggle == 6
                tSeriesnew.lumo = tSeries.lumo;
                tSeriesnew.lumc = tSeries.lumc;
                tSeriesnew.lums = tSeries.lums;
                tSeriesnew.lump = tSeries.lump;
            end
            
            clear tSeries;
            tSeries = tSeriesnew;
            
            % calculate voxel mean for each attention condition.
            tSeriesvoxels = struct2cell(tSeries);
            for thisfield = 1:size(tSeriesvoxels,1)
                dataList(thisfield,:,thisROI) = mean(tSeriesvoxels{thisfield},2);
            end
        elseif isequal(params.cond, '3feature')
            
            % calculate the minimum and maximum sizes of data in terms of
            % condition repetitions for the naturalistic data.
            maxsize = max([size(tSeries.face,1), size(tSeries.orientation,1),size(tSeries.colour,1),size(tSeries.shape,1),size(tSeries.passive,1)]);
            minsize = min([size(tSeries.face,1), size(tSeries.orientation,1),size(tSeries.colour,1),size(tSeries.shape,1),size(tSeries.passive,1)]);
            
            fields = fieldnames(tSeries);
            
            % as within the interaction analyses, we sample with
            % replacement over a specified number of bootstrapped
            % iterations to ensure each condition has the same amount of
            % data.
            for i = 1:numel(fieldnames(tSeries))
                data = tSeries.(fields{i});
                bstrp_size = maxsize-size(data,1);
                data_values = [1:size(data,1)]';
                samp_number = floor(size(data,1)/2);
                
                count = size(data,1) + 1;
                for x = 1:bstrp_size
                    datasamp(x,:) = randsample(data_values,samp_number,'true');
                    tSeries.(fields{i})(count,:) = mean(tSeries.(fields{i})(datasamp(x,:),:));
                    count = count + 1;
                end
                clear datasamp data_values samp_number
            end
            tSeriesvoxels = struct2cell(tSeries);
            for thisfield = 1:size(tSeriesvoxels,1)
                dataList(thisfield,:,thisROI) = mean(tSeriesvoxels{thisfield},2);
            end
        elseif isequal(params.cond, '3x3feature')
            
            % repeat the same bootstrapping procedure as above.
            structfields = fieldnames(tSeries);
            for i = 1:numel(structfields)
                maxsize(i) = max(size(tSeries.(structfields{i}),1));
                minsize(i) = min(size(tSeries.(structfields{i}),1));
            end
            maxsize = max(maxsize);
            minsize = min(minsize);
            
            fields = fieldnames(tSeries);
            
            for i = 1:numel(fieldnames(tSeries))
                data = tSeries.(fields{i});
                bstrp_size = maxsize-size(data,1);
                data_values = [1:size(data,1)]';
                samp_number = floor(size(data,1)/2);
                
                count = size(data,1) + 1;
                for x = 1:bstrp_size
                    datasamp(x,:) = randsample(data_values,samp_number,'true');
                    tSeries.(fields{i})(count,:) = mean(tSeries.(fields{i})(datasamp(x,:),:));
                    count = count + 1;
                end
                clear datasamp data_values samp_number
            end
            
            if params.condition_toggle == 1
                tSeriesnew.vertical = tSeries.vertical;
                tSeriesnew.horizontal = tSeries.horizontal;
                tSeriesnew.diagonal = tSeries.diagonal;
            elseif params.condition_toggle == 2
                tSeriesnew.red = tSeries.red;
                tSeriesnew.green = tSeries.green;
                tSeriesnew.blue = tSeries.blue;
            elseif params.condition_toggle == 3
                tSeriesnew.circular = tSeries.circular;
                tSeriesnew.square = tSeries.square;
                tSeriesnew.triangular = tSeries.triangular;
            end
            
            clear tSeries;
            tSeries = tSeriesnew;
            
            % calculate voxel mean for each attention condition.
            tSeriesvoxels = struct2cell(tSeries);
            for thisfield = 1:size(tSeriesvoxels,1)
                dataList(thisfield,:,thisROI) = mean(tSeriesvoxels{thisfield},2);
            end
            
        end
    end
    
    ROIMean=squeeze(mean(dataList,3)); % calculate mean across ROIs.
    condMean=squeeze(mean(ROIMean)); % calculate mean across conditions.
    
    % fit and remove the condMean from the condition and ROI-specific data.
    for thisROI = 1:length(params.ROIs)
        for thiscond=1:length(params.condlabels)
            fittedAmp=squeeze(dataList(thiscond,:,thisROI))/condMean;
            dataList(thiscond,:,thisROI)=dataList(thiscond,:,thisROI)-condMean.*fittedAmp;
        end
    end
    
    % correlate each pairwise combination of ROIs timeseries for each
    % condition and participant.
    for thiscond=1:size(dataList,1)
        [output.ccMat(thiscond,thispp,:,:),output.pVal(thiscond,thispp,:,:)]...
            = corr(squeeze(dataList(thiscond,:,:)),'Type','Kendall');
    end
end

% global mean correlation coefficient across conditions and participants.
output.globalmean = squeeze(mean(squeeze(mean(output.ccMat))));

% save the correlation coefficients and p-values.
original = output; save(strcat(params.outputdir, 'correlation-pvalues.mat'), 'original');

% select only significant p-values for group correlation matrices.
output.pVal(output.pVal>.05 & output.pVal ~= 1)=NaN; output.pVal(output.pVal<=.05)=1;

% multiply coefficients by p-value significance.
% output.ccMat=output.ccMat.*output.pVal;

%% ------------- FISHER TRANSFORM CORRELATION COEFFICIENTS ------------- %%

% extract the global mean across all voxels and conditions and participants
% for fisher transformation.
output.globalmean = squeeze(mean(squeeze(mean(output.ccMat))));

for thiscond = 1:size(dataList,1) % for each condition in turn:
    
    % calculate the condition-specific mean across voxels and participants.
    output.corr.means(thiscond,:,:)=squeeze(nanmean(squeeze(output.ccMat(thiscond,:,:,:))));
    
    % remove the global mean from this mean data and fisher transform the
    % coefficients.
    output.corr.fisherdata(thiscond,:,:) = atanh(squeeze(output.corr.means(thiscond,:,:))-output.globalmean);
    
    % extract the lower diagonal of data for plotting.
    plotdata = squeeze(output.corr.fisherdata(thiscond,:,:));
    idx = logical(tril(ones(size(plotdata)),-1));
    plotdata_edit(:,thiscond) = plotdata(idx);
end

%% ------------------- CORRELATION ACROSS CONDITIONS ------------------- %%

% for each combination of conditions in turn:
for thiscond1 = 1:size(plotdata_edit,2)
    for thiscond2 = 1:size(plotdata_edit,2)
        
        % correlate the fisher-transformed correlation coefficients across
        % ROIs.
        [r(thiscond1,thiscond2),p(thiscond1,thiscond2)] =...
            corr(plotdata_edit(:,thiscond1),plotdata_edit(:,thiscond2));
    end
end

% save the correlation data across conditions as a .mat file.
save(strcat(params.outputdir, 'correlation-acrossconditions.mat'), 'r','p');

%% ---------------------- REPEATED MEASURES ANOVA ---------------------- %%

[means, stderrs] = deal([]); % initialise storage variables.

% for each condition in turn:
for thiscond = 1:length(params.condlabels)
    
    % calculate the mean and standard error for each condition across ROIs.
    means = [means, mean(plotdata_edit(:,thiscond))];
    stderrs = [stderrs, std(plotdata_edit(:, thiscond)/sqrt(length(plotdata_edit(:,thiscond))))];
end

% create a ROI comparison column.
ppnumber = repmat({'ROIcomparison'}, size(plotdata_edit(:,1),1),1);

% perform a one-way repeated measures ANOVA comparing activation across all
% attention columns.
if isequal(params.cond, 'feature') || isequal(params.cond, 'block') || isequal(params.cond, 'colourxfeature') && params.condition_toggle > 3
    t = table(ppnumber, plotdata_edit(:,1), plotdata_edit(:,2), plotdata_edit(:,3), plotdata_edit(:,4),...
        'VariableNames', {'ROIcomparison', 'Orientation', 'Contrast', 'Shape', 'Passive'});
    Meas = table([1 2 3 4]', 'VariableNames', {'Connectivity'});
    rm = fitrm(t, 'Orientation-Passive~1', 'WithinDesign', Meas);
elseif isequal(params.cond, 'colour')
    t = table(ppnumber, plotdata_edit(:,1), plotdata_edit(:,2), plotdata_edit(:,3), plotdata_edit(:,4),...
        plotdata_edit(:,5), plotdata_edit(:,6), 'VariableNames', {'ROIcomparison',...
        'RedGreenAttention', 'RedGreenPassive', 'BlueYellowAttention',...
        'BlueYellowPassive', 'LuminanceAttention', 'LuminancePassive'});
    Meas = table([1 2 3 4 5 6]', 'VariableNames', {'Connectivity'});
    rm = fitrm(t, 'RedGreenAttention-LuminancePassive~1', 'WithinDesign', Meas);
elseif isequal(params.cond, '3x3feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle < 4
    t = table(ppnumber, plotdata_edit(:,1), plotdata_edit(:,2), plotdata_edit(:,3),...
        'VariableNames', {'ROIcomparison', 'VAR1', 'VAR2', 'VAR3'});
    Meas = table([1 2 3]', 'VariableNames', {'Connectivity'});
    rm = fitrm(t, 'VAR1-VAR3~1', 'WithinDesign', Meas);
elseif isequal(params.cond, '3feature')
    t = table(ppnumber, plotdata_edit(:,1), plotdata_edit(:,2), plotdata_edit(:,3), plotdata_edit(:,4),...
        plotdata_edit(:,5), 'VariableNames', {'ROIcomparison',...
        'Face', 'Orientation', 'Colour', 'Shape', 'Passive'});
    Meas = table([1 2 3 4 5]', 'VariableNames', {'Connectivity'});
end

% perform mauchly's test of sphericity and bonferroni-corrected post-hoc
% tests.
mauchlyoutput = mauchly(rm);
ranovatbl = ranova(rm);
posthoc = multcompare(rm, 'Connectivity', 'ComparisonType', 'bonferroni');

% extract the unique post-hoc significance values for plotting.
if isequal(params.cond, 'feature') || isequal(params.cond, 'block') || isequal(params.cond, 'colourxfeature') && params.condition_toggle > 3
    posthoc_edit = [posthoc(1,:);posthoc(2,:);posthoc(3,:);posthoc(5,:);posthoc(6,:);posthoc(9,:)];
    
    P = nan(numel(means), numel(means));
    P(1,2) = cell2mat(table2cell(posthoc_edit(1,5)));
    P(1,3) = cell2mat(table2cell(posthoc_edit(2,5)));
    P(1,4) = cell2mat(table2cell(posthoc_edit(3,5)));
    P(2,3) = cell2mat(table2cell(posthoc_edit(4,5)));
    P(2,4) = cell2mat(table2cell(posthoc_edit(5,5)));
    P(3,4) = cell2mat(table2cell(posthoc_edit(6,5)));
    
    % specify bar-chart colours.
    colour =[.3 .3 .3; .5 .5 .5; .7 .7 .7; .9 .9 .9];
    
    % repeat the same process for the remaining analysis conditions.
elseif isequal(params.cond, 'colour')
    posthoc_edit = [posthoc(1,:);posthoc(2,:);posthoc(3,:);posthoc(4,:);posthoc(5,:);...
        posthoc(7,:);posthoc(8,:);posthoc(9,:);posthoc(10,:);posthoc(13,:);posthoc(14,:);...
        posthoc(15,:);posthoc(19,:);posthoc(20,:);posthoc(25,:)];
    
    P = nan(numel(means), numel(means));
    P(1,2) = cell2mat(table2cell(posthoc_edit(1,5)));
    P(1,3) = cell2mat(table2cell(posthoc_edit(2,5)));
    P(1,4) = cell2mat(table2cell(posthoc_edit(3,5)));
    P(1,5) = cell2mat(table2cell(posthoc_edit(4,5)));
    P(1,6) = cell2mat(table2cell(posthoc_edit(5,5)));
    P(2,3) = cell2mat(table2cell(posthoc_edit(6,5)));
    P(2,4) = cell2mat(table2cell(posthoc_edit(7,5)));
    P(2,5) = cell2mat(table2cell(posthoc_edit(8,5)));
    P(2,6) = cell2mat(table2cell(posthoc_edit(9,5)));
    P(3,4) = cell2mat(table2cell(posthoc_edit(10,5)));
    P(3,5) = cell2mat(table2cell(posthoc_edit(11,5)));
    P(3,6) = cell2mat(table2cell(posthoc_edit(12,5)));
    P(4,5) = cell2mat(table2cell(posthoc_edit(13,5)));
    P(4,6) = cell2mat(table2cell(posthoc_edit(14,5)));
    P(5,6) = cell2mat(table2cell(posthoc_edit(15,5)));
    
    colour =[.3 .3 .3];
elseif isequal(params.cond, '3x3feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle < 4
    posthoc_edit = [posthoc(1,:);posthoc(2,:);posthoc(4,:)];
    
    P = nan(numel(means), numel(means));
    P(1,2) = cell2mat(table2cell(posthoc_edit(1,5)));
    P(1,3) = cell2mat(table2cell(posthoc_edit(2,5)));
    P(2,3) = cell2mat(table2cell(posthoc_edit(3,5)));
    
    colour =[.3 .3 .3; .5 .5 .5; .7 .7 .7; .9 .9 .9];
elseif isequal(params.cond, '3feature')
    posthoc_edit = [posthoc(1,:);posthoc(2,:);posthoc(3,:);posthoc(4,:);posthoc(6,:);posthoc(7,:);posthoc(8,:);...
        posthoc(11,:); posthoc(12,:);posthoc(16,:)];
    
    P = nan(numel(means), numel(means));
    P(1,2) = cell2mat(table2cell(posthoc_edit(1,5)));
    P(1,3) = cell2mat(table2cell(posthoc_edit(2,5)));
    P(1,4) = cell2mat(table2cell(posthoc_edit(3,5)));
    P(1,5) = cell2mat(table2cell(posthoc_edit(4,5)));
    P(2,3) = cell2mat(table2cell(posthoc_edit(5,5)));
    P(2,4) = cell2mat(table2cell(posthoc_edit(6,5)));
    P(2,5) = cell2mat(table2cell(posthoc_edit(7,5)));
    P(3,4) = cell2mat(table2cell(posthoc_edit(8,5)));
    P(3,5) = cell2mat(table2cell(posthoc_edit(9,5)));
    P(4,5) = cell2mat(table2cell(posthoc_edit(10,5)));
    
    colour =[.3 .3 .3];
end

% plot the data as a bar chart with standard error bars and associated
% significance values.
P(P>.05) = NaN;
PT = P'; lidx = tril(true(size(P)), -1); P(lidx) = PT(lidx);

figure(1)
superbar(means, 'E', stderrs, 'P', P,'BarFaceColor', colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2,'PStarShowNS', false);

% apply condition-specific labels and axis scaling.
if isequal(params.cond, 'feature') || isequal(params.cond, 'block') || isequal(params.cond, 'colourxfeature') && params.condition_toggle > 3
    ylim([-.025 .05])
    set(gca, 'xtick', 1:4, 'xticklabel', {'Orientation', 'Contrast', 'Shape', 'Passive'});
elseif isequal(params.cond, 'colour')
    ylim([-.12 .15]);
    set(gca, 'xtick', 1:6, 'xticklabel', {'RedGreenAttention', 'RedGreenPassive',...
        'BlueYellowAttention', 'BlueYellowPassive', 'LuminanceAttention', 'LuminancePassive'});
    xtickangle(-45)
elseif isequal(params.cond, '3x3feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle < 4
    ylim([-.025 .05])
    set(gca, 'xtick', 1:3, 'xticklabel', params.condlabels);
end

% save the figure as a .pdf file and rmanova analysis output as a .mat
% file.
saveas(gcf, strcat(params.outputdir, sprintf('connectivity_barchart_%s.pdf', params.cond)));
save(strcat(params.outputdir, sprintf('connectivity_rmanova_%s.mat', params.cond)),...
    'mauchlyoutput', 'ranovatbl', 'posthoc_edit');

%% -------------------- GROUP CORRELATION MATRICES --------------------- %%

if isequal(params.expcond, 'original') || isequal(params.cond, 'feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle > 3
    fig = figure(1); [ha] = tight_subplot(1,4,[0.05, 0.05],0.08,0.08);
elseif isequal(params.cond, 'colour')
    fig = figure(1); [ha] = tight_subplot(1,6,[0.05, 0.05],0.08,0.08);
elseif isequal(params.cond, '3x3feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle < 4
    fig = figure(1); [ha] = tight_subplot(1,3,[0.05, 0.05],0.08,0.08);
elseif isequal(params.expcond, 'naturalistic') && isequal(params.cond, '3feature')
    fig = figure(1); [ha] = tight_subplot(1,5,[0.05, 0.05],0.08,0.08);
end

cmap=fireice; % load hot-cold colourmap.

for thiscond = 1:size(dataList,1) % for each condition:
    
    % calculate mean and standard deviation across participants.
    output.corr.means(thiscond,:,:)=squeeze(nanmean(squeeze(original.ccMat(thiscond,:,:,:))));
    output.corr.stds(thiscond,:,:)=squeeze(nanstd(squeeze(original.ccMat(thiscond,:,:,:))));
    
    % subtract global mean and fisher transform.
    output.corr.fisherdata(thiscond,:,:) = atanh(squeeze(output.corr.means(thiscond,:,:))-output.globalmean);
    axes(ha(thiscond));
    
    % retain lower-half of correlation matrix (other values set to NaN for
    % clearer plotting).
    plotdata = squeeze(output.corr.fisherdata(thiscond,:,:));
    plotdata = tril(plotdata); plotdata(plotdata==0)=NaN;
    imagesc(plotdata,'AlphaData',~isnan(plotdata),[-.06,.06]);
    colormap(cmap); colorbar; axis equal;
    set(gca,'TickLength',[0 0]); title(params.condlabels{thiscond});
    xlim([0.5 length(params.ROIs)+0.5]); ylim([0.5 length(params.ROIs)+0.5]);
    set(gca, 'XTickLabel', params.ROIs); set(gca, 'YTickLabel', params.ROIs);
    xticks(1:length(params.ROIs));xtickangle(45)
    
    output.corr.plotdata = squeeze(output.corr.fisherdata(thiscond,:,:));
    output.corr.plotdata = tril(output.corr.plotdata);
    output.corr.plotdata(output.corr.plotdata==0)=NaN;
    output.corr.groupmatrices(thiscond,:,:) = output.corr.plotdata;
end

% save group correlation matrices to a .pdf file.
fig.Color = 'white'; saveas(gca, strcat(params.outputdir, sprintf('timeseries_correlations-fisher.pdf')));

%% ----------------- INDIVIDUAL CORRELATION MATRICES ------------------- %%

% repeat the same process as above on an individual basis.
for thispp = 1:size(original.ccMat,2)
    for thiscond = 1:size(original.ccMat,1)
        
        % extract condition and participant specific unedited correlation coefficients.
        indcorr.data(thiscond,thispp,:,:) = squeeze(original.ccMat(thiscond,thispp,:,:));
        
        % calculate global mean across all conditions for each participant.
        indcorr.ppglobalmean(thispp,:,:) = squeeze(mean(original.ccMat(:,thispp,:,:)));
        
        % subtract global mean from each condition and fisher transform.
        indcorr.fisherdata(thiscond,thispp,:,:) = atanh(squeeze(indcorr.data(thiscond,thispp,:,:))...
            -squeeze(indcorr.ppglobalmean(thispp,:,:)));
        
        % reshape data to retain only upper half of correlation matrix.
        splitdata = squeeze(indcorr.fisherdata(thiscond,thispp,:,:));
        idx = logical(triu(ones(size(splitdata)),1));
        indcorr.fisherdatasplit(thiscond,thispp,:,:) = splitdata(idx)';
    end
    
    % calculate correlations between conditions for each participant.
    [indcorr.ccMat(thispp,:,:), indcorr.pVal(thispp,:,:)] = corr(squeeze(indcorr.fisherdatasplit(:,thispp,:,:))','Type','Kendall');
end

% for each pairwise condition comparison across participants, assess if the correlation
% coefficient is significantly greater than zero.
for thiscond1 = 1:size(indcorr.ccMat,2)
    for thiscond2 = 1:size(indcorr.ccMat,2)
        if thiscond1 == thiscond2
            continue
        end
        [indcorr.wsr.p(thiscond1, thiscond2),~,indcorr.wsr.stats(thiscond1,thiscond2)]...
            = signrank(squeeze(indcorr.ccMat(:,thiscond1,thiscond2)),0,'method','approximate');
    end
end

% save the important analysis output in a struct for later writing to file.
if isequal(params.cond, 'block') || isequal(params.cond, 'feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle > 3
    indcorr.summaryp = [indcorr.wsr.p(1,[2:end]), indcorr.wsr.p(2,[3:end]), indcorr.wsr.p(3,end)];
elseif isequal(params.cond, 'colourxfeature') && params.condition_toggle < 4 || isequal(params.cond, '3x3feature')
    indcorr.summaryp = [indcorr.wsr.p(1,[2:end]), indcorr.wsr.p(2,end)];
elseif isequal(params.cond, '3feature')
    indcorr.summaryp = [indcorr.wsr.p(1,[2:end]), indcorr.wsr.p(2,[3:end]), indcorr.wsr.p(3,[4:end]), indcorr.wsr.p(4,end)];
end
[~,~,~,indcorr.adj_p]=fdr_bh(indcorr.summaryp);

% specify figure for individual correlation matrices.
if isequal(params.expcond, 'original') || isequal(params.cond, 'feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle > 3
    fig = figure(3); [ha] = tight_subplot(1,4,[0.05, 0.05],0.08,0.08);
elseif isequal(params.cond, 'colour')
    fig = figure(3); [ha] = tight_subplot(2,3,[0.05, 0.05],0.08,0.08);
elseif isequal(params.cond, '3x3feature') || isequal(params.cond, 'colourxfeature') && params.condition_toggle < 4
    fig = figure(3); [ha] = tight_subplot(1,3,[0.05, 0.05],0.08,0.08);
elseif isequal(params.expcond, 'naturalistic') && isequal(params.cond, '3feature')
    fig = figure(3); [ha] = tight_subplot(1,5,[0.05, 0.05],0.08,0.08);
end

for thispp = 1:size(indcorr.fisherdata,2)
    for thiscond = 1:size(indcorr.fisherdata,1)
        plotdata = squeeze(indcorr.fisherdata(thiscond,thispp,:,:));
        axes(ha(thiscond));
        
        % retain lower-half of correlation matrix (other values set to NaN for
        % clearer plotting).
        plotdata = tril(plotdata); plotdata(plotdata==0)=NaN;
        imagesc(plotdata,'AlphaData',~isnan(plotdata),[-.1,.1]);
        colormap(cmap); colorbar; axis equal;
        set(gca,'TickLength',[0 0]); title(params.condlabels{thiscond});
        xlim([0.5 length(params.ROIs)+0.5]); ylim([0.5 length(params.ROIs)+0.5]);
        set(gca, 'XTickLabel', params.ROIs); set(gca, 'YTickLabel', params.ROIs);xticks(1:length(params.ROIs));
    end
    
    % save figure to .pdf file.
    fig.Color = 'white'; saveas(gca, strcat(params.outputdir,...
        sprintf('timeseries_correlations-individual_%s.pdf', pps{thispp})));
end

%%  CORRELATION WITH THRESHOLD DETECTION SCORES (COGNITIVE FLEXIBILITY)  %%

% if we are analysing data from our original experiment:
if isequal(params.expcond, 'original')
    
    % load participants' feature-specific 75% correct detection thresholds.
    thresholddata = xlsread('/Users/kirstie/Documents/analysis/colour/feature/connectivity/psychophysics_thresholds_feature.xlsx');
    
    % thresholddata(11,:) = []; % remove data from participant(s) excluded from analysis.
    
    % extract condition-specific correlation coefficients (correlation with passive).
    thresholds.op.data = squeeze(indcorr.ccMat(:,4,1));
    thresholds.cp.data = squeeze(indcorr.ccMat(:,4,2));
    thresholds.sp.data = squeeze(indcorr.ccMat(:,4,3));
    
    % calculate correlations between condition-passive correlation coefficients
    % and threshold detection scores.
    [thresholds.op.ccMat, thresholds.op.pVal] = corr(thresholds.op.data, thresholddata(:,3), 'Type', 'Kendall');
    [thresholds.cp.ccMat, thresholds.cp.pVal] = corr(thresholds.cp.data, thresholddata(:,4), 'Type', 'Kendall');
    [thresholds.sp.ccMat, thresholds.sp.pVal] = corr(thresholds.sp.data, thresholddata(:,5), 'Type', 'Kendall');
end

%% -------------------- EUCLIDEAN DISTANCE ANALYSIS -------------------- %%

% specify ROIs to analyse.
cond.ROIcomps = {'V3AB', 'V1', 'IPS0', 'V4', 'LO1', 'LO2'};

if isequal(params.cond, 'block') || isequal(params.cond, 'feature')
    
    % specify numeric values corresponding to desired condition pairs to analyse.
    comps = [1 2 3 5 6 9];
    
    % specify condition-specific labels for axis naming.
    labels = {'O-C', 'O-S', 'O-P', 'C-S', 'C-P', 'S-P'};
    
elseif isequal(params.cond, 'colour')
    comps = [1 2 4 13 14 25];
    labels = {'RGA-RGP', 'RGA-BYA', 'RGA-LUMA', 'BYA-BYP', 'BYA-LUMA', 'LUMA-LUMP'};
elseif isequal(params.cond, 'colourxfeature')
    
    if params.condition_toggle < 4
        comps = [1 2 4];
    elseif params.condition_toggle > 3
        comps =  [1 2 3 5 6 9];
    end
    if params.condition_toggle == 1
        labels = {'RGO-BYO', 'RGO-LUMO', 'BYO-LUMO'};
    elseif params.condition_toggle == 2
        labels = {'RGC-BYC', 'RGC-LUMC', 'BYC-LUMC'};
    elseif params.condition_toggle == 3
        labels = {'RGS-BYS', 'RGS-LUMS', 'BYS-LUMS'};
    elseif params.condition_toggle == 4
        labels = {'RGO-RGC', 'RGO-RGS', 'RGO-RGP', 'RGC-RGS', 'RGC-RGP', 'RGS-RGP'};
    elseif params.condition_toggle == 5
        labels = {'BYO-BYC', 'BYO-BYS', 'BYO-BYP', 'BYC-BYS', 'BYC-BYP', 'BYS-BYP'};
    elseif params.condition_toggle == 6
        labels = {'LUMO-LUMC', 'LUMO-LUMS', 'LUMO-LUMP', 'LUMC-LUMS', 'LUMC-LUMP', 'LUMS-LUMP'};
    end
elseif isequal(params.expcond, 'naturalistic') && isequal(params.cond, '3feature')
    comps = [1 2 3 4 6 7 8 11 12 16];
    labels = {'F-O', 'F-C', 'F-S', 'F-P', 'O-C', 'O-S', 'O-P', 'C-S', 'C-P', 'S-P'};
elseif isequal(params.cond, '3x3feature')
    comps = [1 2 4];
    if params.condition_toggle == 1
        labels = {'V-H','V-D','H-D'};
    elseif params.condition_toggle == 2
        labels = {'R-G', 'R-B','G-B'};
    elseif params.condition_toggle == 3
        labels = {'C-S', 'C-T','S-T'};
    end
end

close all; % close all figures currently open.

for thisROIcomp = 1:size(cond.ROIcomps,2) % for each ROI in turn:
    
    [condcount,plotcount] = deal(1); % intialise counter variables.
    figure('units','normalized','outerposition',[0 0 1 1])
    
    % for each pairwise condition comparison in turn:
    for thiscond1 = 1:size(params.condlabels,2)
        for thiscond2 = 1:size(params.condlabels,2)
            
            % for each pairwise comparison of conditions, if the conditions are
            % identical, skip to the next loop iteration.
            if thiscond1 == thiscond2, continue, end
            
            for thisiteration = 1:params.iterations % for a specified number of iterations:
                
                % extract the fisher-transformed correlation coefficients
                % corresponding to the pairwise condition data for the
                % desired ROI.
                if isequal(cond.ROIcomps{thisROIcomp},'V1')
                    cond.data(1,:,:) = squeeze(indcorr.fisherdata(thiscond1,:,2:end,1));
                    cond.data(2,:,:) = squeeze(indcorr.fisherdata(thiscond2,:,2:end,1));
                elseif isequal(cond.ROIcomps{thisROIcomp}, 'IPS0')
                    cond.data(1,:,:) = squeeze(indcorr.fisherdata(thiscond1,:,6,1:5));
                    cond.data(2,:,:) = squeeze(indcorr.fisherdata(thiscond2,:,6,1:5));
                elseif isequal(cond.ROIcomps{thisROIcomp}, 'V4')
                    cond.data(1,:,:) = squeeze(indcorr.fisherdata(thiscond1,:,3,[1,2,4:6]));
                    cond.data(2,:,:) = squeeze(indcorr.fisherdata(thiscond2,:,3,[1,2,4:6]));
                elseif isequal(cond.ROIcomps{thisROIcomp}, 'LO1')
                    cond.data(1,:,:) = squeeze(indcorr.fisherdata(thiscond1,:,4,[1,2,3,5,6]));
                    cond.data(2,:,:) = squeeze(indcorr.fisherdata(thiscond2,:,4,[1,2,3,5,6]));
                elseif isequal(cond.ROIcomps{thisROIcomp}, 'LO2')
                    cond.data(1,:,:) = squeeze(indcorr.fisherdata(thiscond1,:,5,[1,2,3,4,6]));
                    cond.data(2,:,:) = squeeze(indcorr.fisherdata(thiscond2,:,5,[1,2,3,4,6]));
                elseif isequal(cond.ROIcomps{thisROIcomp},'V3AB')
                    cond.data(1,:,:) = squeeze(indcorr.fisherdata(thiscond1,:,2,[1,3:6]));
                    cond.data(2,:,:) = squeeze(indcorr.fisherdata(thiscond2,:,2,[1,3:6]));
                end
                
                % create vector representng individual participants for random
                % selection.
                cond.idx = 1:size(original.ccMat,2);
                
                % randomly sample from participants with replacement.
                cond.sample = datasample(cond.idx,length(pps), 'Replace', true);
                
                % extract condition-specific randomly-sampled data.
                cond.datasample(1,:,:) = squeeze(cond.data(1,cond.sample,:));
                cond.datasample(2,:,:) = squeeze(cond.data(2,cond.sample,:));
                
                for thispp = 1:size(cond.datasample,2) % for each randomly-sampled data iteration:
                    
                    cond.ppdata = [squeeze(cond.datasample(1,thispp,:)),...
                        squeeze(cond.datasample(2,thispp,:))];
                    
                    for thiscorr = 1:size(cond.ppdata,1) % for each ROI-pariwse correlation coefficient:
                        
                        r = rand(); % generate a random number between 0 and 1.
                        
                        % flip the position of the condition 1 and condition 2
                        % correlation coefficients 50% of the time to create a
                        % noise (scrambled) dataset.
                        if r < .5
                            cond.noisedata(1,thispp,thiscorr) = cond.datasample(1,thispp,thiscorr);
                            cond.noisedata(2,thispp,thiscorr) = cond.datasample(2,thispp,thiscorr);
                        else
                            cond.noisedata(1,thispp,thiscorr) = cond.datasample(2,thispp,thiscorr);
                            cond.noisedata(2,thispp,thiscorr) = cond.datasample(1,thispp,thiscorr);
                        end
                    end
                end
                
                % calculate mean across participants for each condition (for noise
                % and original datasets).
                cond.datamean(1,thisiteration,:) = mean(squeeze(cond.datasample(1,:,:)));
                cond.datamean(2,thisiteration,:) = mean(squeeze(cond.datasample(2,:,:)));
                
                cond.noisemean(1,thisiteration,:) = mean(squeeze(cond.noisedata(1,:,:)));
                cond.noisemean(2,thisiteration,:) = mean(squeeze(cond.noisedata(2,:,:)));
                
                % calculate the iteration root mean squared error (for noise
                % and original datasets).
                cond.datarmse(thisiteration) = sum((squeeze(cond.datamean(1,thisiteration,:))'-...
                    squeeze(cond.datamean(2,thisiteration,:))').^2);
                cond.noisermse(thisiteration) = sum((squeeze(cond.noisemean(1,thisiteration,:))'-...
                    squeeze(cond.noisemean(2,thisiteration,:))').^2);
            end
            
            % calcuate the root mean squared error for the noise distribution
            % across all iterations.
            cond.noisermsemean = mean(cond.noisermse);
            
            % calculate the percentage of original RMSEs falling below the
            % mean RMSE of the noise distribution.
            cond.percentile(condcount) = sum(cond.datarmse<cond.noisermsemean)/size(cond.noisermse,2);
            
            % extract the outline of the noise distribution for later
            % plotting.
            [f,xi] = ksdensity(cond.noisermse(:),'NumPoints', 50);
            h = histcounts(cond.noisermse(:),50,'Normalization','probability');
            mult = max(h)/max(f);
            
            % if this comparison is a comparison we specified we were
            % interested in:
            if ismember(condcount, comps)
                
                % extract percentile data for alternative plotting.
                comp_percentile(plotcount) = cond.percentile(condcount);
                
                % create a histogram with the distribution of our observed
                % euclidean distances and the outline of the distribution
                % of noise distances.
                if isequal(params.cond, 'original') || isequal(params.cond, 'feature') || isequal(params.cond, 'colour') || isequal(params.cond, 'colourxfeature') && condition_toggle > 3
                    subplot(2, 3, plotcount);
                elseif isequal(params.cond, 'colourxfeature') && condition_toggle < 4 || isequal(params.cond, '3x3feature')
                    subplot(1, 3, plotcount);
                elseif isequal(params.cond, '3feature')
                    subplot(2, 5, plotcount);
                end
                f = f*mult;
                
                histogram(cond.datarmse(:),50, 'FaceColor', [ 0.5843 0.8157 0.9882],...
                    'EdgeColor', [ 0.5843 0.8157 0.9882],'Normalization','probability');
                
                % scale the x-axis in a condition-specific fashion.
                if isequal(params.cond, 'feature')
                    if plotcount == 1
                        set(gca,'XLim',[0 .02]);
                    elseif plotcount == 2
                        set(gca,'XLim',[0 .04]);
                    elseif plotcount == 3
                        set(gca,'XLim',[0 .06]);
                    elseif plotcount == 4
                        set(gca,'XLim',[0 .04]);
                    elseif plotcount == 5
                        set(gca,'XLim',[0 .06]);
                    elseif plotcount == 6
                        set(gca,'XLim', [0 .08]);
                    end
                elseif isequal(params.cond, 'colour')
                    if plotcount == 1
                        set(gca,'XLim',[0 .5]);
                    elseif plotcount == 2
                        set(gca,'XLim',[0 .05]);
                    elseif plotcount == 3
                        set(gca,'XLim',[0 .05]);
                    elseif plotcount == 4
                        set(gca,'XLim',[0 .5]);
                    elseif plotcount == 5
                        set(gca,'XLim',[0 .05]);
                    elseif plotcount == 6
                        set(gca,'XLim', [0 .5]);
                    end
                elseif isequal(params.cond, '3feature')
                    if ismember(plotcount, [1 2 3 7 9 10])
                        set(gca,'XLim',[0 .4]);
                    elseif ismember(plotcount, [5 8])
                        set(gca,'XLim',[0 .1]);
                    elseif plotcount == 4
                        set(gca,'XLim',[0 .06]);
                    elseif plotcount == 6
                        set(gca,'XLim',[0 .1]);
                    end
                elseif isequal(params.cond, 'colourxfeature')
                    if params.condition_toggle == 1
                        if plotcount == 1
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 2
                            set(gca,'XLim',[0 .4]);
                        elseif plotcount == 3
                            set(gca,'XLim',[0 .5]);
                        end
                    elseif params.condition_toggle == 2
                        if plotcount == 1
                            set(gca,'XLim',[0 .4]);
                        elseif plotcount == 2
                            set(gca,'XLim',[0 .6]);
                        elseif plotcount == 3
                            set(gca,'XLim',[0 .4]);
                        end
                    elseif params.condition_toggle == 3
                        if plotcount == 1
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 2
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 3
                            set(gca,'XLim',[0 .3]);
                        end
                    elseif params.condition_toggle == 4
                        set(gca,'XLim', [0 .4]);
                    elseif params.condition_toggle == 5
                        if plotcount == 1
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 2
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 3
                            set(gca,'XLim',[0 .4]);
                        elseif plotcount == 4
                            set(gca,'XLim',[0 .4]);
                        elseif plotcount == 5
                            set(gca,'XLim',[0 .5]);
                        elseif plotcount == 6
                            set(gca,'XLim',[0 .6]);
                        end
                    elseif params.condition_toggle == 6
                        if plotcount == 1
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 2
                            set(gca,'XLim',[0 .3]);
                        elseif plotcount == 3
                            set(gca,'XLim',[0 .4]);
                        elseif plotcount == 4
                            set(gca,'XLim',[0 .2]);
                        elseif plotcount == 5
                            set(gca,'XLim',[0 .2]);
                        elseif plotcount == 6
                            set(gca,'XLim',[0 .2]);
                        end
                    end
                elseif isequal(params.cond, '3x3feature')
                    if params.condition_toggle == 1 || params.condition_toggle == 2
                        if plotcount == 1
                            set(gca, 'XLim', [0 .2]);
                        elseif plotcount == 2
                            set(gca, 'XLim', [0 .2]);
                        elseif plotcount == 3
                            set(gca, 'XLim', [0 .2]);
                        end
                    elseif params.condition_toggle == 3
                        if plotcount == 1
                            set(gca, 'XLim', [0 .2]);
                        elseif plotcount == 2
                            set(gca, 'XLim', [0 .15]);
                        elseif plotcount == 3
                            set(gca, 'XLim', [0 .2]);
                        end
                    end
                end
                
                % specify further formatting variables.
                set(gca,'YLim',[0 .2]);
                set(gca, 'FontSize', 13);
                hold on
                h=plot(xi,f);
                h(1).LineWidth = 1.5;
                h(1).Color = 'k';
                h(1).LineStyle = '- -';
                xlabel(labels{plotcount})
                hold off
                plotcount = plotcount + 1;
            end
            
            condcount = condcount + 1; % increment the ROI pairwise comparison counter variable.
        end
        
        % save all ROI percentile data for alternative plotting.
        data_percentile{thisROIcomp} = comp_percentile;
    end
    
    fig = gcf;
    
    % save the figure and ROI-specific analysis.
    print(fig, '-dpdf', strcat(params.outputdir,...
        sprintf('histogram_%s_condition_correlations_probaility-fisher.pdf', cond.ROIcomps{thisROIcomp})) ,'-bestfit');
    save(strcat(params.outputdir, ...
        sprintf('%s_condition_correlations_data-fisher.mat', cond.ROIcomps{thisROIcomp})),'cond');
    close all
end

%% -------------- ALTERNATIVE EUCLIDEAN DISTANCE PLOTTING -------------- %%

% combine euclidean distance data across ROIs.
plotdata = [];
for thisROI = 1:length(cond.ROIcomps)
    plotdata = [plotdata; data_percentile{thisROI}];
end

% convert overlap values to percentages.
plotdata_percentage = plotdata.*100;

% specify barchart colours.
C = [.15 .3 .45 .6 .75 .9 ;.15 .3 .45 .6 .75 .9 ;.15 .3 .45 .6 .75 .9]';

% plot percentage data as a bar chart with associated significance values.
superbar(plotdata_percentage, 'P', plotdata, 'BarFaceColor', C, 'BarEdgeColor',...
    [0 0 0], 'BarLineWidth', 1, 'PStarColor', [0 0 0],'PStarFontSize', 12,...
    'PLineColor', [0 0 0], 'PLineWidth', 1.2, 'PStarShowNS', false,...
    'PStarShowGT', false, 'PStarThreshold', [0.05, 0.01, 0.001]);
ylim([0 30]); ylabel('Percentage Overlap');
set(gca, 'xtick', 1:6, 'xticklabel', cond.ROIcomps);
hline = refline([0 5]);
fig = gcf;

% save figure to a .pdf file.
print(fig, '-dpdf', strcat(params.outputdir,...
    sprintf('barchart_%s_condition_correlations_probaility-fisher.pdf', params.cond)) ,'-bestfit');

%% ------------------------ SAVE ANALYSIS DATA ------------------------- %%

% load in ROI-specific euclidean distance values.
V1cond = load(strcat(params.outputdir, 'V1_condition_correlations_data-fisher.mat'));
V3ABcond = load(strcat(params.outputdir, 'V3AB_condition_correlations_data-fisher.mat'));
V4cond = load(strcat(params.outputdir, 'V4_condition_correlations_data-fisher.mat'));
LO1cond = load(strcat(params.outputdir, 'LO1_condition_correlations_data-fisher.mat'));
LO2cond = load(strcat(params.outputdir, 'LO2_condition_correlations_data-fisher.mat'));
IPS0cond = load(strcat(params.outputdir, 'IPS0_condition_correlations_data-fisher.mat'));

% save all data to a .mat file.
connectivity.rawcorrelations = original;
connectivity.indconditioncorrs = indcorr;
connectivity.rmsedistance.V1 = V1cond;
connectivity.rmsedistance.V3AB = V3ABcond;
connectivity.rmsedistance.V4 = V4cond;
connectivity.rmsedistance.LO1 = LO1cond;
connectivity.rmsedistance.LO2 = LO2cond;
connectivity.rmsedistance.IPS0 = IPS0cond;

if isequal(params.expcond, 'original')
    connectivity.thresholdcorrs = thresholds;
end

save('/Users/kirstie/Desktop/connectivity_colour-colour-fisher.mat', 'connectivity');

%% --------------------------------------------------------------------- %%