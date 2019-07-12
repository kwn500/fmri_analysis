%% --------------- 2-WAY SUPPORT VECTOR MACHINE DECODING --------------- %%

% classifies pairwise attentional focus from top 100 voxels' from each ROI
% and assesses if classification performance is significantly greater than
% chance.

% KWN 22/01/2019

clear; close all;

% add function directories.
addpath(genpath('/groups/labs/wadelab/toolbox/libsvm/libsvm-3.22/matlab'));
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions');

%% ---------------------- PARAMETERS TO EDIT --------------------------- %%

output.expcond = 'colour'; % project (original, colour, naturalistic).
output.cond = 'colour'; % analysis (block, colour, feature, na).
cond.facetoggle = 1; % 1 yes, 0 no.

% statistical analysis type (parametric, nonparametric).
output.analysis = 'nonparametric';

output.ROIs = {'V1' 'V2', 'V3', 'V4'}; % region-of-interest names.
output.voxels = 100; % number of voxels retained per ROI.

if isequal(output.expcond, 'original') || isequal(output.cond, 'feature') || isequal(output.cond,'feature_event')
    
    % data columns to compare, and comparison-specific names.
    output.comp = {[1,2], [1,3], [2,3]};
    output.compnames = {'OvsC', 'OvsS', 'CvsS'};
    
    % analysis-specific partcipants and data root directory.
    if isequal(output.expcond, 'original')
        participants = {2268 2548 2590 2904 3111 3455 3517 3773 3932 4065 4127 4496}';
        numparticipants = 12;
        
        if isequal(output.cond, 'block')
            output.rootdirectory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/';
        elseif isequal(output.cond, 'event')
            output.rootdirectory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/event/';
        end
    elseif isequal(output.cond, 'feature')
        participants = {2268 2548 2590 2904 3111 3517 3773 4059 4065 4244 4831 4928}';
        numparticipants = 12;
        
        output.rootdirectory = '/scratch/groups/Projects/P1323/colour/fMRI/vista_output/multivariate/feature/';
    elseif isequal(output.cond, 'feature_event')
        participants = {2268 2548 2590 2904 3111 3517 3773 4059 4065 4244 4831 4928}';
        numparticipants = 12;
        
        output.rootdirectory = '/scratch/groups/Projects/P1323/colour/fMRI/vista_output/multivariate/feature_event/';
    end
    
    % condition-specific output figure filename.
    output.graphfilename = sprintf('%s_OCS_twoway-%s', output.cond, output.analysis);
    
elseif isequal(output.expcond, 'colour')
    
    if isequal(output.cond, 'colour')
        output.comp = {[1,2], [1,3], [2,3]};
        output.compnames =  {'RGvsBY', 'RGvsLUM', 'BYvsLUM'};
        
        output.rootdirectory = '/scratch/groups/Projects/P1323/colour/fMRI/vista_output/multivariate/colour/';
        output.graphfilename = sprintf('%s_RGBYLUM_twoway-%s', output.cond, output.analysis);
        
    elseif isequal(output.cond, 'colour-passive')
        output.comp = {[1,3], [1,5], [3,5]};
        output.compnames = {'RGvsBY', 'RGvsLUM', 'BYvsLUM'};
        
        output.rootdirectory = '/scratch/groups/Projects/P1323/colour/fMRI/vista_output/multivariate/colour-passive/';
        output.graphfilename = sprintf('%s_RGBYLUM_passive_twoway-%s', output.cond, output.analysis);
        
    elseif isequal(output.cond, 'colourxfeature')
        if output.analysis_toggle == 1
            output.comp = {[1,2], [1,3], [2,3]};
            output.compnames = {'RGO-RGC', 'RGO-RGS', 'RGC-RGS'};
            output.graphfilename = sprintf('RedGreen-Feature-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 2
            output.comp = {[5,6], [5,7], [6,7]};
            output.compnames = {'BYO-BYC', 'BYO-BYS', 'BYC-BYS'};
            output.graphfilename = sprintf('BlueYellow-Feature-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 3
            output.comp = {[9,10], [9,11], [10,11]};
            output.compnames = {'LUMO-LUMC', 'LUMO-LUMS', 'LUMC-LUMS'};
            output.graphfilename = sprintf('Luminance-Feature-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 4
            output.comp = {[1,5], [1,9], [5,9]};
            output.compnames = {'ORG-OBY', 'ORG-OLUM', 'OBY-OLUM'};
            output.graphfilename = sprintf('Orientation-Colour-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 5
            output.comp =  {[2,6], [2,10], [6,10]};
            output.compnames = {'CRG-CBY', 'CRG-CLUM', 'CBY-CLUM'};
            output.graphfilename = sprintf('Contrast-Colour-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 6
            output.comp = {[3,7], [3,11], [7,11]};
            output.compnames = {'SRG-SBY', 'SRG-SLUM', 'SBY-SLUM'};
            output.graphfilename = sprintf('Shape-Colour-twoway-%s', output.analysis);
        end
        
        output.rootdirectory = '/scratch/groups/Projects/P1323/colour/fMRI/vista_output/multivariate/colourxfeature/';
        participants = {2268 2548 2590 2904 3111 3517 3773 4059 4065 4244 4831 4928}';
        numparticipants = 12;
    end
    
elseif isequal(output.expcond, 'naturalistic')
    
    if isequal(output.cond, '3feature')
        if isequal(cond.facetoggle,1)
            output.comp =  {[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]};
            output.compnames = {'OvsC', 'OvsS', 'OvsF', 'CvsS', 'CvsF', 'SvsF'};
        elseif isequal(cond.facetoggle,0)
            output.comp =  {[1,2], [1,3], [2,3]};
            output.compnames = {'OvsC', 'OvsS', 'CvsS'};
        end
        
        participants = {2548, 2904, 3111, 3517, 4059, 4065, 4127, 4244, 4831, 4833, 4890, 4928,}';
        numparticipants = 12;
        
        output.rootdirectory = '/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3feature/';
        output.graphfilename = sprintf('%s_OCSF_twoway-%s', output.expcond, output.analysis);
    elseif isequal(output.cond, '3x3feature')
        if output.analysis_toggle == 1
            output.comp = {[2,3], [2,4], [3,4]};
            output.compnames =  {'V-H', 'V-D', 'H-D'};
            output.graphfilename = sprintf('orientation3x3-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 2
            output.comp = {[5,6], [5,7], [6,7]};
            output.compnames = {'R-G', 'R-B', 'G-B'};
            output.graphfilename = sprintf('colour3x3-twoway-%s', output.analysis);
        elseif output.analysis_toggle == 3
            output.comp = {[8,9], [8,10], [9,10]};
            output.compnames = {'C-S', 'C-T', 'S-T'};
            output.graphfilename = sprintf('shape3x3-twoway-%s', output.analysis);
        end
        
        participants = {2548, 2904, 3111, 3517, 3773, 4059, 4065, 4127, 4244, 4831, 4833, 4890, 4928, 5006};
        numparticipants = 14;
    end
end

output.datadirectory = [output.rootdirectory,'betasxvarexp/']; % beta storage directory.

% analysis and figure output directories.
output.directory = [output.rootdirectory, 'backprojection/libsvmclassification/twoway/'];
[~,~] = mkdir(output.directory);
output.figuredirectory = [output.rootdirectory, 'backprojection/libsvmclassification/twoway/figures/'];
[~,~] = mkdir(output.figuredirectory);

%% -------------------- PERFORM DECODING ANALYSIS ---------------------- %%

for roinum = 1:length(output.ROIs) % for each ROI:
    
    ROI = output.ROIs{roinum}; % extract ROI name.
    
    % load beta values.
    if isequal(output.expcond, 'naturalistic')
        if isequal(ROI, 'OFA') || isequal(ROI, 'STS')
            load(strcat(output.datadirectory, sprintf('%s/%d/group_extractedBetas_%s.mat', ROI, output.voxels, ROI)));
        else
            load(strcat(output.datadirectory, sprintf('%s/%d/group_extractedbetas_%s.mat', ROI, output.voxels, ROI)));
        end
    else
        load(strcat(output.datadirectory, sprintf('%s/%d/group_extractedBetas_%s.mat', ROI, output.voxels, ROI)));
    end
    
    for participant = 1:length(group_betas) % for each participant:
        
        if isequal(output.expcond, 'original') && isequal(output.cond, 'block') ||...
                isequal(output.expcond, 'colour') && ~isequal(output.cond,'feature_event') ||...
                isequal(output.expcond, 'colour') && ~isequal(output.cond, 'colour-passive') % for RFP experiment data:
            
            output.origdata = group_betas{participant}; % extract participant beta data.
            
            for compnum = 1:length(output.comp) % for each pairwise comparison:
                
                comp = output.comp{compnum}; % extract pair columns numbers.
                
                output.targetdata = output.origdata(:,comp,:); % extract corresponding condition data.
                
                % record repetition, condtion and voxel numbers.
                [output.repnums, output.condnums, output.voxelnums] = size(output.targetdata);
                
                % extract condition-specific data individually.
                output.targetcond1 = squeeze((output.targetdata(:,1,:)));
                output.targetcond2 = squeeze((output.targetdata(:,2,:)));
                
                % vertically stack the two conditions data.
                output.targetreshape = [output.targetcond1; output.targetcond2];
                
                % z-score across voxels.
                output.targetzscore = zscore(output.targetreshape,[],2);
                
                % create column of labels matching order of conditions.
                output.labels = [ones(output.repnums,1); ones(output.repnums,1)*2];
                
                % run decoding with radial basis function kernel and
                % leave-one-out cross validation (number of reps -1).
                output.params = sprintf('-t 2, -v %d -q', size(output.targetcond1,1)-1 );
                acc = libsvmtrain(output.labels, output.targetzscore, output.params);
                
                compacc(compnum) = acc; % store pairwise comparison-specific data.
            end % next pairwise comparison.
            
        elseif isequal(output.expcond, 'naturalistic') || isequal(output.cond, 'event') ||...
                isequal(output.cond, 'feature_event') || isequal(output.cond, 'colour-passive') % for naturalistic experiment:
            
            % extract participant beta data.
            output.origdata = group_betas{participant};
            
            for compnum = 1:length(output.comp) % for each pairwise comparison:
                
                data1 = squeeze(output.origdata(:,output.comp{compnum}(1),:));
                data2 = squeeze(output.origdata(:,output.comp{compnum}(2),:));
                
                data1(~any(~isnan(data1), 2),:)=[];
                data2(~any(~isnan(data2), 2),:)=[];
                
                output.targetdata{1} = data1;
                output.targetdata{2} = data2;
                
                output.targetreshape = [data1;data2];
                output.sizes = [size(data1,1), size(data2,1)];
                
                % repeat same process as above.
                output.targetzscore = zscore(output.targetreshape,[],2);
                output.labels = [ones(size(data1,1),1); ones(size(data2,1),1)*2];
                
                [maxsize, ind] = max(output.sizes);
                
                weights = 1:2;
                weights(ind) = [];
                
                output.weights{ind} = 1;
                output.weights{weights} = maxsize/size(output.targetdata{weights},1);
                
                % perform weighted support vectore machine classification.
                output.params = sprintf('-t 2 -v %d -w1 %.2f -w2 %.2f -q',...
                    size(output.targetdata{ind},1)-1, output.weights{1}, output.weights{2});
                acc = libsvmtrain(output.labels, output.targetzscore, output.params);
                
                compacc(compnum) = acc; % store pairwise comparison-specific data.
            end % next pairwise comparison.
        end
        groupacc(participant,:) = compacc; % store participant-specific data.
    end % next participant.
    
    [analysis.means, analysis.stderrs, analysis.stats, analysis.sig] = deal([]);
    for i = 1:size(groupacc,2) % for each comparison:
        
        % calculate mean and standard error across participants.
        analysis.mean(i) = mean(groupacc(:,i));
        analysis.stderr(i) = std(groupacc(:,i))/sqrt(length(groupacc(:,i)));
        
        % perform parametric or nonparametric analysis versus chance.
        if isequal(output.analysis, 'parametric')
            [analysis.h(i), analysis.p(i), ~, analysis.fullstats(i)] = ttest(groupacc(:,i), (100/2));
            analysis.stat(i) = analysis.fullstats(i).tstat;
        elseif isequal(output.analysis, 'nonparametric')
            [analysis.p(i), analysis.h(i), analysis.fullstats(i)] = signrank(groupacc(:,i),(100/2), 'method', 'approximate');
            analysis.stat(i) = analysis.fullstats(i).zval;
        end
        
        % concatenate data across comparisons.
        analysis.means = [analysis.means, analysis.mean(i)];
        analysis.stderrs = [analysis.stderrs, analysis.stderr(i)];
        analysis.stats = [analysis.stats, analysis.stat(i)];
        analysis.sig = [analysis.sig, analysis.p(i)];
    end % next comparison.
    
    % store all comparison data for each ROI.
    roi.acc{roinum} = groupacc;
    roi.means{roinum} = analysis.means;
    roi.stderrs{roinum} = analysis.stderrs;
    roi.stats{roinum} = analysis.stats;
    roi.sig{roinum} = analysis.sig;
end % next ROI.

% benjamini-hochberg correct p values across all comparisons and ROIs.
[analysis.adjh, ~, ~, analysis.adjp]=fdr_bh(cell2mat(roi.sig),0.05,'dep','yes');

% store participant decoding accuracies, means, standard errors,
% test statistic, significance values and adjusted significance values.
outputtable = [cell2mat(roi.acc); cell2mat(roi.means); cell2mat(roi.stderrs); cell2mat(roi.stats); cell2mat(roi.sig); analysis.adjp];
rowheaders = [participants; 'mean'; 'stderr'; output.analysis;'sig'; 'adjsig']; % corresponding row headers.
outputtable = [rowheaders, num2cell(outputtable)];

% for each ROI, create column headers matching comparison names.
count = 1;
for i = 1:length(output.ROIs)
    for x = 1:length(output.compnames)
        colheaders{count} = strcat(output.ROIs{i}, output.compnames{x});
        count = count + 1;
    end
end
colheaders = ['-----------', colheaders];
outputtable = [colheaders; outputtable];

cell2csv(strcat(output.directory, output.graphfilename, '.csv'), outputtable); % write data to file.

count = 1;
for i = 1:size(roi.sig,2) % for each ROI:
    
    % plot data with error bars and significance values.
    superbar(cell2mat(roi.means(i)), 'E', cell2mat(roi.stderrs(i)), 'P', analysis.adjp(count:count+ size(roi.sig{1},2)-1),...
        'BarFaceColor', [.5 .5 .5], 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1.5,...
        'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0], 'PStarFontSize', 16,...
        'PStarShowGT', false, 'PStarShowNS', false);
    
    set(gca, 'xtick', 1:length(output.comp), 'xticklabel', output.compnames);
    ylim([30 100]);
    xlabel('Conditions of Interest');
    ylabel('Classification Accuracy (%)');
    
    hline = refline([0 50]); % add chance performance line.
    hline.Color = 'r'; hline.LineWidth = 1; hline.LineStyle = '--';
    
    % write figure to file.
    saveas(gcf, strcat(output.figuredirectory, output.graphfilename, sprintf('-%s-face.pdf', output.ROIs{i})));
    
    count = count + size(roi.sig{1},2); % increment ROI x comparison counter.
end % next ROI.

close all; % close all figures.

% this needs making condition-specific- combine roi-specific p-values
% across all conditions and experiments.
adj_reshape = [{analysis.adjp(1:3)},{analysis.adjp(4:6)},{analysis.adjp(7:9)},{analysis.adjp(10:12)}];

% save important analysis data to .mat file.
multivariate_block_twoway.data = roi.acc;
multivariate_block_twoway.means = roi.means;
multivariate_block_twoway.stderr = roi.stderrs;
multivariate_block_twoway.stats = roi.stats;
multivariate_block_twoway.p = roi.sig;
multivariate_block_twoway.adj_p = adj_reshape;

save(strcat(output.rootdirectory, output.graphfilename, '-twoway.mat'), 'multivariate_block_twoway');

%% --------------------------------------------------------------------- %%