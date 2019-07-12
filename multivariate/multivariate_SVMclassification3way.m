%% --------------- 3-WAY SUPPORT VECTOR MACHINE DECODING --------------- %%

% Attempts to classify attentional-focus from top 100 voxels' for each ROI,
% then assess if significantly greater than chance.

% KWN 22/01/2019

clear; close all;

% add function directories.
addpath(genpath('/groups/labs/wadelab/toolbox/libsvm/libsvm-3.22/matlab'));
addpath('/scratch/groups/Projects/P1323/code/fmri_analysis/functions');

%% ---------------------- PARAMETERS TO EDIT --------------------------- %%

output.expcond = 'original'; % project (original, colour, naturalistic).
output.cond = 'block'; % analysis (block, colour, feature, na).

% statistical analysis type (parametric, nonparametric).
output.analysis = 'nonparametric';
output.analysissize = 3; % threeway or fourway (naturalistic).

% specify colourxfeature interaction to analyse (1-6).
output.analysis_toggle = 1;

output.ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2', 'A1'}; % ROI names.
output.voxels = 100; % number of voxels retained per ROI.

% analysis-specific participants.
if isequal(output.expcond, 'original')
    participants = {2268 2548 2590 2904 3111 3455 3517 3773 3932 4065 4127 4496}';
elseif isequal(output.expcond, 'colour')
    participants = {2268 2548 2590 2904 3111 3517 3773 4059 4065 4244 4831 4928}';
elseif isequal(output.expcond, 'naturalistic')
    participants = {2548, 2904, 3111, 3517, 4059, 4065, 4127, 4244, 4831, 4833, 4890, 4928}';
end

% specify some default values for looping when we're not processing
% colourxfeature interaction data.
comparison = {1};
comparison_name = {'NA'};

% condition-specific output figure filename.
if isequal(output.expcond, 'original') || isequal(output.cond, 'feature') || isequal(output.cond, 'feature_event')
    output.graphfilename = sprintf('%s_OCS_threeway-%s', output.cond, output.analysis);
elseif isequal(output.cond, 'colour')
    output.graphfilename = sprintf('%s_RGBYLUM_threeway-%s', output.cond, output.analysis);
elseif isequal(output.cond, 'colour-passive')
    output.graphfilename = sprintf('%s_RGBYLUM_passive_threeway-%s', output.cond, output.analysis);
elseif isequal(outputl.cond, 'colourxfeature')
    if output.analysis_toggle == 1
        comparison = [1 2 3];
        comparison_name = 'RedGreen-Feature';
    elseif output.analysis_toggle == 2
        comparison = [5 6 7];
        comparison_name = 'BlueYellow-Feature';
    elseif output.analysis_toggle == 3
        comparison = [9 10 11];
        comparison_name = 'Luminance-Feature';
    elseif output.analysis_toggle == 4
        comparison = [1 5 9];
        comparison_name = 'Orientation-Colour';
    elseif output.analysis_toggle == 5
        comparison = [2 6 10];
        comparison_name = 'Contrast-Colour';
    elseif output.analysis_toggle == 6
        comparison = [3 7 11];
        comparison_name = 'Shape-Colour';
    end
    output.graphfilename  = sprintf('%s-%s', comparison_name, output.analysis);
elseif isequal(output.expcond, 'naturalistic')
    if isequal(output.cond, '3feature')
        if output.analysissize == 4
            output.graphfilename = sprintf('%s_OCSF_fourway-%s', output.expcond, output.analysis);
        elseif output.analysissize == 3
            output.graphfilename = sprintf('%s_OCSF_threeway-%s', output.expcond, output.analysis);
        end
    elseif isequal(output.cond, '3x3feature')
        if output.analysis_toggle == 1
            comparison = [2 3 4];
            comparison_name = 'orientation3x3';
        elseif output.analysis_toggle == 2
            comparison = [5 6 7];
            comparison_name = 'colour3x3';
        elseif output.analysis_toggle == 3
            comparison = [8 9 10];
            comparison_name = 'shape3x3';
        end
        output.graphfilename  = sprintf('%s-%s', comparison_name, output.analysis);
    end
end

% condition-specific data root directory.
if isequal(output.expcond, 'original') && isequal(output.cond, 'block') || isequal(output.expcond, 'colour') || isequal(output.cond, 'event')
    output.rootdirectory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/', output.expcond, output.cond);
elseif isequal(output.expcond, 'original') && isequal(output.cond, 'event')
    output.rootdirectory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/event/';
elseif isequal(output.cond, 'feature_event')
    output.rootdirectory = '/scratch/groups/Projects/P1323/colour/fMRI/vista_output/multivariate/feature_event/';
elseif isequal(output.expcond, 'naturalistic')
    output.rootdirectory = sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3feature/');
end

% beta, analysis and figure output directories.
output.datadirectory = [output.rootdirectory,'betasxvarexp/']; % beta storage directory.
if output.analysissize == 3
    output.directory = [output.rootdirectory, 'libsvmclassification/threeway/']; [~,~] = mkdir(output.directory);
    output.figuredirectory = [output.rootdirectory, 'libsvmclassification/threeway/figures/'];  [~,~] = mkdir(output.figuredirectory);
elseif output.analysissize == 4
    output.directory = [output.rootdirectory, 'libsvmclassification/fourway/']; [~,~] = mkdir(output.directory);
    output.figuredirectory = [output.rootdirectory, 'libsvmclassification/fourway/figures/'];  [~,~] = mkdir(output.figuredirectory);
end

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
        data = group_betas;
    end
    
    for participant = 1:length(group_betas) % for each participant:
        
        if ~isequal(output.expcond, 'naturalistic') && isequal(output.cond, 'block') || isequal(output.expcond, 'colour') && ~isequal(output.cond,'feature_event') % for RFP experiment data:
            
            % extract participant beta data.
            output.origdata = group_betas{participant};
            
            if isequal(output.cond, 'colour-passive')
                
                % extract corresponding condition data.
                output.targetdata = output.origdata(:,[1 3 5],:);
                
            elseif isequal(output.cond, 'colourxfeature')
                output.targetdata = output.origdata(:,comparison,:);
                
            else
                % extract corresponding condition data.
                output.targetdata = output.origdata(:,[1 2 3],:);
            end
            
            % record repetition, condtion and voxel numbers.
            [output.repnums, output.condnums, output.voxelnums] = size(output.targetdata);
            
            % extract condition-specific data individually.
            output.targetcond1 = squeeze((output.targetdata(:,1,:)));
            output.targetcond2 = squeeze((output.targetdata(:,2,:)));
            output.targetcond3 = squeeze((output.targetdata(:,3,:)));
            
            % vertically stack the three conditions data.
            output.targetreshape = [output.targetcond1; output.targetcond2; output.targetcond3];
            
            % z-score across voxels.
            output.targetzscore = zscore(output.targetreshape,[],2);
            
            % create column of labels matching order of conditions.
            output.labels = [ones(output.repnums,1); ones(output.repnums,1)*2; ones(output.repnums, 1)*3];
            
            % run decoding with radial basis function kernel and
            % leave-one-out cross validation (number of reps -1).
            output.params = sprintf('-t 2, -v %d -q', size(output.targetcond1,1)-1 );
            acc = libsvmtrain(output.labels, output.targetzscore, output.params);
            
            groupacc(participant,:) = acc; % store participant-specific data.
            
            % else, if we are analysing naturalistic or event-related data:
        elseif isequal(output.expcond, 'naturalistic') || isequal(output.cond, 'event') || isequal(output.cond,'feature_event')
            
            % perform the same process as above and extract the
            % condition-specific data, rehaping the data in accordance with
            % the three-way or four-way analysis we have specified.
            output.origdata = group_betas{participant};
            
            if ~isequal(output.cond, '3x3feature')
                for conddata = 1:output.analysissize
                    output.targetdata{conddata} = squeeze(output.origdata(:,conddata,:));
                    output.targetdata{conddata}(~any(~isnan(output.targetdata{conddata}), 2),:)=[];
                end
            elseif isequal(output.cond, '3x3feature')
                for cond =  1:size(comparison_number,2)
                    output.targetdata{cond} = squeeze(output.origdata(:,comparison_number(cond),:));
                    output.targetdata{cond}(~any(~isnan(output.targetdata{cond}),2),:) = [];
                end
            end
            
            if output.analysissize == 4
                output.targetreshape = [output.targetdata{1}; output.targetdata{2};...
                    output.targetdata{3}; output.targetdata{4}];
                output.sizes = [size(output.targetdata{1},1), size(output.targetdata{2},1),...
                    size(output.targetdata{3},1), size(output.targetdata{4},1)];
            elseif output.analysissize == 3
                output.targetreshape = [output.targetdata{1}; output.targetdata{2};...
                    output.targetdata{3}];
                output.sizes = [size(output.targetdata{1},1), size(output.targetdata{2},1),...
                    size(output.targetdata{3},1)];
            end
            
            % calculate the maximum condition size, and specify the
            % remaining condition weights accordingly.
            [maxsize, ind] = max(output.sizes);
            
            weights = 1:output.analysissize;
            weights(ind) = [];
            output.weights{ind} = 1;
            
            for weightind = weights
                output.weights{weightind} = maxsize/size(output.targetdata{weightind},1);
            end
            
            % z-score across voxels.
            output.targetzscore = zscore(output.targetreshape,[],2);
            
            % perform weighted svm classification.
            if output.analysissize == 4
                % create column of labels matching order of conditions.
                output.labels = [ones(size(output.targetdata{1},1),1); ones(size(output.targetdata{2},1),1)*2;...
                    ones(size(output.targetdata{3},1),1)*3; ones(size(output.targetdata{4},1),1)*4];
                
                % run decoding with radial basis function kernel and leave-one-out cross validation (number of reps -1).
                output.params = sprintf('-t 2 -v %d -w1 %.2f -w2 %.2f -w3 %.2f -w4 %.2f -q',...
                    size(output.targetdata{ind},1)-1, output.weights{1}, output.weights{2}, output.weights{3}, output.weights{4});
                
            elseif output.analysissize == 3
                % create column of labels matching order of conditions.
                output.labels = [ones(size(output.targetdata{1},1),1); ones(size(output.targetdata{2},1),1)*2;...
                    ones(size(output.targetdata{3},1),1)*3];
                
                % run decoding with radial basis function kernel and
                % leave-one-out cross validation (number of reps -1).
                output.params = sprintf('-t 2 -v %d -w1 %.2f -w2 %.2f -w3 %.2f -q',...
                    size(output.targetdata{ind},1)-1, output.weights{1}, output.weights{2}, output.weights{3});
            end
            acc = libsvmtrain(output.labels, output.targetzscore, output.params);
            
            groupacc(participant,:) = acc; % store participant-specific data.
        end
    end % next participant.
    
    % calculate mean and standard error across participants.
    analysis.mean = mean(groupacc);
    analysis.stderr = std(groupacc)/sqrt(length(groupacc));
    
    % perform parametric or nonparametric analysis versus chance.
    if isequal(output.analysis, 'parametric')
        [analysis.h, analysis.p, ~, analysis.fullstats] = ttest(groupacc,(100/size(output.targetdata,2)));
        analysis.stat = analysis.fullstats.tstat;
    elseif isequal(output.analysis, 'nonparametric')
        [analysis.p, analysis.h, analysis.fullstats] = signrank(groupacc,(100/size(output.targetdata,2)), 'method', 'approximate');
        analysis.stat = analysis.fullstats.zval;
    end
    
    % store all comparison data for each ROI.
    roi.acc{roinum} = groupacc;
    roi.means{roinum} = analysis.mean;
    roi.stderrs{roinum} = analysis.stderr;
    roi.stats{roinum} = analysis.stat;
    roi.sig{roinum} = analysis.p;
end % next ROI.

% benjamini-hochberg correct p values across all comparisons and ROIs.
[analysis.adjh, ~, ~, analysis.adjp]=fdr_bh(cell2mat(roi.sig),0.05,'dep','yes');

% store participant decoding accuracies, means, standard errors,
% test statistic, significance values and adjusted significance values.
outputtable = [cell2mat(roi.acc); cell2mat(roi.means); cell2mat(roi.stderrs);...
    cell2mat(roi.stats); cell2mat(roi.sig); analysis.adjp];
rowheaders = [participants; 'mean'; 'stderr'; output.analysis;'sig'; 'adjsig']; % corresponding row headers.
outputtable = [rowheaders, num2cell(outputtable)];
colheaders = output.ROIs; colheaders = ['-----------', colheaders];  % corresponding column headers.
outputtable = [colheaders; outputtable];

cell2csv(strcat(output.directory, output.graphfilename, '.csv'), outputtable); % write data to file.

% plot data with error bars and significance values.
superbar(cell2mat(roi.means), 'E', cell2mat(roi.stderrs), 'P', analysis.adjp,...
    'BarFaceColor', [.5 .5 .5], 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1.5,...
    'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0], 'PStarFontSize', 16,...
    'PStarShowGT', false, 'PStarShowNS', false);

set(gca, 'xtick', 1:length(output.ROIs), 'xticklabel', output.ROIs);
ylim([20 80]);
xlabel('Regions of Interest');
ylabel('Classification Accuracy (%)');

hline = refline([0 100/size(output.targetdata,2)]); % add chance performance line.
hline.Color = 'r'; hline.LineWidth = 1; hline.LineStyle = '--';

% write figure to file.
saveas(gcf, strcat(output.figuredirectory, output.graphfilename, '.pdf'));

close all; % close all figures.

multivariate_block_threeway.data = roi.acc;
multivariate_block_threeway.means = roi.means;
multivariate_block_threeway.stderr = roi.stderrs;
multivariate_block_threeway.stats = roi.stats;
multivariate_block_threeway.p = roi.sig;
multivariate_block_threeway.adj_p = analysis.adjp;

save(strcat(output.rootdirectory, output.graphfilename), 'multivariate_block_threeway');

%% --------------------------------------------------------------------- %%