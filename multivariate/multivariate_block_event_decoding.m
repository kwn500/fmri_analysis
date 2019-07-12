%% -------------- BLOCK VERSUS EVENT MULTIVARIATE DECODING ------------- %%  

% performs three-way support vector machine classification, training the
% model on block data and training on event (or vice versa) to assess the
% generalisability of data from attention to bottom-up stimulus driven
% effects. 

% 10/7/19 KWN 

clear; clc;

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

ROIs = {'V1', 'V3AB', 'V4','LO1', 'LO2', 'A1'}; % specify ROIs to analyse.
voxels = 100; % specify number of voxels retained per ROI.

% specify participants consistent across the two conditions.
pps = [2268 2548 2590 2904 3111 3455 3517 3773 3932 4065 4127 4496];

% specify train-test order (1 = train attention blocks, 2 = train event). 
analysis = 2; 

% specify the condition labels for each set of condition data. 
attn.conds = [1,2,3];
event.conds = [1,2,3];

% specify data directories for both block and event data. 
attn.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/betasxvarexp/';
event.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/event/betasxvarexp/';

% specify output directory. 
output.directory = '/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block-event/libsvmclassification/';

%% ------------- RUN SUPPORT VECTOR MACHINE CLASSIFICATION ------------- %% 

for thisROI = 1:length(ROIs) % for each ROI in turn:
    ROI = ROIs{thisROI};
    
    for thispp = 1:length(pps) % for each participant in turn:
        pp = pps(thispp);
        
        % load in the block and event related betas. 
        load(strcat(attn.directory, sprintf('%s/100/', ROI), sprintf('R%d_extractedBetas_%s.mat', pp, ROI)));
        attn.data{thispp} = final_betas; clear final_betas;
        load(strcat(event.directory, sprintf('%s/100/', ROI), sprintf('R%d_extractedBetas_%s.mat', pp, ROI)));
        event.data{thispp} = final_betas; clear final_betas;

        % for each condition in the attention block data, remove any rows
        % containing nan values. 
        for attncond = 1:length(attn.conds)
            attncond_data{attncond} = squeeze(attn.data{thispp}(:,attncond,:));
            attncond_data{attncond}(isnan(attncond_data{attncond})) =[];
        end
        attn.data{thispp} = attncond_data; clear attncond_data
        
        % perform the same process for the event data. 
        for eventcond = 1:length(event.conds)
            eventcond_data{eventcond} = squeeze(event.data{thispp}(:,eventcond,:));
            eventcond_data{eventcond}(~any(~isnan(eventcond_data{eventcond}),2),:)=[];
        end
        event.data{thispp} = eventcond_data; clear eventcond_data
        
        % stack the data for each attention/feature condition vertically
        attn.stackdata{thispp} = [cell2mat(squeeze(attn.data{thispp}(:,1,:)));...
            cell2mat(squeeze(attn.data{thispp}(:,2,:)));cell2mat(squeeze(attn.data{thispp}(:,3,:)))];
        event.stackdata{thispp} = [cell2mat(squeeze(event.data{thispp}(:,1,:)));...
            cell2mat(squeeze(event.data{thispp}(:,2,:)));cell2mat(squeeze(event.data{thispp}(:,3,:)))];
        
        % z-score this data across voxels. 
        attn.zscoredata{thispp} = zscore(attn.stackdata{thispp},[],2);
        event.zscoredata{thispp} = zscore(event.stackdata{thispp},[],2);
        
        % create matrices with labels corresponding to the order of
        % conditions. 
        attn.labels{thispp} = [ones(1,size(attn.data{thispp}{1},1)),...
            ones(1,size(attn.data{thispp}{2},1))*2,ones(1,size(attn.data{thispp}{3},1))*3]';
        event.labels{thispp} = [ones(1,size(event.data{thispp}{1},1)),...
            ones(1,size(event.data{thispp}{2},1))*2,ones(1,size(event.data{thispp}{3},1))*3]';
        
        % calculate the amount of entries for each attention/feature
        % condition, and extract the location of the maximum size.
        attn.size{thispp} = [size(attn.data{thispp}{1},1), size(attn.data{thispp}{2},1), size(attn.data{thispp}{3},1)];
        [attn.maxsize{thispp},attn.sizeidx{thispp}] = max(attn.size{thispp});
        
        event.size{thispp} = [size(event.data{thispp}{1},1), size(event.data{thispp}{2},1), size(event.data{thispp}{3},1)];
        [event.maxsize{thispp}, event.sizeidx{thispp}] = max(event.size{thispp});
        
        % calculate the weights (matching the number of entries for each
        % condition) for the later weighted svm analysis. 
        attn.weights{thispp} = 1:size(attn.data{thispp},2);
        event.weights{thispp} = 1:size(event.data{thispp},2);
        
        attn.weights{thispp}(attn.sizeidx{thispp}) = [];
        event.weights{thispp}(event.sizeidx{thispp}) = [];
        
        attn.outputweights{thispp}(attn.sizeidx{thispp}) = 1;
        event.outputweights{thispp}(event.sizeidx{thispp}) = 1;
        
        for weightidx = attn.weights{thispp}
            attn.outputweights{thispp}(weightidx) = attn.maxsize{thispp}/size(cell2mat(attn.data{thispp}(weightidx)),1);
        end
        
        for weightidx = event.weights{thispp}
            event.outputweights{thispp}(weightidx) = event.maxsize{thispp}/size(cell2mat(event.data{thispp}(weightidx)),1);
        end
        
        % if we have first specified to train on attention data, perform
        % weighted support vector machine classification and train the
        % model to assess it's generalisation with the event condition
        % data. 
        if analysis == 1
            model = libsvmtrain(attn.labels{thispp}, attn.zscoredata{thispp}, sprintf('t -2 -w1 %.2f -w2 %.2f -w3 %.2f -q',...
                attn.outputweights{thispp}(1), attn.outputweights{thispp}(2), attn.outputweights{thispp}(3)));
            [label,acc,decision] = svmpredict(event.labels{thispp}, event.zscoredata{thispp}, model, ...
                sprintf('t -2 -w1 %.2f -w2 %.2f -w3 %.2f -q', event.outputweights{thispp}(1), event.outputweights{thispp}(2), event.outputweights{thispp}(3)));
        
        % perform the reverse analysis if specified to train on event data
        % first.
        elseif analysis == 2
            model = libsvmtrain(event.labels{thispp}, event.zscoredata{thispp}, sprintf('t -2 -w1 %.2f -w2 %.2f -w3 %.2f -q',...
                event.outputweights{thispp}(1), event.outputweights{thispp}(2), event.outputweights{thispp}(3)));
            [label,acc,decision] = svmpredict(attn.labels{thispp}, attn.zscoredata{thispp}, model,...
                sprintf('t -2 -w1 %.2f -w2 %.2f -w3 %.2f -q', attn.outputweights{thispp}(1), attn.outputweights{thispp}(2), attn.outputweights{thispp}(3)));
        end
        
        classification.acc{thispp} = acc(1);
    end
    
    % calculate the mean and standard error classification accuracies for
    % later plotting.
    output.mean = mean(cell2mat(classification.acc));
    output.stderr = std(cell2mat(classification.acc))/sqrt(size(classification.acc,2));
    
    % perform a one-sample wilcoxon signed-rank test of classification
    % versus chance.
    [output.p, ~, output.fullstats] = signrank(cell2mat(classification.acc),(100/size(attn.conds,2)), 'method', 'approximate');
    
    % store important values from this analysis. 
    ROIoutput.acc{thisROI} = classification.acc;
    ROIoutput.mean{thisROI} = output.mean;
    ROIoutput.stderr{thisROI} = output.stderr;
    ROIoutput.p{thisROI} = output.p;
    ROIoutput.stat{thisROI} = output.fullstats.zval;
end

% benjamini-hochberg significance values across ROIs. 
[~, ~, ~, ROIoutput.adjp]=fdr_bh(cell2mat(ROIoutput.p),0.05,'dep','yes');

% format classification data into header with column names and write to
% csv.
output.table = [num2cell(pps)', ROIoutput.acc{1}', ROIoutput.acc{2}', ROIoutput.acc{3}', ROIoutput.acc{4}', ROIoutput.acc{5}', ROIoutput.acc{6}'];
output.table = [output.table;[{'mean'}, ROIoutput.mean]];
output.table = [output.table;[{'stderr'}, ROIoutput.stderr]];
output.table = [output.table;[{'stat'}, ROIoutput.stat]];
output.table = [output.table;[{'p'}, ROIoutput.p]];
output.table = [output.table;[{'adj_p'}, num2cell(ROIoutput.adjp)]];
output.table = [[{'participants'}, ROIs]; output.table];

if analysis == 1
    cell2csv(strcat(output.directory,'attention-naturalistic_decoding.csv'), output.table); % write data to file.
elseif analysis == 2
    cell2csv(strcat(output.directory,'event-attention_decoding.csv'), output.table); % write data to file.
end

% plot this data as a bar chart with associated significance values and 
% save as .pdf. 
superbar(cell2mat(ROIoutput.mean), 'E', cell2mat(ROIoutput.stderr), 'P', ROIoutput.adjp,...
    'BarFaceColor', [.5 .5 .5], 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1.5,...
    'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0], 'PStarFontSize', 16,...
    'PStarShowGT', false, 'PStarShowNS', false);

set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs);
ylim([20 80]);
xlabel('Regions of Interest');
ylabel('Classification Accuracy (%)');

hline = refline([0 100/size(attn.conds,2)]); % add chance performance line.
hline.Color = 'r'; hline.LineWidth = 1; hline.LineStyle = '--';

if analysis == 1
    saveas(gcf, strcat(output.directory, 'attention-event_decoding.pdf'));
elseif analysis == 2
    saveas(gcf, strcat(output.directory, 'event-attention_decoding.pdf'));
end

% save analysis output values to .mat file. 
multivariate_eventattn_threeway.data = ROIoutput.acc;
multivariate_eventattn_threeway.means = ROIoutput.mean;
multivariate_eventattn_threeway.stderr = ROIoutput.stderr;
multivariate_eventattn_threeway.stats = ROIoutput.stat;
multivariate_eventattn_threeway.p = ROIoutput.p;
multivariate_eventattn_threeway.adj_p = ROIoutput.adjp;

save('/scratch/groups/Projects/P1323/original/fMRI/vista_output/output/multivariate_eventattn_threeway.mat', 'multivariate_eventattn_threeway');

%% --------------------------------------------------------------------- %%