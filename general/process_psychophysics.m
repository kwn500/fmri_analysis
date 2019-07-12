%% ---------------------- PROCESS PSYCHOPHYSICS ------------------------ %%

% updated script which loads in participants data from equal and unequal
% runs, extracts response and reaction time variables, removes data from
% participants without data from all conditions and calculates loglinear
% d', then performs wilcoxon signed ranks and repeated measures anovas.

% KWN 3/3/2019

clear; clc;
addpath(genpath('/Users/kirstie/Documents/code/fmri_analysis/functions'));

%% --------------------------- EXTRACT DATA ---------------------------- %%

runtype = {'equal', 'unequal'}; % specify run conditions.
conditions = {'orientation', 'contrast', 'shape', 'all'};
responses = {'Hit', 'Miss', 'False Alarm', 'Correct Rejection'};
outputdir = '/Users/kirstie/Documents/analysis/original/psychophysics/output/';

for thisrun = 1:length(runtype) % for each run condition:
    
    run = runtype{thisrun}; % specify the run name.
    
    % extract the filenames corresponding to all runs for all participants.
    filedir = dir(sprintf('/Users/kirstie/Documents/analysis/original/psychophysics/%s/*.csv', run)); 
    filenames = {filedir.name};
    ppcount = 1; % initialise participant counter variable.
    
    for thispp = 1:2:length(filenames) % for each set of participant runs (2 per participant).
        
        % load in data across the two runs for each participant.
        run1 = readtable(sprintf('/Users/kirstie/Documents/analysis/original/psychophysics/%s/%s',...
            run, filenames{thispp}));
        run2 = readtable(sprintf('/Users/kirstie/Documents/analysis/original/psychophysics/%s/%s',...
            run, filenames{thispp+1}));
        
        % extract only the columns of interest for further processing.
        run1edit = table2cell(run1(:,[2,7,12,17,20,21,22,23]));
        run2edit = table2cell(run2(:,[2,7,12,17,20,21,22,23]));
        
        % combine the data into a singular matrix.
        data(thisrun,ppcount,:,:) = [run1edit;run2edit];
        
        % increment participant counter variable.
        ppcount = ppcount +1;
    end
end

%% -------------------- PROCESS OUTLIERS/CLEAN DATA -------------------- %%

for thisrun = 1:size(data,1) % for each run in turn:
    for thispp = 1:size(data,2) % for each participant in turn:
        
        % intialise outlier counter variable.
        rtoutliercount = 0;
        
        for thiscond = 1:size(conditions,2) % for each condtion in turn:
            cond = conditions{thiscond};
            
            % extract the condition-specific data.
            ind{thisrun, thispp, thiscond} = find(contains(squeeze(data(thisrun,thispp,:,1)),cond));
            conddata{thisrun,thispp,thiscond} = squeeze(data(thisrun,thispp,ind{thisrun,thispp,thiscond},:));
            
            tempdp = squeeze(data(thisrun,thispp,ind{thisrun,thispp,thiscond},:));
          
            % remove any row of reaction time smaller than 150ms (this is
            % not thought to represent a plausible response).
            temprt = cell2mat(squeeze(data(thisrun,thispp,ind{thisrun,thispp,thiscond},6)));
            outlier = temprt<150;
            temprt(outlier) = [];
            if sum(outlier) > 0
                rtoutliercountpp(thiscond, thisrun,thispp) = rtoutliercount + 1;
            end
            
            % remove these same outlier rows from the participant's
            % corresponding dprime data.
            tempdp(outlier,:) = [];
            conddata_rt{thisrun,thispp,thiscond} = temprt;
            conddata{thisrun,thispp,thiscond} = tempdp;
        end
    end
end

%% --------------------------- BOOTSTRAPPING --------------------------- %%

for thisrun = 1:size(conddata,1) % for each run in turn:
        for thiscond = 1:size(conddata,3) % for each condition in turn:
            count = 0;
            for thispp = 1:size(conddata,2) % for each participant in turn:
            
            % combine data across runs for each condition.
            if ~isempty(cell2mat(ind(thisrun,thispp,thiscond)))
                if count == 0
                    conddata_bstrp{thiscond} = conddata{thisrun,thispp,thiscond};
                else
                    conddata_bstrp{thiscond} = [conddata_bstrp{thiscond}; conddata{thisrun,thispp,thiscond}];
                end
                count = count + 1;
            end
        end
    end
end

bstrp_samp = 50; % specify the number of required bootstrapped samples.

for thisrun = 1:size(conddata,1) % for each run in turn:
    for thiscond = 1:size(conddata,3) % for each condition in turn:
        for thispp = 1:size(conddata,2) % for each participant in turn:
            
            % if there is no data corresponding to this particular
            % condition.
            if isempty(conddata{thisrun,thispp,thiscond})
                
                % extract the condition-specific data across all runs,
                % created earlier.
                dataselect = 1:size(conddata_bstrp{thiscond},1);
                
                % take a sample of this data with replacement. 
                sample = datasample(dataselect,bstrp_samp, 'Replace',true);
                
                % add this boostrapped data to the corresponding condition
                % and run position. 
                conddata{thisrun,thispp,thiscond} = conddata_bstrp{thiscond}(sample,:);
                conddata_rt{thisrun,thispp,thiscond} = cell2mat(conddata_bstrp{thiscond}(sample,6));
            end
        end
    end
end

%% ------------------------- CALCULATE D PRIME ------------------------- %%

for thisrun = 1:size(data,1) % for each run in turn:
    for thispp = 1:size(data,2) % for each participant in turn:
        
        for thiscond = 1:size(conditions,2) % for each condition:
              for thisresponse = 1:size(responses,2) % for each response type in turn:
                response = responses{thisresponse};
                
                % extract the data and timing information corresponding to
                % each response type.
                responseind{thisrun,thispp,thiscond,thisresponse} = find(contains(squeeze(conddata{thisrun,thispp,thiscond}(:,5)), response));
                responseindrt{thisrun,thispp,thiscond,thisresponse} = find(contains(squeeze(conddata{thisrun,thispp,thiscond}(:,5)), response));
                responsedata{thisrun,thispp,thiscond,thisresponse} = squeeze(conddata{thisrun,thispp,thiscond}(responseind{thisrun,thispp,thiscond,thisresponse},:));
                responsedata_rt{thisrun,thispp,thiscond,thisresponse} = squeeze(conddata_rt{thisrun,thispp,thiscond}(responseind{thisrun,thispp,thiscond,thisresponse}));
              end

            % calculate the number of signal and noise trials and increment
            % by 1.
            signal = (size(responsedata{thisrun,thispp,thiscond,1},1)+size(responsedata{thisrun,thispp,thiscond,2},1)) + 1;
            noise = (size(responsedata{thisrun,thispp,thiscond,3},1)+size(responsedata{thisrun,thispp,thiscond,4},1)) + 1;
            
            % calculate the number of hit and false alarm responses and
            % increment by 0.5.
            hit = (size(responsedata{thisrun,thispp,thiscond,1},1))+0.5;
            fa = (size(responsedata{thisrun,thispp,thiscond,3},1))+0.5;
            
            % calculate loglinear d prime.
            dprime(thisrun,thispp,thiscond) = norminv(hit/signal)-norminv(fa/noise);
        end
    end
end

%% ----------------------- EQUAL-UNEQUAL DPRIME ------------------------ %%

for thiscond = 1:size(conditions,2)+1 % for each condition in turn:
    
    % either extract the condition-specific data, or if analysing all
    % conditions, extract all data for each run. 
    if thiscond ~= 5
        equal_dp = squeeze(dprime(1,:,thiscond));
        unequal_dp = squeeze(dprime(2,:,thiscond));
    elseif thiscond == 5
        equal_dp = mean(squeeze(dprime(1,:,:)),2);
        unequal_dp = mean(squeeze(dprime(2,:,:)),2);
    end
    
    % calculate the difference between equal and unequal dprime using the
    % appropriate parametric or non-parametric test. 
    if thiscond == 1 || thiscond == 4
        [dp_p, ~ ,dpstats] = signrank(equal_dp, unequal_dp, 'method', 'approximate');
    else
        [~,dp_p,~,dpstats] = ttest(equal_dp, unequal_dp);
    end
    
    % calculate the mean and standard error of equal and unequal
    % distributions for later plotting.
    plotmatrix = [mean(equal_dp), mean(unequal_dp)];
    stderrmatrix = [std(equal_dp)/sqrt(length(equal_dp)),std(unequal_dp)/sqrt(length(unequal_dp))];
    
    % extract the important statistical tests results and perform 
    % shapiro-wilk normality tests.
    % we plot this data on a bar-chart with significance asterisks and save
    % the figure to a file.
    if thiscond ~=5
        makefigure(plotmatrix,stderrmatrix,dp_p,[0,3],1:2,runtype,outputdir,sprintf('equal-unequal_dp_%s', conditions{thiscond}));
        dprimeoutput.equalunequal.(conditions{thiscond}).p = dp_p;
        dprimeoutput.equalunequal.(conditions{thiscond}).stats = dpstats;
        dprimeoutput.equalunequal.(conditions{thiscond}).equalmean = mean(equal_dp);
        dprimeoutput.equalunequal.(conditions{thiscond}).unequalmean = mean(unequal_dp);
        
        [~,dprimeoutput.equalunequal.(conditions{thiscond}).equalswp] = swtest(equal_dp);
        [~,dprimeoutput.equalunequal.(conditions{thiscond}).unequalswp] = swtest(unequal_dp);
        
    elseif thiscond == 5
        makefigure(plotmatrix,stderrmatrix,dp_p,[0,3],1:2,runtype,outputdir,'equal-unequal_dp_overall');
        dprimeoutput.equalunequal.overall.p = dp_p;
        dprimeoutput.equalunequal.overall.stats = dpstats;
        dprimeoutput.equalunequal.overall.equalmean = mean(equal_dp);
        dprimeoutput.equalunequal.overall.unequalmean = mean(unequal_dp);
        
        [~,dprimeoutput.equalunequal.overall.equalswp] = swtest(equal_dp);
        [~,dprimeoutput.equalunequal.overall.unequalswp] = swtest(unequal_dp);
    end
end

% perform benjamini-hochberg correction on the resulting p-values across
% all conditions. 
dprime_equalunequal_p = [dprimeoutput.equalunequal.all.p,dprimeoutput.equalunequal.orientation.p,...
    dprimeoutput.equalunequal.contrast.p,dprimeoutput.equalunequal.shape.p];
[~, ~, ~, adj_p] = fdr_bh(dprime_equalunequal_p);
dprimeoutput.equalunequal.adjplabels = {'all', 'orientation', 'contrast', 'shape'};
dprimeoutput.equalunequal.adjp = adj_p;


%% -------------------- EQUAL-UNEQUAL REACTION TIME -------------------- %%

% repeat the same process as before with the reaction time data.
for thiscond = 1:size(conditions,2)+1
    
    if thiscond ~= 5
        equal_rt = cellfun(@nanmean,squeeze(conddata_rt(1,:,:,:)));
        equal_rt = equal_rt(:,thiscond);
        
        unequal_rt = cellfun(@nanmean,squeeze(conddata_rt(2,:,:,:)));
        unequal_rt = unequal_rt(:,thiscond);
        
    elseif thiscond == 5
        equal_rt = nanmean(cellfun(@nanmean,squeeze(conddata_rt(1,:,:,:))),2);
        unequal_rt = nanmean(cellfun(@nanmean,squeeze(conddata_rt(2,:,:,:))),2);
    end
    
    if thiscond == 2
        [rt_p, ~ ,rtstats] = signrank(equal_rt, unequal_rt, 'method', 'approximate');
    else 
        [~,rt_p,~,rtstats] = ttest(equal_rt, unequal_rt);
    end
    
    plotmatrix = [nanmean(equal_rt), nanmean(unequal_rt)];
    stderrmatrix = [nanstd(equal_rt)/sqrt(length(equal_rt)),nanstd(unequal_rt)/sqrt(length(unequal_rt))];
    
    if thiscond ~=5
        makefigure(plotmatrix,stderrmatrix,rt_p,[0,500],1:2,runtype,outputdir,sprintf('equal-unequal_rt_%s', conditions{thiscond}));
        rtoutput.equalunequal.(conditions{thiscond}).p = rt_p;
        rtoutput.equalunequal.(conditions{thiscond}).stats = rtstats;
        rtoutput.equalunequal.(conditions{thiscond}).equalmean = nanmean(equal_rt);
        rtoutput.equalunequal.(conditions{thiscond}).unequalmean = nanmean(unequal_rt);
        
        [~,rtoutput.equalunequal.(conditions{thiscond}).equalswp] = swtest(equal_rt);
        [~,rtoutput.equalunequal.(conditions{thiscond}).unequalswp] = swtest(unequal_rt);
    elseif thiscond == 5
        makefigure(plotmatrix,stderrmatrix,rt_p,[0,500],1:2,runtype,outputdir,'equal-unequal_rt_overall');
        rtoutput.equalunequal.overall.p = rt_p;
        rtoutput.equalunequal.overall.stats = rtstats;
        rtoutput.equalunequal.overall.equalmean = nanmean(equal_rt);
        rtoutput.equalunequal.overall.unequalmean = nanmean(unequal_rt);
        
        [~,rtoutput.equalunequal.overall.equalswp] = swtest(equal_rt);
        [~,rtoutput.equalunequal.overall.unequalswp] = swtest(unequal_rt);
    end
end

rt_equalunequal_p = [rtoutput.equalunequal.all.p,rtoutput.equalunequal.orientation.p,...
    rtoutput.equalunequal.contrast.p,rtoutput.equalunequal.shape.p];
[~, ~, ~, adj_p] = fdr_bh(rt_equalunequal_p);
rtoutput.equalunequal.adjplabels = {'all', 'orientation', 'contrast', 'shape'};
rtoutput.equalunequal.adjp = adj_p;

%% ------------------- SELECTIVE-NONSELECTIVE DPRIME ------------------- %%

for thisrun = 1:length(runtype) % for each run in turn:
    
    % extract the selective and non-selective d prime data.
    selective_dp = nanmean(squeeze(dprime(thisrun,:,[1:3])),2);
    nonselective_dp = squeeze(dprime(thisrun,:,4));
    
    % perform a wilcoxon signed-rank tests comparing selective and
    % non-selective data. 
    [dp_p, ~ ,dpstats] = signrank(selective_dp, nonselective_dp, 'method', 'approximate');
    
    % calculate the mean and standard error of the selective and
    % non-selective data for later plotting. 
    plotmatrix = [nanmean(selective_dp), nanmean(nonselective_dp)];
    stderrmatrix = [nanstd(selective_dp)/sqrt(length(nonselective_dp)),nanstd(selective_dp)/sqrt(length(nonselective_dp))];
    
    % plot the data as a bar chart with significance asterisks and save t
    % to a file.
    makefigure(plotmatrix,stderrmatrix,dp_p,[0,3],1:2,{'selective', 'non-selective'},...
        outputdir,sprintf('selective-nonselective_dp_overall-%s', runtype{thisrun}));
    
    % extract rthe important analysis output statistics and perform
    % shapiro-wilk normality tests. 
    dprimeoutput.selectivenonselective.(runtype{thisrun}).p = dp_p;
    dprimeoutput.selectivenonselective.(runtype{thisrun}).stats = dpstats;
    dprimeoutput.selectivenonselective.(runtype{thisrun}).selectivemean = mean(selective_dp);
    dprimeoutput.selectivenonselective.(runtype{thisrun}).nonselectivemean = mean(nonselective_dp);
    
    [~,dprimeoutput.selectivenonselective.(runtype{thisrun}).selectiveswp] = swtest(selective_dp);
    [~,dprimeoutput.selectivenonselective.(runtype{thisrun}).nonselectiveswp] = swtest(nonselective_dp);
    
end

%% ---------------- SELECTIVE-NONSELECTIVE REACTION TIME --------------- %%

% repeat the same process as above for the selective versus non-selective
% reaction time data.
for thisrun = 1:length(runtype)
    
    selective_rt = nanmean(cellfun(@nanmean,squeeze(conddata_rt(thisrun,:,[1:3],:))),2);
    nonselective_rt = (cellfun(@nanmean,squeeze(conddata_rt(thisrun,:,4,:))))';
    
    [~,rt_p,~,rtstats] = ttest(selective_rt, nonselective_rt);

    plotmatrix = [nanmean(selective_rt), nanmean(nonselective_rt)];
    stderrmatrix = [nanstd(selective_rt)/sqrt(length(nonselective_rt)),nanstd(selective_rt)/sqrt(length(nonselective_rt))];
    
    makefigure(plotmatrix,stderrmatrix,rt_p,[0,500],1:2,{'selective', 'non-selective'},outputdir,sprintf('selective-nonselective_rt_overall-%s', runtype{thisrun}));
    
    rtoutput.selectivenonselective.(runtype{thisrun}).p = rt_p;
    rtoutput.selectivenonselective.(runtype{thisrun}).stats = rtstats;
    rtoutput.selectivenonselective.(runtype{thisrun}).selectivemean = nanmean(selective_rt);
    rtoutput.selectivenonselective.(runtype{thisrun}).nonselectivemean = nanmean(nonselective_rt);
    
    [~,rtoutput.selectivenonselective.(runtype{thisrun}).selectiveswp] = swtest(selective_rt);
    [~,rtoutput.selectivenonselective.(runtype{thisrun}).nonselectiveswp] = swtest(nonselective_rt);
end

%% -------------- ORIENTATION, CONTRAST, SHAPE D PRIME ----------------- %%

for thisrun = 1:length(runtype) % for each run in turn:
    
    % extract the feature-specific d prime data. 
    orientation_dp = squeeze(dprime(thisrun,:,1));
    contrast_dp = squeeze(dprime(thisrun,:,2));
    shape_dp = squeeze(dprime(thisrun,:,3));
    
    % perform a one-way repeated measures anova and associated
    % bonferroni-corrected post-hoc tests. 
    ppnumber = repmat({'participant'}, size(orientation_dp,2),1);
    t = table(ppnumber, orientation_dp', contrast_dp', shape_dp',...
        'VariableNames', {'participant', 'Orientation', 'Contrast', 'Shape'});
    Meas = table([1 2 3]', 'VariableNames', {'Loglinear_dprime'});
    rm = fitrm(t, 'Orientation-Shape~1', 'WithinDesign', Meas);
    dp_mauchly = mauchly(rm);
    dp_ranovatbl = ranova(rm);
    dp_posthoc = multcompare(rm, 'Loglinear_dprime', 'ComparisonType', 'bonferroni');
    
    % extract the correct-p value for examination in line with the output
    % of mauchly sphericity tests. 
    if cell2mat(table2cell(dp_mauchly(1,4))) <.05
        dp_p = cell2mat(table2cell(dp_ranovatbl(1,6)));
    else
        dp_p = cell2mat(table2cell(dp_ranovatbl(1,5)));
    end
    
    % calculate the mean and standard error of each condition for later
    % plotting.
    plotmatrix = [mean(orientation_dp), mean(contrast_dp), mean(shape_dp)];
    stderrmatrix = [nanstd(orientation_dp)/sqrt(length(orientation_dp)),...
        nanstd(contrast_dp)/sqrt(length(contrast_dp)), nanstd(shape_dp)/sqrt(length(shape_dp))];
    
    % plot the data as a bar chart with significance asterisks and save to
    % file.
    makefigure(plotmatrix,stderrmatrix,dp_p,[0,3],1:3,{'Orientation', 'Contrast', 'Shape'},...
        outputdir,sprintf('ocs_dp-%s', runtype{thisrun}), dp_posthoc);
    
    % store the important analysis data and perform shapiro wilk normality
    % tests.
    dprimeoutput.ocs.(runtype{thisrun}).mauchly = dp_mauchly;
    dprimeoutput.ocs.(runtype{thisrun}).rmanova = dp_ranovatbl;
    dprimeoutput.ocs.(runtype{thisrun}).posthoc = dp_posthoc;
    dprimeoutput.ocs.(runtype{thisrun}).omean = nanmean(orientation_dp);
    dprimeoutput.ocs.(runtype{thisrun}).cmean = nanmean(contrast_dp);
    dprimeoutput.ocs.(runtype{thisrun}).smean = nanmean(shape_dp);
    
    [~,dprimeoutput.ocs.(runtype{thisrun}).orientationswp] = swtest(orientation_dp);
    [~,dprimeoutput.ocs.(runtype{thisrun}).contrastswp] = swtest(contrast_dp);
    [~,dprimeoutput.ocs.(runtype{thisrun}).shapeswp] = swtest(shape_dp);
end

%% ----------- ORIENTATION, CONTRAST, SHAPE REACTION TIME -------------- %%

% repeat the same process for the orientation, contrast and shape reaction
% time data. 
for thisrun = 1:length(runtype)
    
    orientation_rt = cellfun(@nanmean,squeeze(conddata_rt(thisrun,:,1,:)))';
    contrast_rt = cellfun(@nanmean,squeeze(conddata_rt(thisrun,:,2,:)))';
    shape_rt = cellfun(@nanmean,squeeze(conddata_rt(thisrun,:,3,:)))';
    
    ppnumber = repmat({'participant'}, size(orientation_rt,1),1);
    t = table(ppnumber, orientation_rt, contrast_rt, shape_rt,...
        'VariableNames', {'participant', 'Orientation', 'Contrast', 'Shape'});
    Meas = table([1 2 3]', 'VariableNames', {'Reaction_Time'});
    rm = fitrm(t, 'Orientation-Shape~1', 'WithinDesign', Meas);
    rt_mauchly = mauchly(rm);
    rt_ranovatbl = ranova(rm);
    rt_posthoc = multcompare(rm, 'Reaction_Time', 'ComparisonType', 'bonferroni');
    
    if cell2mat(table2cell(rt_mauchly(1,4))) <.05
        rt_p = cell2mat(table2cell(rt_ranovatbl(1,6)));
    else
        rt_p = cell2mat(table2cell(rt_ranovatbl(1,5)));
    end
    
    plotmatrix = [mean(orientation_rt), mean(contrast_rt), mean(shape_rt)];
    stderrmatrix = [nanstd(orientation_rt)/sqrt(length(orientation_rt)),...
        nanstd(contrast_rt)/sqrt(length(contrast_rt)), nanstd(shape_rt)/sqrt(length(shape_rt))];
    
    makefigure(plotmatrix,stderrmatrix,rt_p,[0,500],1:3,{'Orientation', 'Contrast', 'Shape'},...
        outputdir,sprintf('ocs_rt-%s', runtype{thisrun}), rt_posthoc);
    
    rtoutput.ocs.(runtype{thisrun}).mauchly = rt_mauchly;
    rtoutput.ocs.(runtype{thisrun}).rmanova = rt_ranovatbl;
    rtoutput.ocs.(runtype{thisrun}).posthoc = rt_posthoc;
    rtoutput.ocs.(runtype{thisrun}).omean = nanmean(orientation_rt);
    rtoutput.ocs.(runtype{thisrun}).cmean = nanmean(contrast_rt);
    rtoutput.ocs.(runtype{thisrun}).smean = nanmean(shape_rt);
    
    [~,rtoutput.ocs.(runtype{thisrun}).orientationswp] = swtest(orientation_rt);
    [~,rtoutput.ocs.(runtype{thisrun}).contrastswp] = swtest(contrast_rt);
    [~,rtoutput.ocs.(runtype{thisrun}).shapeswp] = swtest(shape_rt);
end

% save the d prime and reaction time output across all analyses to a .mat
% file.
save(strcat(outputdir,'psychophysics_analysis_output.mat'), 'dprimeoutput', 'rtoutput');

%%---------------------------------------------------------------------- %%