%% ---------- COMMUNICATION SUBSPACE REDUCED RANK REGRESSION ----------- %% 

% peforms reduced rank regression with cross validation on timeseries data 
% from our fMRI experiment using code from: https://github.com/joao-semedo/
% communication-subspace 

% 11/7/19 KWN

clear;clc;
SET_CONSTS
addpath(genpath('/scratch/sg3/P1323/code/fmri_analysis/functions'));

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %% 

% specify participants and conditions to analyse. 
participants = {'R2268','R2548', 'R2590', 'R2904', 'R3111', 'R3455',...
    'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

conditions = {'orientation', 'contrast', 'shape', 'passive'};

ROIs = {'V1', 'LO1'}; % specify ROIs to analyse.

%% ------------------ PERFORM REDUCED RANK REGRESSION ------------------ %%

for thiscond = 1:length(conditions) % for each condition in turn:
    for thispp = 1:length(participants) % for each participant in turn:
        participant = participants{thispp};
        
        % load in the ROI-specific data from both ROIs. 
        V1_ROI = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/processed_timeseries/%s/%s_%s_timeseries.mat',...
            ROIs{1}, participant, ROIs{1}));
        V2_ROI = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/processed_timeseries/%s/%s_%s_timeseries.mat',...
            ROIs{2}, participant, ROIs{2}));
        
        % extract the condition-specific data. 
        X = V1_ROI.tSeries.(conditions{thiscond});
        Y_V2 = V2_ROI.tSeries.(conditions{thiscond});
           
        % reshape and z-score across voxels.
        X_r = reshape(X,20,(length(X)/20),50);
        Y_V2_r = reshape(Y_V2,20,(length(X)/20),50);
        
        X_r = zscore(X_r);
        Y_V2_r= zscore(Y_V2_r);
         
        X = reshape(X_r,20*(length(X)/20), 50);
        Y_V2 = reshape(Y_V2_r,20*(length(X)/20), 50);
       
        % specify some default values for analysis and cross-validation. 
        numDimsUsedForPrediction = 1:10;
        cvNumFolds = 10;
        cvOptions = statset('crossval');
        
        % use parallel processing if available.
        cvOptions.UseParallel = true;
        
        % specify regression method to be used.
        regressMethod = @ReducedRankRegress;
        
        % provide train and test set inputs, it then fits the model to the
        % train set and uses it to predict the test set, reporting the
        % model's test performance using normalised squared error as the
        % performance metric. 
        cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
            (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
            numDimsUsedForPrediction, 'LossMeasure', 'NSE');
        
        % run the cross-validation routine. 
        cvl = crossval(cvFun, Y_V2, X, ...
            'KFold', cvNumFolds, ...
            'Options', cvOptions);
        
        % store the mean loss and standard error of the mean across folds
        % cross validation routines.
        cvLoss = [mean(cvl); std(cvl)/sqrt(cvNumFolds)];
        
        % compute the optimal dimensionality for the regression model. 
        optDimReducedRankRegress = ModelSelect...
            (cvLoss, numDimsUsedForPrediction);
        
        % plot the reduced rank regression cross-validation results. 
        x = numDimsUsedForPrediction;
        y = 1-cvLoss(1,:);
        e = cvLoss(2,:);
        errorbar(x, y, e, 'o--', 'Color', COLOR(V2,:), ...
            'MarkerFaceColor', COLOR(V2,:), 'MarkerSize', 10)
        xlabel('Number of predictive dimensions')
        ylabel('Predictive performance')
        
        % store this dimensionality for each condition and participant.
        optdimensions(thiscond,thispp) = optDimReducedRankRegress;
    end
end

%% ------- RUN ONE-WAY REPEATED MEASURES ANOVA ACROSS CONDITIONS ------- %% 

pp_number = repmat({'participant'},1,length(participants));

% create the data table, with participant, orientation, contrast and shape columns named accordingly. 
t = table(pp_number', optdimensions(1,:)', optdimensions(2,:)',...
    optdimensions(3,:)',optdimensions(4,:)',...
    'VariableNames', {'participant', 'Orientation', 'Contrast', 'Shape', 'Passive'});

% specify an overall name for the orientation, contrast and shape columns.
Meas = table([1 2 3 4]', 'VariableNames', {'Optimal_Dimensions'});

% fit a repeated measures model to this data.
rm = fitrm(t, 'Orientation-Passive~1', 'WithinDesign', Meas);

% run the repeated measures anova.
ranovatbl = ranova(rm);  
anova.full_output = ranovatbl;

% run associated mauchly's test of sphericity.
anova.mauchlytbl = mauchly(rm);
anova.mauchlyp = anova.mauchlytbl{1,4};

% save the important statistics to individually-named variables.
if anova.mauchlyp<.05
    anova.p = ranovatbl{1,6}; anova.f = ranovatbl{1,4}; anova.df_m = ranovatbl{1,2}; anova.df_r = ranovatbl{2,2};
else
    anova.p = ranovatbl{1,5}; anova.f = ranovatbl{1,4}; anova.df_m = ranovatbl{1,2}; anova.df_r = ranovatbl{2,2};
end

% run bonferroni-corrected post-hoc tests. 
anova.posthoc = multcompare(rm, 'Optimal_Dimensions', 'ComparisonType', 'bonferroni');

% save the data to a .mat file.
save(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/communication_subspace/reduced_rank_%s-%s.mat', ROIs{1}, ROIs{2}),...
    'optdimensions', 'anova');

%% --------------------------------------------------------------------- %%