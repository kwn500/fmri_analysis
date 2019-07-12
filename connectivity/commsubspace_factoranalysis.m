%% -------------- COMMUNICATION SUBSPACE FACTOR ANALYSIS --------------- %% 

% peforms factor analysis with cross validation on timeseries data from our
% fMRI experiment using code from: https://github.com/joao-semedo/
% communication-subspace 

% 11/7/19 KWN

clear;clc;
SET_CONSTS
addpath(genpath('/scratch/sg3/P1323/code/fmri_analysis/functions'));

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %% 

% specify participants and conditions to analyse. 
participants = {'R2268','R2548', 'R2590', 'R2904', 'R3111',...
    'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

conditions = {'orientation', 'contrast', 'shape', 'passive'};

ROI = 'V1'; % specify ROI to analyse.

%% ---------------------- PERFORM FACTOR ANALYSIS ---------------------- %%

for thiscond = 1:length(conditions) % for each condition in turn:
    for thispp = 1:length(participants) % for each participant in turn:
        participant = participants{thispp};
        
        % load in the ROI of interest data.
        V1_ROI = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/vista_output/multivariate/block/processed_timeseries/%s/%s_%s_timeseries.mat',...
            ROI, participant, ROI));
      
        % extract the condition-specific data.
        X = V1_ROI.tSeries.(conditions{thiscond});
        
        % specify some default values for analysis and cross-validation. 
        q = 0:30;
        cvNumFolds = 10;
        cvOptions = statset('crossval');
        
        % use parallel processing if available.
        cvOptions.UseParallel = true;
        
        % return the cumulative shared variance explained. 
        cvLoss= CrossValFa(X, q, cvNumFolds, cvOptions);
        
        % compute optimal factor analysis dimensionality.
        qOpt = FactorAnalysisModelSelect(cvLoss, q);
        
        % store this dimensionality for each condition and participant. 
        optdimensions(thiscond,thispp) = qOpt;
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
save('/scratch/groups/Projects/P1323/original/fMRI/vista_output/communication_subspace/factor_analysis.mat',...
    'optdimensions', 'anova');
%% --------------------------------------------------------------------- %%