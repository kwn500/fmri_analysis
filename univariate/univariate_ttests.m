%% ---------- VERSUS PASSIVE & VERSUS ZERO UNIVARIATE T-TESTS ---------- %% 

% runs t-tests to compare activation in attention versus passive conditions
% and assesses if BOLD signal modulation is significantly greater than zero
% for all experiments. 

% 10/7/19 KWN

clear; clc

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %% 

ROIs = {'V1', 'V3AB', 'V4', 'LO1', 'LO2', 'A1'}; % specify the ROIs to analyse. 
exp_condition = 'original'; % specify the experiment to analyse. 

% specify the analysis conditions to process, for the naturalistic
% 3x3feature analysis, we specify the condition name (e.g. orientation3x3).
conditions = {'block', 'event'};

% for naturalistic, state whether to analyse face data (1 = yes, 0 = no). 
face_toggle = 1;

%% ------------------------ RUN T-TEST ANALYSES ------------------------ %% 

for thiscond = 1:length(conditions) % for each condition in turn:
    condition = conditions{thiscond};
    
    for thisROI = 1:length(ROIs) % for each ROI in turn:
        ROI = ROIs{thisROI};
        
        % load in the ROI and condition-specific univariate beta-values.
        if ~isequal(exp_condition, 'naturalistic')
            load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/group/betas/%s_%s_univariate_betas.mat',...
                condition, condition, ROI, condition));
        elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
            load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3feature/STS_OFA/group/%s_univariate_betas.mat', ROI));
        elseif isequal(exp_condition, 'naturalistic') && ~isequal(condition, '3feature')
            load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3x3feature/%s/group/%s_%s_univariate_betas.mat',...
                condition, ROI, condition));
        end
        
        % calculate the mean univariate betas across all attention
        % conditions, and extract the corresponding passive beta values.
        if ~isequal(exp_condition, 'naturalistic')
            attnbetas(thiscond,thisROI,:) = mean(betas_pp(:,[1:3]),2);
            passivebetas(thiscond,thisROI,:) = betas_pp(:,4);
            
        elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
            if face_toggle == 1
                attnbetas(thiscond,thisROI,:) = mean(betas_pp(:,1:4),2);
            elseif face_toggle == 0
                attnbetas(thiscond,thisROI,:) = mean(betas_pp(:,2:4),2);
            end
            passivebetas(thisROI,:) = betas_pp(:,5);
        elseif isequal(exp_condition, 'naturalistic') && ~isequal(condition, '3feature')
            attnbetas(thiscond,thisROI,:) = mean(betas_pp,2);
        end
        
        if isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature') || isequal(exp_condition, 'original') || isequal(exp_condition, 'colour')
            % calculate the difference in BOLD signal modulation between the
            % attention and passive conditions.
            difference(thiscond,thisROI,:) = squeeze(attnbetas(thiscond,thisROI,:))-squeeze(passivebetas(thiscond,thisROI,:));
            
            % perform a wilcoxon signed-rank test comparing activation across
            % attention and passive conditions.
            [passivep(thiscond,thisROI),~,passivestats] = signrank(squeeze(attnbetas(thiscond,thisROI,:)),...
                squeeze(passivebetas(thiscond,thisROI,:)),'method', 'approximate');
            passivestat(thiscond,thisROI) = passivestats.zval;
            
            % perform a wilcoxon signed-rank test comparing activation in the
            % attention condition versus zero.
            [p(thiscond, thisROI),~,stats] = signrank(squeeze(attnbetas(thiscond,thisROI,:)), 0, 'method', 'approximate');
            stat(thiscond,thisROI) = stats.zval;
            
        elseif isequal(exp_condition, 'naturalistic') && ~isequal(condition, '3feature')
            [p(thiscond, thisROI),~,stats] = signrank(squeeze(attnbetas(thiscond,thisROI,:)), 0, 'method', 'approximate');
            stat(thiscond,thisROI) = stats.zval;
        end
    end
    
    % perform benjamini-hochberg correction on the significance values
    % across ROIs for the versus-zero and versus-passive conditions.
    [~,~,~,adj_p(thiscond,:)]=fdr_bh(p(thiscond,:));
    [~,~,~,passiveadj_p(thiscond,:)]=fdr_bh(passivep(thiscond,:));
    
    % calculate the mean difference and standard error for the attention
    % versus passive difference. 
    diffmean(thiscond,:) = mean(squeeze(difference(thiscond,:,:)),2);
    diffstd(thiscond,:) = std(squeeze(difference(thiscond,:,:)),[],2)/sqrt(length(difference(thiscond,:,:)));
    
    % plot these difference values for each ROI on a bar chart with
    % associated versus-passive significance values.
    figure(thiscond)
    superbar(diffmean(thiscond,:), 'E', diffstd(thiscond,:), 'P', passiveadj_p(thiscond,:),...
        'BarFaceColor', [.5 .5 .5], 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1.5,...
        'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0], 'PStarFontSize', 16,...
        'PStarShowGT', false, 'PStarShowNS', false);
    set(gca, 'xtick', 1:length(ROIs), 'xticklabel', ROIs);
    ylim([-.25 .25]);
    xlabel('Regions of Interest');
    ylabel('Attention-Passive Difference');
    
    % save the figure to a .pdf file. 
    if ~isequal(exp_condition, 'naturalistic')
        saveas(gcf, strcat(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/univariate/%s/',...
            condition, condition), sprintf('attn-passive_%s', condition), '.pdf'));
    else
        saveas(gcf, strcat(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/univariate/3feature/'),...
            sprintf('attn-passive_naturalistic'), '.pdf'));
    end
end

% save the important analysis values for the versus-zero and versus-passive
% analyses to associated struct.
univariate_wsr.vszero.order = conditions;
univariate_wsr.vszero.p = p;
univariate_wsr.vszero.adj_p = adj_p;
univariate_wsr.vszero.zval = stat;
univariate_wsr.vszero.attndata = attnbetas;

univariate_wsr.vspassive.order = conditions;
univariate_wsr.vspassive.p = passivep;
univariate_wsr.vspassive.adj_p = passiveadj_p;
univariate_wsr.vspassive.zval = passivestat;
univariate_wsr.vspassive.attndata = attnbetas;
univariate_wsr.vspassive.passivedata = passivebetas;

% save these analysis structs to a .mat file. 
if ~isequal(exp_condition, 'naturalistic')
    save('/scratch/groups/Projects/P1323/original/fMRI/vista_output/output/univariate_wsr.mat', 'univariate_wsr');
elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
    save(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/output/3feature/univariate_wsr.mat'), 'univariate_wsr');
elseif isequal(exp_condition, 'naturalistic') && ~isequal(condition, '3feature')
    save(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/output/3x3feature/%s_univariate_wsr.mat', condition), 'univariate_wsr');
end

%% --------------------------------------------------------------------- %%