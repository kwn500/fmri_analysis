%% ----- COMPARE ATTENTION VERSUS PASSIVE CORRELATION COEFFICIENTS ----- %%

% loads in condition-specific fisher-transformed correlation coefficients
% for each individual and assess difference between attention and passive
% conditions.

% still not particularly efficient, should be re-written. 

% 11/7/19 KWN

clear; clc;

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

exp_condition = 'colour'; % specify overall experiment to analyse. 
condition = 'interaction'; % specify condition to analyse. 

% for the interaction analysis, specify the interaction condition to
% analyse.
analysis_condition = 'luminance_feature'; 

% load in the condition-specific data.
if isequal(exp_condition, 'original')
    data = '/Users/kirstie/Documents/writeup/P1323feature/final/connectivity-fisher.mat';
elseif isequal(exp_condition, 'colour') && isequal(condition, 'feature')
    data = '/Users/kirstie/Documents/analysis/%s/final/feature/colour-featureconnectivity-fisher.mat';
elseif isequal(exp_condition, 'colour') && isequal(condition, 'colour-passive')
    data = '/Users/kirstie/Documents/analysis/colour/final/colour/colour-passiveconnectivity-fisher.mat';
elseif isequal(exp_condition, 'colour') && isequal(condition, 'interaction')
    data = sprintf('/Users/kirstie/Documents/analysis/colour/final/interaction/connectivity_interaction-fisher_%s.mat', analysis_condition);
elseif isequal(exp_condition, 'naturalistic')
    data = '/Users/kirstie/Documents/analysis/naturalistic/final/3feature/connectivity_naturalistic3feature-fisher.mat';
end
load(data);

%% ----------- PERFORM ATTENTION-PASSIVE DIFFERENCE ANALYSIS ----------- %%

if isequal(exp_condition, 'original') || isequal(condition, 'feature') || strfind(analysis_condition, 'feature') ~= 0
    
    % extract the overall atention data) (mean across all all attention conditions)
    % and the passive data. 
    attndata = squeeze(connectivity.indconditioncorrs.fisherdatasplit(1:3,:,:,:));
    attndata = squeeze(mean(attndata,1));
    
    passivedata = squeeze(connectivity.indconditioncorrs.fisherdatasplit(4,:,:,:));
    
    % for each attention condition individually, extract the
    % condition-specific data. 
    for thiscond = 1:3
        data{thiscond} = squeeze(connectivity.indconditioncorrs.fisherdatasplit(thiscond,:,:,:));
    end
    
    for thispp = 1:size(attndata,1) % for each participant in turn:
        
        % calculate the euclidean distance between the overall attention
        % data and passive data. 
        attnrmse(thispp) = sum((attndata(thispp,:)-passivedata(thispp,:)).^2);
        
        % for each attention condition in turn, calculate the rmse with the
        % passive data. 
        for thiscond = 1:3
            rmse(thiscond,thispp) = sum((data{thiscond}(thispp,:)-passivedata(thispp,:)).^2);
        end
    end
    
    % perform parametric and non-parametric t-tests for the overall attention
    % -passive distance data. 
    [wsr.attnp,~,wsr.attnstats] = signrank(attnrmse,0);
    [~,t.attnp,~,t.attnstats] = ttest(attnrmse);
    
    % for each attention condition in turn, repeat the same t-tests for
    % attention versus passive distance data. 
    for thiscond = 1:3
        [wsr.p{thiscond},~,wsr.stat{thiscond}] = signrank(squeeze(rmse(thiscond,:)),0);
        [~,t.p{thiscond},~,t.,stat{thiscond}] = ttest(squeeze(rmse(thiscond,:)));
    end
    
    % store the important analysis variables for writing to a .mat file. 
    connectivity.attentionpassive.wsr = wsr;
    connectivity.attentionpassive.t = t;
    connectivity.attentionpassive.attn.rmse = attnrmse;
    connectivity.attentionpassive.rmse = rmse;
    
% repeat the same process for the colour analysis. 
elseif isequal(condition, 'colour-passive') || strfind(analysis_condition, 'colour') ~= 0
    
    for thiscond = 1:6
        data{thiscond} = squeeze(connectivity.indconditioncorrs.fisherdatasplit(thiscond,:,:,:));
    end
    
    % calculate the difference between each chromatic attention and passive
    % condition. 
    for thispp = 1:size(rgadata,1)
        count = 1;
        for thiscond = 1:2:6
            rmse(count, thispp) = sum((data{thiscond}(thispp,:)-rgpdata(data{thiscond+1},:)).^2);
            count = count + 1;
        end
    end
    
    for thiscond = 1:3
        [wsr.p{thiscond},~,wsr.stat{thiscond}] = signrank(squeeze(rmse(thiscond,:)),0);
        [~,t.p{thiscond},~,t.stat{thiscond}] = ttest(squeeze(rmse(thiscond,:)));
    end
    
    connectivity.attentionpassive.wsr = wsr;
    connectivity.attentionpassive.t = t;
    connectivity.attentionpassive.rmse = rmse;
    
% repeat the same process for the naturalistic data. 
elseif isequal(exp_cond, 'naturalistic')
    
    pdata = squeeze(connectivity.indconditioncorrs.fisherdatasplit(5,:,:,:));
    
    for thiscond = 1:4
        data{thiscond} = squeeze(connectivity.indconditioncorrs.fisherdatasplit(thiscond,:,:,:));
    end
    
    for thispp = 1:size(fdata,1)
        for thiscond = 1:4
            rmse(thiscond,thispp) = sum((data{thiscond}(thispp,:)-pdata(thispp,:)).^2);
        end
    end
    
    for thiscond = 1:4
        [wsr.p{thiscond},~,wsr.stat{thiscond}] = signrank(squeeze(rmse(thiscond,:)),0);
        [~,t.p{thiscond},~,t.stat{thiscond}] = ttest(squeeze(rmse(thiscond,:)));
    end
    
    connectivity.attentionpassive.wsr = wsr;
    connectivity.attentionpassive.t = t;
    connectivity.attentionpassive.rmse = frmse;
end

% write the data to a .mat file. 
save('/Users/kirstie/Desktop/connectivity.mat', 'connectivity');
%% --------------------------------------------------------------------- %%