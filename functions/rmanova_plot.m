function rmanova_plot(input_matrix, input_ROI, analysis_condition, toggle, stats_output_dir, figure_output_dir)

% performs one-way repeated measures ANOVA with associated post-hoc tests
% and plot the results as a bar-chart with error bars and associated
% significance values before saving to a .pdf figure. 

% 11/7/19 KWN

%% ---------------- RUN ONE-WAY REPEATED MEASURES ANOVA ---------------- %%

% select only the orientation, contrast and shape (or RG, BY and
% luminance) columns for further analysis. 
if isequal(analysis_condition, 'colourxfeature')
    plot_matrix = input_matrix; stderr_matrix = input_matrix;
else
    plot_matrix = input_matrix(:,1:end-1); stderr_matrix = input_matrix(:,1:end-1);
end

% if we have specified an analysis on normalised data:
if isequal(toggle, 'normalised') 
    
    % subtract the participant-specific passive activation from the
    % remaining three columns of interest. 
    passive = repmat(input_matrix(:,4),1,3);
    normalised = plot_matrix - passive;
    plot_matrix = normalised;
end

% calculate the overall mean and standard error across all participants for each of the conditions.
plot_matrix = nanmean(plot_matrix); stderr_matrix = nanstd(stderr_matrix) / sqrt(length(stderr_matrix));
 
% create a column the length of the number of participants, containing the 'participant' string. 
pp_number = repmat({'participant'}, size(input_matrix,1),1); 

% if this is a 'feature' analysis, either from the original or colour
% experiment: 
if isequal(analysis_condition, 'block') || isequal(analysis_condition, 'feature') || ~isempty(strfind(condition_name,'feature'));  
    
    % create the data table, with participant, orientation, contrast and shape columns named accordingly. 
    t = table(pp_number, input_matrix(:,1), input_matrix(:,2), input_matrix(:,3),...
        'VariableNames', {'participant', 'Orientation', 'Contrast', 'Shape'});
    
    % specify an overall name for the orientation, contrast and shape columns. 
    Meas = table([1 2 3]', 'VariableNames', {'Amplitude_Modulation'}); 
    
    % fit a repeated measures model to this data. 
    rm = fitrm(t, 'Orientation-Shape~1', 'WithinDesign', Meas);
    
    % specify the corresponding x axis labels for later plotting. 
    xtick_labels = {'Orientation', 'Contrast', 'Shape'};
    
    % specify corresponding bonferroni post-hoc comparison names for later
    % writing to a .csv file header. 
    headers = {'rm_p', 'rm_f', 'rmdf_m', 'rmdf_r', 'Bo_vs_c', 'Bo_vs_s', 'Bc_vs_s'};
    
% otherwise, if this is an analysis comparing the three colour conditions. 
elseif isequal(analysis_condition, 'colour') || ~isempty(strfind(condition_name,'colour'));  
    
    % repeat the same process as above using RG, BY and LUM variable names
    % as oppose to orientation, contrast and shape. 
    t = table(pp_number, input_matrix(:,1), input_matrix(:,2), input_matrix(:,3),...
        'VariableNames', {'participant', 'RG', 'BY', 'LUM'});

    Meas = table([1 2 3]', 'VariableNames', {'Amplitude_Modulation'});  
    rm = fitrm(t, 'RG-LUM~1', 'WithinDesign', Meas); 
    xtick_labels = {'RG', 'BY', 'LUM'};
    headers = {'rm_p', 'rm_f', 'rmdf_m', 'rmdf_r', 'Brg_vs_by', 'Brg_vs_lum', 'Bby_vs_lum'};
end 
    
% run the one-way repeated measures anova.    
ranovatbl = ranova(rm);  

% save the important statistics to individually-named variables. 
p = ranovatbl{1,5}; f = ranovatbl{1,4}; df_m = ranovatbl{1,2}; df_r = ranovatbl{2,2}; 

% run bonferroni-corrected post-hoc tests to look for relationships between the three conditions. 
posthoc = multcompare(rm, 'Amplitude_Modulation', 'ComparisonType', 'bonferroni');

% create a square NaN matrix with a length the number of conditions (3). 
P = nan(numel(plot_matrix), numel(plot_matrix)); 

% add the contrast-specific bonferroni-corrected p values to the relevant position in this NaN matrix. 
P(1,2) = posthoc{1,5}; P(1,3) = posthoc{2,5}; P(2,3) = posthoc{4,5};

% combine the important rm anova statistics and post-hoc tests.
rm_output = [p, f, df_m, df_r, posthoc{1,5}, posthoc{2,5}, posthoc{4,5}]; 

% save this data in a ROI-specific .csv file.
csvwrite_with_headers(strcat(stats_output_dir, input_ROI,'_', analysis_condition,...
    '_', toggle,'.csv'), rm_output, headers); 

%% --------------------------- PLOT THE DATA --------------------------- %% 

% if the main effect p-value was not significant, leave all remaining
% p-values as NaN (they will not be plotted). 
if p > .05
    P = nan(numel(plot_matrix), numel(plot_matrix));
  
else % otherwise, if the main effect p-value was significant:    
    
    % if any p-value in this matrix is greater than .05, set to NaN to 
    % make significance bar plotting nicer. this will also leave the
    % remaining significant post-hoc comparison p-values for plotting. 
    P(P>.05) = NaN; 
end

% mirror these specified p-values so they are matched on both sides of the
% diagonal. 
PT = P'; lidx = tril(true(size(P)), -1); P(lidx) = PT(lidx);

C = [.3 .3 .3; .5 .5 .5; .7 .7 .7]; % specify the bar colours (grayscale). 

fig = gcf;

% plot the means of the univariate betas, with corresponding standard
% errors and significance values. 
superbar(plot_matrix, 'E', stderr_matrix, 'P', P, 'BarFaceColor', C, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2);

% specify the analysis-specific y axis limits (normalised are much smaller
% as we remove the passive activation and hence the amplitude of the betas
% decrease). 
if isequal(toggle, 'normalised')
    ylim([-0.3, 2])
else
    ylim([-0.6, 2])
end

% set condition-specific x-axis labels. 
set(gca, 'xtick', 1:3, 'xticklabel', xtick_labels); 

% save the figure as a high-resolution .pdf, in high resolution and a
% specified size across all iterations. 
set(fig,'Units','Inches'); pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf','-r200', '-bestfit', strcat(figure_output_dir, input_ROI, '_', analysis_condition, '_', toggle, '.pdf'));

close(); % close the figure. 
end

%% --------------------------------------------------------------------- %%