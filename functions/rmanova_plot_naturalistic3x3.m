function [p,f, df_m, df_r, means, stderr, mauchlyp] = rmanova_plot_naturalistic3x3(input_matrix, stats_output_dir, figure_output_dir, input_ROI, toggle, adj_p)

% performs one-way repeated measures ANOVA with associated post-hoc tests
% and plot the results as a bar-chart with error bars and associated
% significance values before saving to a .pdf figure. 

% 11/7/19 KWN

%% ---------------- RUN ONE-WAY REPEATED MEASURES ANOVA ---------------- %%

% calculate the overall mean and standard error across all participants for each of the conditions.
plot_matrix = mean(input_matrix(:,[1 2 3]));
analysis_matrix = input_matrix(:,[1 2 3]);
stderr_matrix = std(input_matrix(:,[1 2 3])) / sqrt(length(input_matrix));

means = plot_matrix;
stderr = stderr_matrix;

if isequal(toggle, 'orientation')
    condlabels = {'Vertical', 'Horizontal', 'Diagonal'};
    posthoclabels = {'B-VH', 'B-VD', 'B-HD'};
    outputlabels = {'Orientation3x3feature'};
elseif isequal(toggle, 'colour')
    condlabels = {'Red', 'Green', 'Blue'};
    posthoclabels = {'B-RG', 'B-RB', 'B-GB'};
    outputlabels = {'Colour3x3feature'};
elseif isequal(toggle, 'shape')
    condlabels = {'Circular', 'Square', 'Triangular'};
    posthoclabels = {'B-CS', 'B-CT', 'B-ST'};
    outputlabels = {'Shape3x3feature'};
end

% create a column the length of the number of participants, containing the 'participant' string.
pp_number = repmat({'participant'}, size(input_matrix,1),1);

% create the data table, with participant, orientation, contrast and shape columns named accordingly.
t = table(pp_number, analysis_matrix(:,1), analysis_matrix(:,2), analysis_matrix(:,3),...
    'VariableNames', ['participant', condlabels]);

% specify an overall name for the orientation, contrast and shape columns.
Meas = table([1 2 3]', 'VariableNames', {'Amplitude_Modulation'});

% fit a repeated measures model to this data.
rm = fitrm(t,char(strcat(condlabels(1),'-',condlabels(3),'~1')), 'WithinDesign', Meas);

% specify the corresponding x axis labels for later plotting.
xtick_labels = condlabels;

% specify corresponding bonferroni post-hoc comparison names for later
% writing to a .csv file header.
headers = ['rm_p', 'rm_f', 'rmdf_m', 'rmdf_r', posthoclabels];

% run the one-way repeated measures anova.
ranovatbl = ranova(rm);

mauchlytbl = mauchly(rm);
mauchlyp = mauchlytbl{1,4};

if mauchlyp<.05
    
    % save the important statistics to individually-named variables.
    p = ranovatbl{1,6}; f = ranovatbl{1,4}; df_m = ranovatbl{1,2}; df_r= ranovatbl{2,2};
else
    p = ranovatbl{1,5}; f = ranovatbl{1,4}; df_m = ranovatbl{1,2}; df_r = ranovatbl{2,2};
end

% run bonferroni-corrected post-hoc tests to look for relationships between the three conditions.
posthoc = multcompare(rm, 'Amplitude_Modulation', 'ComparisonType', 'bonferroni');

% create a square NaN matrix with a length the number of conditions (3).
P = nan(numel(plot_matrix), numel(plot_matrix));

% add the contrast-specific bonferroni-corrected p values to the relevant position in this NaN matrix.
P(1,2) = posthoc{1,5}; P(1,3) = posthoc{2,5}; P(2,3) = posthoc{4,5};

% combine the important rm anova statistics and post-hoc tests.
rm_output = [p, f, df_m, df_r, posthoc{1,5}, posthoc{2,5}, posthoc{4,5}];

% save this data in a ROI-specific .csv file.
csvwrite_with_headers(char(strcat(stats_output_dir, input_ROI,'_',outputlabels,'.csv')), rm_output, headers);

%% --------------------------- PLOT THE DATA --------------------------- %% 

% if the main effect p-value was not significant, leave all remaining
% p-values as NaN (they will not be plotted).
if adj_p > .05
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
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2, 'PStarShowGT', false);

ylim([-1, 2])

% set condition-specific x-axis labels.
set(gca, 'xtick', 1:3, 'xticklabel', xtick_labels);

% save the figure as a high-resolution .pdf, in high resolution and a
% specified size across all iterations.
set(fig,'Units','Inches'); pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf','-r200', '-bestfit', char(strcat(figure_output_dir, input_ROI, '_', outputlabels, '.pdf')));

close(); % close the figure.

%% --------------------------------------------------------------------- %%
