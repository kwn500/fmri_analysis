function rmanova_plot_naturalistic(input_matrix, stderr_matrix, toggle, stats_output_dir, figure_output_dir, input_ROI)

% performs one-way repeated measures ANOVA with associated post-hoc tests
% and plot the results as a bar-chart with error bars and associated
% significance values before saving to a .pdf figure. 

% 11/7/19 KWN

%% ---------------- RUN ONE-WAY REPEATED MEASURES ANOVA ---------------- %%

plot_matrix = input_matrix(:,1:end-1); stderr_matrix = stderr_matrix(:,1:end-1);

% calculate the overall mean and standard error across all participants for each of the conditions.
plot_matrix = mean(plot_matrix); stderr_matrix = std(stderr_matrix) / sqrt(length(stderr_matrix));

% create a column the length of the number of participants, containing the 'participant' string.
pp_number = repmat({'participant'}, size(input_matrix,1),1);

% create the data table, with participant, orientation, contrast and shape columns named accordingly.
t = table(pp_number, input_matrix(:,1), input_matrix(:,2), input_matrix(:,3), input_matrix(:,4),...
    'VariableNames', {'participant', 'Faces', 'Orientation', 'Contrast', 'Shape'});

% specify an overall name for the orientation, contrast and shape columns.
Meas = table([1 2 3 4]', 'VariableNames', {'Amplitude_Modulation'});

% fit a repeated measures model to this data.
rm = fitrm(t, 'Faces-Shape~1', 'WithinDesign', Meas);

% specify the corresponding x axis labels for later plotting.
xtick_labels = {'Faces', 'Orientation', 'Contrast', 'Shape'};

% specify corresponding bonferroni post-hoc comparison names for later
% writing to a .csv file header.
headers = {'rm_p', 'rm_f', 'rmdf_m', 'rmdf_r', 'BF-O', 'BF-C', 'BF-S', 'BO-C', 'BO-S', 'BC-S'};

% run the one-way repeated measures anova.
ranovatbl = ranova(rm);

% save the important statistics to individually-named variables.
p = ranovatbl{1,5}; f = ranovatbl{1,4}; df_m = ranovatbl{1,2}; df_r = ranovatbl{2,2};

% run bonferroni-corrected post-hoc tests to look for relationships between the three conditions.
posthoc = multcompare(rm, 'Amplitude_Modulation', 'ComparisonType', 'bonferroni');

% create a square NaN matrix with a length the number of conditions (3).
P = nan(numel(plot_matrix), numel(plot_matrix));

% add the contrast-specific bonferroni-corrected p values to the relevant position in this NaN matrix.
P(1,2) = posthoc{1,5}; P(1,3) = posthoc{2,5}; P(1,4) = posthoc{3,5};
P(2,3) = posthoc{5,5}; P(2,4) = posthoc{6,5};
P(3,4) = posthoc{9,5};
%
% % combine the important rm anova statistics and post-hoc tests.
rm_output = [p, f, df_m, df_r, posthoc{1,5}, posthoc{2,5}, posthoc{3,5}, posthoc{5,5}, posthoc{6,5}, posthoc{9,5}];

% save this data in a ROI-specific .csv file.
csvwrite_with_headers(strcat(stats_output_dir, input_ROI,'_',...
    toggle,'.csv'), rm_output, headers);

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
    ylim([-1, 2])
else
    ylim([-1, 2])
end

% set condition-specific x-axis labels.
set(gca, 'xtick', 1:4, 'xticklabel', xtick_labels);

% save the figure as a high-resolution .pdf, in high resolution and a
% specified size across all iterations.
set(fig,'Units','Inches'); pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf','-r200', '-bestfit', strcat(figure_output_dir, input_ROI, '_', toggle, '.pdf'));

close(); % close the figure.
end

%% --------------------------------------------------------------------- %%