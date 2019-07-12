function makefigure(plotmatrix,stderrmatrix,p,yscale,xtickvals,xticklabels,outputdir,name,posthoc)

% takes in data to plot, error bars, significance values, figure formatting
% values and an output directory and plots the data as a bar chart, before
% writing to a .pdf file. 

% 11/7/19 KWN

%% ----------------------------- PLOT DATA ----------------------------- %%

% create significance value matrix with post-hoc values provided there is a
% significant main effect identified.
P = nan(numel(plotmatrix), numel(plotmatrix));
if p < .05
    if size(plotmatrix,2) == 2
        P(1,2) = p; 
    elseif size(plotmatrix,2) == 3
        P(1,2) = posthoc{1,5}; P(1,3) = posthoc{2,5}; P(2,3) = posthoc{4,5};
    end
end
PT = P'; lidx = tril(true(size(P)), -1); P(lidx) = PT(lidx);
 
% specify the bar colour within the matrix. 
if size(plotmatrix,2) == 2
    colour = [.3 .3 .3; .5 .5 .5];
elseif size(plotmatrix,2) == 3
    colour =[.3 .3 .3; .5 .5 .5; .7 .7 .7];
end

% plot the data with error bars and associated significance values. 
superbar(plotmatrix, 'E', stderrmatrix,...
    'P', P,'BarFaceColor', colour, 'BarEdgeColor', [0 0 0], 'BarLineWidth', 1,...
    'ErrorbarRelativeWidth', 0.25, 'ErrorbarLineWidth', 1.5, 'PStarColor', [0 0 0],...
    'PStarFontSize', 12, 'PLineColor', [0 0 0], 'PLineWidth', 1.2,'PStarShowNS', false);
ylim(yscale)

% save the data to a .pdf file. 
set(gca, 'xtick', xtickvals, 'xticklabel', xticklabels); 
saveas(gcf, strcat(outputdir,name), 'pdf');
end

%% --------------------------------------------------------------------- %%