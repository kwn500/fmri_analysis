%% ----- SAVES BACKPROJECTION FIGURE WITHOUT GRIDLINES (WORKAROUND) ---- %%

% workaround to help improve figure quality from linux- opens figure and
% resaves without grid and removes the small white lines on top of the
% figure. 

% 10/7/19 KWN

%% -------------------- SAVE BACKPROJECTION FIGURE --------------------- %% 

ROI = 'V4'; % specify ROI to analyse.

% specify data storage directory.
root_dir = sprintf('/Users/kirstie/Documents/analysis/colour/final/interaction/backprojection_figures/%s/fig_files/', ROI);

% load in backprojection figure. 
fig = openfig(strcat(root_dir, ...
    sprintf('groupluminance-shape%s_surface_polar_thresh1.70_annulus.fig', ROI)));

grid off % remove gridlines. 

% remove x and y axis labels. 
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% remove colourbar. 
colorbar('off')

% save as figure. 
export_fig
close
%% --------------------------------------------------------------------- %%