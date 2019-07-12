%% -------- PLOT SVM BACKPROJECTION PARAMETER MAPS IN FLAT VIEW -------- %%

% Loads in the previously-created support vector parameter maps indicating
% each voxels 'preference' for a particular visual feature, and for each
% participant in turn, we load up the relevant mrvista session along with
% the relevant ROIs and transform this data to a flat window for writing to
% image files.

% for now, this scripts must be run section-by-section (ctrl+enter), with
% the operator setting the publish figure parameters, publishing the figure
% and switching hemispheres between sections.

% 10/08/2018 KWN

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

clear; clc; close all; % general housekeeping.

% specify the participants we wish to analyse. 
participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};

% specify the parameter map we wish to plot.
parametermap = 'R3773OrientationContrastbetas-V1';

exp_condition = 'original'; % specify the overall experiment we are analysing.
condition = 'block'; % specify the particular analysis we wish to investigate.

% specify the ROIs we wish to overlay on the parameter map. Currently, when
% we publish a flat figure, it joins the dorsal and ventral components of
% a 'combined' ROI (which is not the case when we examine the individual
% ROI in the flat or gray view). For now, we tackle this by plotting the
% individual dorsal and ventral ROIs (which are identical to the combined
% versions).
ROIs{1} = '01_Combined_V1_KWN';

%% ---------------- LOAD PARTICIPANT-SPECIFIC FLAT VIEW ---------------- %%

for i = 1:length(participants) % for each participant in turn:
    
    participant = participants{i}; % extract the participant number we are currently analysing. 
    
    % specify the locations of the condition- and participant-specific mrvista
    % and ROI directories.
    
    % for one participant, their MT ROI is associated with a different
    % high-resolution T1 structural scan, and so if we are processing data from
    % this participant, and analysing the MT+ ROI, we specify the corresponding
    % mrvista session and ROI directory.
    
    % so, first, we check if any of the ROIs we wish to present are the MT+
    % ROI, then if we are presenting the MT+ ROI, we change the file paths
    % accordingly.
    MT_ROI_checker = strfind(ROIs, 'MT+');
    
    if isequal(participant, 'R3111') && sum(cellfun(@isempty, MT_ROI_checker)) == length(ROIs)-1
        mrvista_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/MT', exp_condition, participant, condition);
        ROI_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs/MT/', exp_condition, participant);
        
    else % for all other instances:
        mrvista_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s', exp_condition, participant, condition);
        ROI_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs/', exp_condition, participant);
    end
    
    % specify the desired output location for the images and the
    % corresponding filename.
    output_dir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/backprojection/betas', exp_condition, condition);
    filename = sprintf('%s_%s', participant, parametermap);
    
    % index into the condition- and participant- specific mrvista directory and
    % open a gray 3-view window.
    cd(mrvista_dir);
    VOLUME{1} = mrVista('3view');
    
    % ensure we are viewing the 1st scan of the GLMs data type.
    VOLUME{1} = selectDataType(VOLUME{1},4);
    VOLUME{1} = setCurScan(VOLUME{1},1);
    VOLUME{1} = refreshScreen(VOLUME{1},0);
    
    % load in the specific parameter map of interest.
    VOLUME{1} = loadParameterMap(VOLUME{1}, strcat('/betas/difference/', parametermap, '.mat'));
    VOLUME{1}= refreshScreen(VOLUME{1},0);
    
    % for each ROI we wish to overlay, load the specific ROI into the gray view
    % window in turn.
    for i  = 1:length(ROIs)
        ROI = ROIs{i};
        VOLUME{1}=loadROI(VOLUME{1},strcat(ROI_dir,ROI),[],[],1,0);
    end
    
    % open a flat view window corresponding to the ROIs we wish to present.
    % Unless we are looking at the MT+ ROI, this flat view will be the default
    % 'Flat_KWN' view.
    openFlatWindow('Flat_KWN');
    
    % ensure again we are viewing the first scan of the GLMs data type within
    % this window.
    FLAT{1} = selectDataType(FLAT{1}, 4);
    FLAT{1} = setCurScan(FLAT{1},1);
    
    % transform the parameter map of interest from the gray to the flat view
    % window.
    FLAT{1}= vol2flatParMap(VOLUME{1},FLAT{1},1,1);
    FLAT{1}= setDisplayMode(FLAT{1},'map');
    FLAT{1}= refreshScreen(FLAT{1},0);
    
    % transform the ROIs of interest from the gray to the flat view window.
    FLAT{1}=vol2flatAllROIs(VOLUME{1},FLAT{1});
    FLAT{1} = refreshScreen(FLAT{1},0);
    
    % specify a bicolour colourmap and manually clip the scale of this map.
    FLAT{1}.ui.mapMode=setColormap(FLAT{1}.ui.mapMode, 'coolhotCmap');
    FLAT{1} = setClipMode(FLAT{1}, 'map', [1.5 -1.5]);
    FLAT{1}= refreshScreen(FLAT{1},0);

    % specify we wish to plot the ROIs on the flat view figures as coloured
    % outlines ('filled perimeter')- this is very temperamental, but the
    % below currently works. I think the word 'test' is arbitrary, but it
    % doesn't seem to work without a string here, so for now, it works, so
    % i'm not going to tamper with it anymore!
    FLAT{1} = viewSet(FLAT{1},'roidrawmethod', 'filled perimeter');
    FLAT{1} = viewSet(FLAT{1},'filledperimeter', 'test');
    FLAT{1}= refreshScreen(FLAT{1},0);
    
    % clip the colour map at -1 1, for a better range of colours in the
    % figure. 
    FLAT{1} = viewSet(FLAT{1},'mapclip', [-3 3]);
    
    % set the threshold for displaying a voxels' activation very small,
    % so that data from the majority of voxels is displayed- helping us to
    % more clearly identify patterns. the data we are displaying are not
    % p-values, so this is okay. 
    FLAT{1} = viewSet(FLAT{1},'mapwin', [.00001, -.00001]);
    FLAT{1}= refreshScreen(FLAT{1},0);
    
    % select the left hemisphere flat map and publish the figure with the
    % default parameter options. 
    selectButton(FLAT{1}.ui.sliceButtons,1); FLAT{1}=refreshScreen(FLAT{1},2);
    FLAT{1} = publishFigure(FLAT{1});
    
    % save this figure as a high-resolution .pdf image. 
    fig = gcf;
    set(fig,'Units','Inches'); pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('-dpdf','-r200', strcat(output_dir, sprintf('%s_%s_LH', participant, parametermap)));
    
    % switch to the right hemisphere flat map view, and also publish this
    % figure with the same default options. 
    selectButton(FLAT{1}.ui.sliceButtons,2); FLAT{1}=refreshScreen(FLAT{1},2);
    FLAT{1} = publishFigure(FLAT{1});
    
    % save this right-hemisphere figure with a corresponding hemisphere
    % name, as a high-resolution .pdf image. 
    fig = gcf;
    set(fig,'Units','Inches'); pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('-dpdf','-r200', strcat(output_dir, sprintf('%s_%s_RH', participant, parametermap)));
    close all
end % continue for the next participant. 
    
%% --------------------------------------------------------------------- %%