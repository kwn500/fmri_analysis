function extractpRFdata(prf_dir, exp_condition, participantnumber, participant, ROIanalysis, outputROI)

% index into the participant-specific pRF mrvista directory.
prf_directory = prf_dir{participantnumber};

cd(prf_directory); % index into this mrvista pRF session directory.

% there is one participant, which (for an unknown reason), means we
% have to manually specify the anatomy everytime, so here, for this
% participant we open a mrvista gray 3-view window with the
% corresponding high-resolution structual scan specified.
if isequal(participant, 'R3111')
    anat = '/groups/labs/wadelab/data/pRF/Himmel/TemporalContrastpRF/pRF/anatomies/R3111/vAnatomy.nii.gz';
    VOLUME{1} = mrVista('3view', anat);
else
    % for all other participants, we open a mrvista gray view window
    % and it loads each participants' corresponding high-resolution T1
    % anatomical scan.
    VOLUME{1} = mrVista('3view');
end

% some participants have previously collected pRF data from another PhD
% student, if we are analysing these participants, they have multiple
% analyses contained within their sessions, so on a
% participant-by-participant basis, we ensure we select the 'average
% across all conditions' datatype.
if ~isempty(strfind(prf_directory, 'CombinedSessions'))
    if ~isequal(participant,'R3111') && ~isequal(participant, 'R3455')
        VOLUME{1}.curDataType = 3; VOLUME{1}=refreshScreen(VOLUME{1});
    elseif isequal(participant, 'R3111')
        VOLUME{1}.curDataType = 5; VOLUME{1}=refreshScreen(VOLUME{1});
    elseif isequal(participant, 'R3455')
        VOLUME{1}.curDataType = 6; VOLUME{1}=refreshScreen(VOLUME{1});
    end
    retinotopic_model = strcat(prf_directory, '/Gray/Average_AllConditions/');
else
    % otherwise, for all other participants, select the default average
    % retinotopic model.
    VOLUME{1}.curDataType = 4; VOLUME{1}=refreshScreen(VOLUME{1});
    retinotopic_model = strcat(prf_directory,'/Gray/Averages/');
end
model = dir(strcat(retinotopic_model, '*fFit*'));

% load this retinotopic model into the mrvista session.
VOLUME{1} = rmSelect(VOLUME{1}, 2, model.name);
VOLUME{1} = rmLoadDefault(VOLUME{1});

% specify the participant-specific ROI storage directory.
if isequal(participant, 'R3111') && isequal(outputROI, 'MT+')
    ROI_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/MT/ROIs/', exp_condition, participant);
else
    ROI_directory = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/anatomy/ROIs/', exp_condition, participant);
end

% specify the full path to the ROI of interest and load this ROI into
% the mrvista window.
ROI = strcat(ROI_directory, ROIanalysis);
VOLUME{1}=loadROI(VOLUME{1},ROI,[],[],1,0);
VOLUME{1}=refreshScreen(VOLUME{1},0);

% extract the coordinates for this ROI (in terms of their location in
% the whole brain data).
ROIcoords = viewGet(VOLUME{1}, 'roiindices');

% extract the retinotopic model data.
model = viewGet(VOLUME{1}, 'rmmodel'); model = model{1};

% from this retinotopic model, extract each voxels' eccentricity and
% polar angle data (for all voxels across the brain).
ROIecc = rmGetVoxelData('ecc', VOLUME{1});
ROIpol= rmGetVoxelData('pol', VOLUME{1});
ROIx = rmGetVoxelData('x', VOLUME{1});
ROIy = rmGetVoxelData('y', VOLUME{1});

% specify and create a pRF data output directory within the
% participant-specific ROI directory.
output_directory = strcat(ROI_directory, 'pRFdata/'); [~,~] = mkdir(output_directory);

% save this extracted ROI-specific eccentricity, polar angle and
% co-ordinate data within this specified directory.
save(strcat(output_directory, sprintf('%s%s_pRF_data.mat', participant, outputROI)), 'ROIecc', 'ROIpol', 'ROIcoords', 'ROIx', 'ROIy');

close all % close all mrvista windows.

end
%% --------------------------------------------------------------------- %%