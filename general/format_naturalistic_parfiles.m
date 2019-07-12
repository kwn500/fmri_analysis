%% ------------------- FORMAT NATURALISTIC PARFILES -------------------- %%

% loads in run-specific parfiles from the naturalistic experiment. we
% remove the data corresponding to the two-second cue periods and recoded
% the numeric codes corresponding to the face and passive conditions. 

% KWN 9/7/2019

clear; clc;

%% 

% specify experimental condition to analyse.
experiment_condition = '3feature'; 

participant = 'R5006'; % specify participant to analyse. 

% specify participant-specific parfile directory.
if isequal(experiment_condition, '3feature')
    parfiles = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/logs/3featureblock/original/*.par', participant));
elseif isequal(experiment_condition, '3x3feature')
    parfiles = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/logs/3x3featureblock/original/*.par', participant));
end

parfile_names = {parfiles.name}; % extract parfile names. 

for i = 1:length(parfile_names) % for each parfile in turn:
    
    % load in run-specific parfile data.
    if isequal(experiment_condition, '3feature')
        fid = fopen(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/logs/3featureblock/original/%s',...
            participant, parfile_names{i}), 'rt');
    elseif isequal(experiment_condition, '3x3feature')
        fid = fopen(sprintf('/scratch/groups/Projects/P1361/fMRI/%s/logs/3x3featureblock/original/%s',...
            participant, parfile_names{i}), 'rt');
    end
    
    % extract timing, numeric code and string-label information respectively. 
    C = textscan(fid, '%d %d %s');
    timings = C{1};
    codes = C{2};
    labels = C{3};
    
    % find and remove the data corresponding to cue periods. 
    cues = find(codes==1);
    codes(cues) = [];
    labels(cues) = [];
    
    % re-label conditions corresponding to face and passive conditions.
    if isequal(experiment_condition, '3feature')
        codes(codes == 5) = 1;
        codes(codes == 0) = 5;
    elseif isequal(experiment_condition, '3x3feature')
        codes(codes == 0) = 1;
    end
    
    % specify the new block TRs (with no cue periods). 
    timings = 0:30:270;
    
    % combine these new timings with the re-labelled block data. 
    newdata = [num2cell(timings)' num2cell(double(codes)) labels];
    
    % write this new parfile matrix to a run-specific .par file. 
    outputdir = sprintf('/scratch/groups/Projects/P1361/fMRI/%s/logs/%sblock/edited/',...
        participant, experiment_condition);
    [~,~]= mkdir(outputdir);
    outputfile = strcat(outputdir, parfile_names{i}(:,1:end-4), '_edited.par');
    dlmcell(outputfile, newdata)
end

%% --------------------------------------------------------------------- %%