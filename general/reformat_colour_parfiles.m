%% --------------------- RE-FORMAT COLOUR PARFILES --------------------- %% 

% reads in colour fMRI parfiles, and re-codes into chromatic attention (as
% oppose to orientation, contrast and shape-specific attention conditions)
% and chromatic passive conditions. 

% KWN 9/7/2019

clear; clc;

%%

participant = 'R4928'; % specify participant. 

% extract colour fMRI parfiles for each run. 
files = dir(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/mrvista/colourxfeature/Stimuli/Parfiles/*.par', participant));
filenames = {files.name};

% for each run parfile in turn:
for thisFile = 1:length(filenames)
    
    % load in this parfile information (timing and event code). 
    fid = fopen(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/mrvista/colourxfeature/Stimuli/Parfiles/%s',...
        participant,filenames{thisFile}), 'rt');
    C = textscan(fid, '%d %d');
    
    % extract the timing and coding information as seperate variables. 
    TRs = double(C{1});
    codes = double(C{2});
    
    for thisTR = 1:size(TRs,1) % for each TR in turn:
        
        % if this TR is an orientation, contrast or shape condition, 
        % combine this into a singular chromatic attention condition. 
        if ismember(codes(thisTR), [1 2 3])
            newcodes(thisTR) = 1;
            labels{thisTR} = 'RGAttention';
        elseif ismember(codes(thisTR), [5 6 7])
            newcodes(thisTR) = 3;
            labels{thisTR} = 'BYAttention';
        elseif ismember(codes(thisTR), [9 10 11])
            newcodes(thisTR) = 5;
            labels{thisTR} = 'LumAttention';
        elseif codes(thisTR) == 4
            newcodes(thisTR) = 2;
            labels{thisTR} = 'RGPassive';
        elseif codes(thisTR) == 8
            newcodes(thisTR) = 4;
            labels{thisTR} = 'BYPassive';
        elseif codes(thisTR) == 12
            newcodes(thisTR) = 6;
            labels{thisTR} = 'LumPassive';
        elseif codes(thisTR) == 13
            newcodes(thisTR) = 7;
            labels{thisTR} = 'Interblock';
        end
    end
         
    % combine the timing and coding data as well as string condition
    % descriptors to form the new parfile structure. 
    parfiledata = [num2cell(TRs), num2cell(newcodes'), labels'];
    
    % write this parfile to a run-specific .par file. 
    outputdir = sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/mrvista/colour-passive/Stimuli/Parfiles/', participant); 
    outputfile = strcat(outputdir, 'colour-passive',filenames{thisFile}(15:end));
    dlmcell(outputfile, parfiledata)
end

%% --------------------------------------------------------------------- %%