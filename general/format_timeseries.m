%% --------- REFORMAT PARTICIPANT TIMESERIES (TRS X CONDITION) --------- %%

% For each ROI and participant, loads in the participant's timeseries and
% parfiles. We reformat the parfiles so that we have an attentional
% condition for every TR (as oppose to blocks). We then extract the
% timeseries across voxels for each TR for each feature. These modified
% timeseries are saved as a struct.
% we can also take the timeseries for the top 100 voxels (in terms of
% percentage variance explaned) for each participant and GLM.

% 05/08/2018 KWN

clear; clc; % general housekeeping.

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

exp_condition = 'original'; % specify experimental data to analyse.
condition = 'block'; % specify type of analysis to perform.

% specify the list of participant numbers to analyse.
if isequal(exp_condition, 'original');
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455',...
        'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};
elseif isequal(exp_condition, 'colour');
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3517',...
        'R3773', 'R4059', 'R4065', 'R4244', 'R4831', 'R4928'};
elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
    participants = {'R2548', 'R2904', 'R3111', 'R3517', 'R4059', 'R4065',...
        'R4127', 'R4244', 'R4831', 'R4833', 'R4890', 'R4928'}';
elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3x3feature')
    participants = {'R2548', 'R2904', 'R3111', 'R3517', 'R3773', 'R4059',...
        'R4065', 'R4127', 'R4244','R4831', 'R4833', 'R4890', 'R4928', 'R5006'};
end

% specify the ROIs to analyse.
ROIs = {'V1','V2', 'V3', 'V3AB', 'V4', 'LO1', 'LO2', 'MT+', 'A1', 'IPS0'};

% specify a processed timeseries output directory.
if ~isequal(exp_condition, 'naturalistic')
    outputdir = (sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/processed_timeseries/',...
        exp_condition, condition));
elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
    outputdir = ('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/processed_timeseries/');
elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3x3feature')
    outputdir = ('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3x3feature/processed_timeseries/');
end

voxels = 50; % specify the number of voxels we wish to extract.

%% ------------------------- PROCESS TIMESERIES ------------------------ %%

for ROInum = 1:length(ROIs) % for each ROI in turn:
    ROI = ROIs{ROInum}; % extract the ROI-specific name.
    
    for ppnum = 1:length(participants) % for each participant in turn:
        participant = participants{ppnum}; % extract the participant number.
        
        % load the participant- and ROI- specific timeseries data.
        if ~isequal(exp_condition, 'naturalistic')
            load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/timeseries/%s/%s_tSeries_%s.mat',...
                exp_condition, condition, ROI, participant, ROI));
            load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/vista_output/multivariate/%s/varexp/%s/%s_varexp_%s.mat',...
                exp_condition, condition, ROI, participant, ROI));
        elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
            load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/timeseries/%s/%s_tSeries_%s.mat',...
                ROI, participant, ROI));
            load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/varexp/%s/%s_varexp_%s.mat',...
                ROI, participant, ROI));
        elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3x3feature')
            load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3x3feature/timeseries/%s/%s_tSeries_%s.mat',...
                ROI, participant, ROI));
            load(sprintf('/scratch/groups/Projects/P1361/fMRI/vista_output/multivariate/3x3feature/varexp/%s/%s_varexp_%s.mat',...
                ROI, participant, ROI));
        end
        
        % sort the percentage variance each voxel explains of the GLM in
        % order from highest to lowest.
        [sorted_varexp_vals, sorted_varexp_ind] = sort(varexp,'descend');
        
        % record the values of the top 100 voxels in terms of percentage
        % variance explained.
        sorted_varexp_values = sorted_varexp_vals(1:voxels);
        
        % extract the indicies of the top 100 voxels in the original
        % (unsorted) array.
        sorted_varexp_indicies = sorted_varexp_ind(1:voxels);
        
        % extract only the timeseries data for these top 100 voxels for
        % further processing.
        extractedData = extractedData(:,sorted_varexp_indicies);
        
        % load the participant-specific parfile directory.
        if ~isequal(exp_condition, 'naturalistic')
            parfiledir = (sprintf('/scratch/groups/Projects/P1323/%s/fMRI/%s/mrvista/%s/Stimuli/Parfiles/',...
                exp_condition, participant, condition));
            parfiles = dir(strcat(parfiledir, '*.par')); parfiles = {parfiles.name};
        else
            if isequal(participant, 'R3517');
                if isequal(condition, '3feature')
                    parfiledir = (sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/session1-allscans/Stimuli/Parfiles/', participant));
                else
                    parfiledir = (sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/session1-allscans/multiple/Stimuli/Parfiles/', participant));
                end
            else
                if isequal(condition, '3feature')
                    parfiledir = (sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/Stimuli/Parfiles/', participant));
                else
                    parfiledir = (sprintf('/scratch/groups/Projects/P1361/fMRI/%s/mrvista/multiple/Stimuli/Parfiles/', participant));
                end
            end
            parfiles = dir(strcat(parfiledir, '*.par')); parfiles = {parfiles.name};
        end
        
        counter = 1; % establish a counter variable to track the number of parfiles processed.
        data = []; % create an empty data storage array.
        
        for i = 1:length(parfiles) % for each parfile in turn:
            
            % load in the parfile data.
            [timing, event, labels] = textread(strcat(parfiledir, parfiles{i}),'%f %d %s');
            pdata{i} = [timing,event];
            
            if i > 1 % if this is not the first parfile we have processed:
                
                % multiply the TRs by number of TRs x number of scans, to
                % keep a running total of TRs across all scans to match the
                % formatting of the timeseries data, which is concatenated
                % across scans.
                if isequal(exp_condition, 'original');
                    pdata{i}(:,1) = pdata{i}(:,1) + 393 * counter;
                elseif isequal(exp_condition, 'colour');
                    pdata{i}(:,1) = pdata{i}(:,1) + 297 * counter;
                elseif isequal(exp_condition, 'naturalistic')
                    pdata{i}(:,1) = pdata{i}(:,1) + 300 * counter;
                end
                counter = counter+1; % increment the parfile counter variable.
            end
            % add this edited parfile data to the end of the two-column
            % array of parfile timing and event data.
            data = [data; pdata{i}];
        end % continue to the next parfile.
        
        % specify a list of TRs (in 3s) matching the length of the TRs
        % across all scans.
        if ~isequal(exp_condition, 'naturalistic')
            fullTRs = [0:3:data(end,1)+8]';
        else
            fullTRs = [0:2:data(end,1)+29]';
        end
        
        % specify a column matrix of zeros matching the size of the number
        % of TRs across all scans.
        events = zeros(length(fullTRs),1);
        
        % initialise a counter variable to track the number of TRs
        % processed matching the parfiles.
        datacounter = 1;
        
        for i = 1:length(fullTRs) % for each of the full-scan TRs in turn:
            
            % try to extract the index of the TR in the original parfile
            % data matching the TR in the full timing data.
            ind = find(fullTRs(i,1) == data(:,1));
            
            if ~isempty(ind) % if this index is not empty:
                if i ~= 1 % if this is not the first index we have extracted:
                    datacounter = datacounter+1; % increment the TR counter variable.
                end
                % add the attention event data to the relevant position in
                % the full TR parfile data (this will be a new event).
                fullTRs(i,2) = data(datacounter,2);
                
            elseif isempty(ind) % otherwise, if this index is empty:
                % add the attention event data to the relevant position in
                % the full TR parfile data (this will be the same as the
                % last index identified).
                fullTRs(i,2) = data(datacounter,2);
            end
        end % continue to the next TR.
        
        % extract the indicies of the feature-specific TRs in this new
        % formatted full TR dataset, and use these indicies to extract the
        % feature-specific TR data across all voxels.
        if isequal(exp_condition, 'original') || isequal(condition, 'feature');
            
            oind = find(fullTRs(:,2) == 1); orientation = extractedData(oind,:);
            cind = find(fullTRs(:,2) == 2); contrast = extractedData(cind,:);
            sind = find(fullTRs(:,2) == 3); shape = extractedData(sind,:);
            pind = find(fullTRs(:,2) == 4); passive = extractedData(pind,:);
            ibind = find(fullTRs(:,2) == 5); interblock = extractedData(ibind,:);
            
            % save these new feature-specific timeseries x voxel data matricies
            % in a single struct.
            tSeries.orientation = orientation;
            tSeries.contrast = contrast;
            tSeries.shape = shape;
            tSeries.passive = passive;
            tSeries.interblock = interblock;
            
        elseif isequal(exp_condition, 'colour') && isequal(condition, 'colour');
            
            % perform the same process for the red-green, blue-yellow and
            % luminance data.
            rgind = find(fullTRs(:,2) == 1); rg = extractedData(rgind,:);
            byind = find(fullTRs(:,2) == 2); by = extractedData(byind,:);
            lumind = find(fullTRs(:,2) == 3); lum = extractedData(lumind,:);
            
            tSeries.rg = rg;
            tSeries.by = by;
            tSeries.lum = lum;
            
        elseif isequal(exp_condition, 'colour') && isequal(condition, 'colourxfeature');
            
            rgoind = find(fullTRs(:,2) == 1); rgo = extractedData(rgoind,:);
            rgcind = find(fullTRs(:,2) == 2); rgc = extractedData(rgcind,:);
            rgsind = find(fullTRs(:,2) == 3); rgs = extractedData(rgsind,:);
            rgpind = find(fullTRs(:,2) == 4); rgp = extractedData(rgpind,:);
            byoind = find(fullTRs(:,2) == 5); byo = extractedData(byoind,:);
            bycind = find(fullTRs(:,2) == 6); byc = extractedData(bycind,:);
            bysind = find(fullTRs(:,2) == 7); bys = extractedData(bysind,:);
            bypind = find(fullTRs(:,2) == 8); byp = extractedData(bypind,:);
            lumoind = find(fullTRs(:,2) == 9); lumo = extractedData(lumoind,:);
            lumcind = find(fullTRs(:,2) == 10); lumc = extractedData(lumcind,:);
            lumsind = find(fullTRs(:,2) == 11); lums = extractedData(lumsind,:);
            lumpind = find(fullTRs(:,2) == 12); lump = extractedData(lumpind,:);
            
            tSeries.rgo = rgo; tSeries.rgc = rgc; tSeries.rgs = rgs; tSeries.rgp = rgp;
            tSeries.byo = byo; tSeries.byc = byc; tSeries.bys = bys; tSeries.byp = byp;
            tSeries.lumo = lumo; tSeries.lumc = lumc; tSeries.lums = lums; tSeries.lump = lump;
            
        elseif isequal(exp_condition, 'colour') && isequal(condition, 'colour-passive')
            
            rgaind = find(fullTRs(:,2) == 1); rga = extractedData(rgaind,:);
            rgpind = find(fullTRs(:,2) == 2); rgp = extractedData(rgpind,:);
            byaind = find(fullTRs(:,2) == 3); bya = extractedData(byaind,:);
            bypind = find(fullTRs(:,2) == 4); byp = extractedData(bypind,:);
            lumaind = find(fullTRs(:,2) == 5); luma = extractedData(lumaind,:);
            lumpind = find(fullTRs(:,2) == 6); lump = extractedData(lumpind,:);
            
            tSeries.rga = rga;
            tSeries.rgp = rgp;
            tSeries.bya = bya;
            tSeries.byp = byp;
            tSeries.luma = luma;
            tSeries.lump = lump;
            
        elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3feature')
            
            faceind = find(fullTRs(:,2) == 1); face = extractedData(faceind,:);
            oind = find(fullTRs(:,2) == 2); orientation = extractedData(oind,:);
            cind = find(fullTRs(:,2) == 3); colour = extractedData(cind,:);
            sind = find(fullTRs(:,2) == 4); shape = extractedData(sind,:);
            pind = find(fullTRs(:,2) == 5); passive = extractedData(pind,:);
            
            tSeries.face = face;
            tSeries.orientation = orientation;
            tSeries.colour = colour;
            tSeries.shape = shape;
            tSeries.passive = passive;
            
        elseif isequal(exp_condition, 'naturalistic') && isequal(condition, '3x3feature')
            
             pind = find(fullTRs(:,2) == 1); passive = extractedData(pind,:);
            vind = find(fullTRs(:,2) == 2); vertical = extractedData(vind,:);
            hind = find(fullTRs(:,2) == 3); horizontal = extractedData(hind,:);
            dind = find(fullTRs(:,2) == 4); diagonal = extractedData(dind,:);
            rind = find(fullTRs(:,2) == 5); red = extractedData(rind,:);
            gind = find(fullTRs(:,2) == 6); green = extractedData(gind,:);
            bind = find(fullTRs(:,2) == 7); blue = extractedData(bind,:);
            cind = find(fullTRs(:,2) == 8); circular = extractedData(cind,:);
            sind = find(fullTRs(:,2) == 9); square = extractedData(sind,:);
            tind = find(fullTRs(:,2) == 10); triangular = extractedData(tind,:);
            faceind = find(fullTRs(:,2) == 11); face = extractedData(faceind,:);

            tSeries.passive = passive;
            tSeries.vertical = vertical;
            tSeries.horizontal = horizontal;
            tSeries.diagonal = diagonal;
            tSeries.red = red;
            tSeries.green = green;
            tSeries.blue = blue;
            tSeries.circular = circular;
            tSeries.square = square;
            tSeries.triangular = triangular;
            tSeries.face = face;
        end
        
        % create an ROI-specific directory in the output
        % processed-timeseries directory.
        ROIdir = strcat(outputdir, sprintf('%s/', ROI)); [~,~] = mkdir(ROIdir);
        
        % save this feature-specific timeseries struct with a participant-
        % and ROI- specific filename in the ROI-specific output file
        % directory.
        save(strcat(ROIdir, sprintf('%s_%s_timeseries', participant, ROI)), 'tSeries');
    end % continue for the next participant.
end % continue for the next ROI.

%% --------------------------------------------------------------------- %%