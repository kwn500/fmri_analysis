%% ------------ REFORMAT BLOCK PARFILES TO 'EVENT' PARFILES ------------ %%

% loads in parfiles coding experimental block design and codes the event
% occuring within a single TR period.

% KWN 9/7/2019

clear; clc; % general housekeeping.

%%

experiment_condition = 'original'; % specify experiment to analyse.

% specify list of participants to analyse, their corresponding experimental
% numbers and the numbers of their valid fMRI runs.
% we also specify the full list of TRs we wish to specify our data across
% (this does not include interblock/cue periods).
% we also create a matrix containing the TRs of interblock periods and
% their corresponding numeric label.
if isequal(experiment_condition, 'original')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3455', 'R3517', 'R3773', 'R3932', 'R4065', 'R4127', 'R4496'};
    number = {'21', '04', '02', '25', '22', '16', '13', '01', '10', '19', '03', '17'};
    scans = {1:6, 1:6, 1:6, 1:6, [1 2 3 4 7], 1:6, 1:7, 1:6, 1:6, 1:6, 1:6, 1:6};
    
    trs = [9:3:21, 33:3:45, 57:3:69, 81:3:93, 105:3:117, 129:3:141, 153:3:165, 177:3:189, 201:3:213, 225:3:237,...
        249:3:261, 273:3:285, 297:3:309, 321:3:333, 345:3:369, 372:3:381];
    
    ib_trs = [0 24 28 72 96 120 144 168 192 216 240 264 288 312 336 360 384];
    ib_data =  [ib_trs;ones(1,length(ib_trs))*6]';
    
elseif isequal(experiment_condition, 'colour')
    participants = {'R2268', 'R2548', 'R2590', 'R2904', 'R3111', 'R3517', 'R3773', 'R4059', 'R4065', 'R4244', 'R4831', 'R4928'};
    number = {'12', '15', '04', '07', '08', '05', '01', '13', '06', '14', '11', '09'};
    scans = {1:8, 1:8, 1:7, 1:8, 1:7, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8};
    
    trs = [9:3:21, 33:3:45, 57:3:69, 81:3:93, 105:3:117, 129:3:141, 153:3:165, 177:3:189, 201:3:213, 225:3:237,...
        249:3:261, 273:3:285];
    
    ib_trs = [0 24 28 72 96 120 144 168 192 216 240 264 288];
    ib_data =  [ib_trs;ones(1,length(ib_trs))*6]';
end

%%
for thispp = 1:length(participants) % for each participant in turn:
    
    participant = participants{thispp};
    participantnum = number{thispp};
    
    for thisscan = scans{thispp} % for each valid scan in turn:
        
        % load in the run-specific parfile data and remove entries
        % corresponding to interblock periods.
        if isequal(experiment_condition, 'original')
            data = load(sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/logs/scan_0%d/vista/mrvista_event_log_pp_%s_scan_0%d_interblock&fix-5.par',...
                participant, thisscan, participantnum, thisscan));
            
            data(data(:,2)==5,:) = [];
        elseif isequal(experiment_condition, 'colour')
            data = load(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/logs/run_0%d/mrvista/mrvista_feature_event_pp_%s_run_0%d.txt',...
                participant, thisscan, participantnum, thisscan));
            
            data(data(:,2)==0,:) = [];
            data(data(:,2)==1,:) = [];
            
            % for the colour analysis, we also load in corresponding .csv
            % files which outline the stimulus-feature changes occuring at
            % each TR.
            changedata =  readtable(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/logs/run_0%d/%s_%s_run_0%d_data_edit.csv',...
                participant, thisscan, participantnum, participant, thisscan));
            
            ochange = cell2mat(table2cell(changedata(:,3))); ochange = [data(:,1), ochange];
            cchange = cell2mat(table2cell(changedata(:,4))); cchange = [data(:,1), cchange];
            schange = cell2mat(table2cell(changedata(:,5))); schange = [data(:,1), schange];
        end
        
        trcounter = 1; % initialise TR counter.
        
        for thistr = 1:length(trs) % for each TR in turn:
            
            lowtr = trs(thistr); % specify the lower TR (the TR we wish to create an event for).
            
            % if this is the last TR of the scan, specify the end point of
            % the scan as 3 TRs on from this, otherwise, select the next
            % TR in sequence for processing.
            if thistr == length(trs)
                hightr = trs(end) + 3;
            else
                hightr = trs(thistr+1);
            end
            
            % find the indicies corresponding to the TR period to analyse.
            if isequal(experiment_condition, 'original')
                ind = find(data(:,1) >= lowtr & data(:,1) <hightr);
            elseif isequal(experiment_condition, 'colour')
                ind = find(ochange(:,1) >= lowtr & ochange(:,1) <hightr);
            end
            
            for thisind = 1:length(ind) % for this extracted TR data:
                
                % extract the event numbers occuring within this event.
                if isequal(experiment_condition, 'original')
                    vals =  data(ind,2);
                elseif isequal(experiment_condition, 'colour')
                    ovals = ochange(ind,2);
                    cvals = cchange(ind,2);
                    svals = schange(ind,2);
                    
                    vals = [sum(ovals), sum(cvals), sum(svals)];
                end
                
                
                if isequal(experiment_condition, 'original')
                    if length(vals) == 2
                        if sum(ismember(vals,[1 1])) == 2 % if both events were orientation:
                            
                            % code this event as an orientation event for this 3s TR.
                            finaldata(trcounter,:) = [lowtr, 1]; 
                            
                            % repeat this process for all other combinations of event.
                        elseif sum(ismember(vals,[2 2])) == 2 % contrast.
                            finaldata(trcounter,:) = [lowtr, 2];
                        elseif sum(ismember(vals,[3 3])) == 2 % shape.
                            finaldata(trcounter,:) = [lowtr, 3];
                        elseif sum(ismember(vals,[4 4])) == 2 % no change.
                            finaldata(trcounter,:) = [lowtr,4];
                        elseif sum(ismember(vals,[1 4])) == 2 % orientation.
                            finaldata(trcounter,:) = [lowtr, 1];
                        elseif sum(ismember(vals,[2 4])) == 2 % contrast.
                            finaldata(trcounter,:) = [lowtr,2];
                        elseif sum(ismember(vals,[3 4])) == 2 % shape.
                            finaldata(trcounter,:) = [lowtr,3];
                        else
                            % mismatch (two different no change events occuring in a TR).
                            finaldata(trcounter,:) = [lowtr,6]; 
                        end
                    elseif length(vals) >= 3 % if three events occured within this TR:
                        finaldata(trcounter,:) = [lowtr,6]; % code this event as a mismatch.
                    end
                    
                elseif isequal(experiment_condition, 'colour')
                    % perform the same process as above, coding events in
                    % respect to the two events occuring within a single
                    % TR. 
                    if sum(vals) == 0
                        finaldata(trcounter,:) = [lowtr, 4]; % no change
                    elseif vals(1) == 1 && vals(2) == 0 && vals(3) == 0
                        finaldata(trcounter,:) = [lowtr, 1]; % orientation
                    elseif vals(2) == 1 && vals(1) == 0 && vals(3) == 0
                        finaldata(trcounter,:) = [lowtr, 2]; % contrast
                    elseif vals(3) == 1 && vals(1) == 0 && vals(2) == 0
                        finaldata(trcounter,:) = [lowtr, 3]; % shape
                    else
                        finaldata(trcounter,:) = [lowtr,5]; % multiple/mismatch
                    end
                end
                
                trcounter = trcounter + 1;
            end
            
            % change the numeric label corresponding to interblock events.
            if isequal(experiment_condition, 'original')
                indchange = find(finaldata(:,2) == 6);
                finaldata(indchange,2) = 5;
            end
            
            % combine the recoded event and interblock data.
            finaldata = [finaldata;ibdata];
            
            % sort the data into chronological order. 
            finaldata = sortrows(finaldata);
            
            % write this data to a run-specific .par file.
            if isequal(experiment_condition, 'original')
                dlmwrite(sprintf('/scratch/groups/Projects/P1323/original/fMRI/%s/mrvista/event/Stimuli/Parfiles/mrvista_event_log_scan_0%d_3tr-new.par',...
                    participant, thisscan), finaldata, 'delimiter', '\t');
            elseif isequal(experiment_condition, 'colour')
                dlmwrite(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/%s/logs/run_0%d/mrvista/mrvista_feature_event_%s_run_0%d_edited.par',...
                    participant, thisscan, participant, thisscan), finaldata, 'delimiter', '\t');
            end
            
            clear finaldata
        end
    end
end
  
%% --------------------------------------------------------------------- %%