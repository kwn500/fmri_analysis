%% --------------------- CREATE CONN DESIGN FILES ---------------------- %%

% loads in mr-vista parfiles and restructures timing and condition
% information to produce conn-toolbox compatible design files.

% 10/7/2019 KWN

clear; clc;
addpath('/scratch/sg3/P1323/code/fmri_analysis/functions');

%% -------------------- SPECIFY ANALYSIS PARAMETERS -------------------- %%

exp_condition = 'original'; % specify overall experiment to analyse.
condition = 'block'; % specify condition to analyse.

if isequal(exp_condition, 'original')
    % specify original participant number (in parfile name).
    orig_ppnum = {'21', '04', '02', '25', '22', '16', '13', '01', '10', '19', '17'};
    
    % specify corresponding participant r-number.
    rnumber = [2268 2548 2590 2904 3111 3455 3517 3773 3932 4065 4496];
    
    % specify associated number of scans per participant.
    valid_runs = {1:6;1:6; 1:6; 1:6; [1 2 3 4 7]; 1:6; 1:7; 1:6; 1:6; 1:6; 1:6};
    
elseif isequal(exp_condition, 'colour')
    rnumber = [2268 2548 2590 2904 3111 3517 3773 4059 4065 4244 4831 4928];
    valid_runs = {1:8, 1:8, 1:7, 1:8, 1:7, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8, 1:8};
elseif isequal(exp_condition, 'naturalistic')
    rnumber = [2548 2904 3111 3517 4059 4065 4127 4244 4831 4833 4890 4928];
    valid_runs = {1:8, 1:8, 1:6, 1:8, 1:7, 1:8, 1:8, 1:7, 1:7, 1:8, 1:8, 1:6};
end

% specify new conn participant numbers.
new_ppnum = 1:length(rnumber);

% specify output file headings.
headers = {'condition_name', 'subject_number', 'session_number', 'onsets', 'durations'};

%% ------------------------ PROCESS DESIGN FILES ----------------------- %%

for i = 1:length(new_ppnum) % for each participant in turn:
    
    % specify the conn participant-specific output directory.
    if ppnum < 10
        if isequal(exp_condition, 'naturalistic')
            outputdir = sprintf('/scratch/groups/Projects/P1361/fMRI/conn_output/%s/participant_0%d/design/',condition, ppnum);
        else
            outputdir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/R0%d/conn/', exp_condition, rnumber(i));
        end
    else
        if isequal(exp_condition, 'naturalistic')
            outputdir = sprintf('/scratch/groups/Projects/P1361/fMRI/conn_output/%s/participant_%d/design/', condition, ppnum);
        else
            outputdir = sprintf('/scratch/groups/Projects/P1323/%s/fMRI/R%d/conn/', exp_condition, rnumber(i));
        end
    end
    [~,~] = mkdir(outputdir);
    
    % identify the run-specific parfiles.
    if ~isequal(exp_condition, 'naturalistic')
        parfiles = dir(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/R%d/mrvista/block/Stimuli/Parfiles/*.par',...
            exp_condition, rnumber(i)));
    elseif isequal(exp_condition, 'naturalistic')
        if isequal(condition, '3feature')
            if rnumber(i) == 3517
                parfiles = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/session1-allscans/Stimuli/Parfiles/*.par', rnumber(i)));
            else
                parfiles = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/Stimuli/Parfiles/*.par', rnumber(i)));
            end
        elseif isequal(condition, '3x3feature')
            if rnumber(i) == 3517
                parfiles = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/session1-allscans/multiple/Stimuli/Parfiles/*.par', rnumber(i)));
            else
                parfiles = dir(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/multiple/Stimuli/Parfiles/*.par', rnumber(i)));
            end
        end
    end
    
    parfilenames = {parfiles.name};
    
    for x = 1:length(valid_runs{i}) % for each run in turn:
        run = valid_runs{i}(x);
        
        % load in the run-specific parfile.
        if isequal(condition, 'colour-passive')
            [timing, event, labels] = textread(sprintf('/scratch/groups/Projects/P1323/colour/fMRI/R%d/mrvista/%s/Stimuli/Parfiles/%s',...
                rnumber(i), condition, parfilenames{x}),'%f %d %s');
            parfile = [timing, event];
        elseif isequal(condition, 'naturalistic')
            if isequal(condition, '3feature')
                if rnumber(i) == 3517
                    [timing, event, labels] = textread(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/session1-allscans/Stimuli/Parfiles/%s',...
                        rnumber(i), parfilenames{x}),'%f %d %s');
                else
                    [timing, event, labels] = textread(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/Stimuli/Parfiles/%s',...
                        rnumber(i), parfilenames{x}),'%f %d %s');
                end
            elseif isequal(condition, '3x3feature')
                if rnumber(i) == 3517
                    [timing, event, labels] = textread(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/session1-allscans/multiple/Stimuli/Parfiles/%s',...
                        rnumber(i), parfilenames{x}),'%f %d %s');
                else
                    [timing, event, labels] = textread(sprintf('/scratch/groups/Projects/P1361/fMRI/R%d/mrvista/multiple/Stimuli/Parfiles/%s',...
                        rnumber(i), parfilenames{x}),'%f %d %s');
                end
            end
            parfile = [timing, event];
        else
            parfile = load(sprintf('/scratch/groups/Projects/P1323/%s/fMRI/R%d/mrvista/%s/Stimuli/Parfiles/%s',...
                exp_condition, rnumber(i), condition, parfilenames{x}));
        end
        
        if isequal(condition, 'block') || isequal(condition, 'feature')
            
            % initalise timing and condition label variables.
            [interblock, orientation, contrast, shape, passive] = deal({i, run});
            [interblockonset, orientationonset, contrastonset, shapeonset, passiveonset] = deal([]);
            
            % store each parfile entry in a matrix containing only information
            % from that condition.
            for event = 1:length(parfile)
                if parfile(event,2) == 5
                    interblockonset = [interblockonset, parfile(event,1)];
                elseif parfile(event,2) == 4
                    passiveonset = [passiveonset, parfile(event,1)];
                elseif parfile(event,2) == 3
                    shapeonset = [shapeonset, parfile(event,1)];
                elseif parfile(event,2) == 2
                    contrastonset = [contrastonset, parfile(event,1)];
                elseif parfile(event,2) == 1
                    orientationonset = [orientationonset, parfile(event,1)];
                end
            end
            
            % add condition-labels to this segregated timing information.
            interblock = ['INTERBLOCK', interblock];
            orientation = ['ORIENTATION', orientation];
            contrast = ['CONTRAST', contrast];
            shape = ['SHAPE', shape];
            passive = ['PASSIVE', passive];
            
            % combine the condition labels, onsets and durations into singular
            % matrices for writing to file.
            interblock = [interblock, interblockonset, ones(1,length(interblockonset))*9];
            orientation = [orientation, orientationonset, ones(1,length(orientationonset))*15];
            contrast = [contrast, contrastonset, ones(1,length(contrastonset))*15];
            shape = [shape, shapeonset, ones(1,length(shapeonset))*15];
            passive = [passive, passiveonset, ones(1,length(passiveonset))*15];
            
            % specify column headers and combine information to a singular
            % matrix.
            data = [headers; interblock; orientation; contrast; shape; passive];
            
        % repeat the same process for the remaining conditions.
        elseif isequal(condition, 'colour-passive')
            
            [interblock, rgattention, rgpassive, byattention, bypassive, lumattention, lumpassive] = deal({ppnum, run});
            [interblockonset, rgattentiononset, rgpassiveonset, byattentiononset,...
                bypassiveonset, lumattentiononset, lumpassiveonset] = deal([]);
            
            for event = 1:length(parfile)
                if parfile(event,2) == 7
                    interblockonset = [interblockonset, parfile(event,1)];
                elseif parfile(event,2) == 1
                    rgattentiononset = [rgattentiononset, parfile(event,1)];
                elseif parfile(event,2) == 2
                    rgpassiveonset = [rgpassiveonset, parfile(event,1)];
                elseif parfile(event,2) == 3
                    byattentiononset = [byattentiononset, parfile(event,1)];
                elseif parfile(event,2) == 4
                    bypassiveonset = [bypassiveonset, parfile(event,1)];
                elseif parfile(event,2) == 5
                    lumattentiononset = [lumattentiononset, parfile(event,1)];
                elseif parfile(event,2) == 6
                    lumpassiveonset = [lumpassiveonset, parfile(event,1)];
                end
            end
            
            rgattention = ['RG-ATTENTION', rgattention];
            rgpassive = ['RG-PASSIVE', rgpassive];
            byattention = ['BY-ATTENTION', byattention];
            bypassive = ['BY-PASSIVE', bypassive];
            lumattention = ['LUM-ATTENTION', lumattention];
            lumpassive = ['LUM-PASSIVE', lumpassive];
            
            rgattention = [rgattention, rgattentiononset, ones(1,length(rgattentiononset))*15];
            rgpassive = [rgpassive, rgpassiveonset, ones(1,length(rgpassiveonset))*15];
            byattention = [byattention, byattentiononset, ones(1,length(byattentiononset))*15];
            bypassive = [bypassive, bypassiveonset, ones(1,length(bypassiveonset))*15];
            lumattention = [lumattention, lumattentiononset, ones(1,length(lumattentiononset))*15];
            lumpassive = [lumpassive, lumpassiveonset, ones(1,length(lumpassiveonset))*15];
            
            rgpassive{4} = [rgpassive{4}, NaN, NaN]; rgpassive{5} = [rgpassive{5}, NaN, NaN];
            bypassive{4} = [bypassive{4}, NaN, NaN]; bypassive{5} = [bypassive{5}, NaN, NaN];
            lumpassive{4} = [lumpassive{4}, NaN, NaN]; lumpassive{5} = [lumpassive{5}, NaN, NaN];
            
            data = [headers; rgattention; rgpassive; byattention; bypassive; lumattention; lumpassive];
            
        elseif isequal(condition, 'colourxfeature')
            
            [interblock, rgo, rgc, rgs, rgp, byo, byc, bys, byp, lumo, lumc, lums, lump] = deal({ppnum, run});
            [interblockonset, rgo_onset, rgc_onset, rgs_onset, rgp_onset,...
                byo_onset, byc_onset, bys_onset, byp_onset, lumo_onset,...
                lumc_onset, lums_onset, lump_onset] = deal([]);
            
            for event = 1:length(parfile)
                if parfile(event,2) == 13
                    interblockonset = [interblockonset, parfile(event,1)];
                elseif parfile(event,2) == 1
                    rgo_onset = [rgo_onset, parfile(event,1)];
                elseif parfile(event,2) == 2
                    rgc_onset = [rgc_onset, parfile(event,1)];
                elseif parfile(event,2) == 3
                    rgs_onset = [rgs_onset, parfile(event,1)];
                elseif parfile(event,2) == 4
                    rgp_onset = [rgp_onset, parfile(event,1)];
                elseif parfile(event,2) == 5
                    byo_onset = [byo_onset, parfile(event,1)];
                elseif parfile(event,2) == 6
                    byc_onset = [byc_onset, parfile(event,1)];
                elseif parfile(event,2) == 7
                    bys_onset = [bys_onset, parfile(event,1)];
                elseif parfile(event,2) == 8
                    byp_onset = [byp_onset, parfile(event,1)];
                elseif parfile(event,2) == 9
                    lumo_onset = [lumo_onset, parfile(event,1)];
                elseif parfile(event,2) == 10
                    lumc_onset = [lumc_onset, parfile(event,1)];
                elseif parfile(event,2) == 11
                    lums_onset = [lums_onset, parfile(event,1)];
                elseif parfile(event,2) == 12
                    lump_onset = [lump_onset, parfile(event,1)];
                end
            end
            
            rgo = ['RG-ORIENTATION', rgo];
            rgc = ['RG-CONTRAST', rgc];
            rgs = ['RG-SHAPE', rgs];
            rgp = ['RG-PASSIVE', rgp];
            byo = ['BY-ORIENTATION', byo];
            byc = ['BY-CONTRAST', byc];
            bys = ['BY-SHAPE', bys];
            byp = ['BY-PASSIVE', byp];
            lumo = ['LUM-ORIENTATION', lumo];
            lumc = ['LUM-CONTRAST', lumc];
            lums = ['LUM-SHAPE', lums];
            lump = ['LUM-PASSIVE', lump];
            
            rgo = [rgo, rgo_onset, ones(1,length(rgo_onset))*15];
            rgc = [rgc, rgc_onset, ones(1,length(rgc_onset))*15];
            rgs = [rgs, rgs_onset, ones(1,length(rgs_onset))*15];
            rgp = [rgp, rgp_onset, ones(1,length(rgp_onset))*15];
            byo = [byo, byo_onset, ones(1,length(byo_onset))*15];
            byc = [byc, byc_onset, ones(1,length(byc_onset))*15];
            bys = [bys, bys_onset, ones(1,length(bys_onset))*15];
            byp = [byp, byp_onset, ones(1,length(byp_onset))*15];
            lumo = [lumo, lumo_onset, ones(1,length(lumo_onset))*15];
            lumc = [lumc, lumc_onset, ones(1,length(lumc_onset))*15];
            lums = [lums, lums_onset, ones(1,length(lums_onset))*15];
            lump = [lump, lump_onset, ones(1,length(lump_onset))*15];
            
            data = [headers; rgo; rgc; rgs; rgp; byo; byc; bys; byp; lumo; lumc; lums; lump];
            
        elseif isequal(condition, '3feature')
            
            [orientation, colour, shape, passive, face] = deal({ppnum, run});
            [orientationonset, colouronset, shapeonset, passiveonset, faceonset] = deal([]);
            
            for event = 1:length(parfile)
                if parfile(event,2) == 2
                    orientationonset = [orientationonset, parfile(event,1)];
                elseif parfile(event,2) == 3
                    colouronset = [colouronset, parfile(event,1)];
                elseif parfile(event,2) == 4
                    shapeonset = [shapeonset, parfile(event,1)];
                elseif parfile(event,2) == 5
                    passiveonset = [passiveonset, parfile(event,1)];
                elseif parfile(event,2) == 1
                    faceonset= [faceonset, parfile(event,1)];
                end
            end
            
            orientation = ['ORIENTATION', orientation];
            colour = ['COLOUR', colour];
            shape = ['SHAPE', shape];
            passive = ['PASSIVE', passive];
            face = ['FACE', face];
            
            orientation = [orientation, orientationonset, ones(1,length(orientationonset))*30];
            colour = [colour, colouronset, ones(1,length(colouronset))*30];
            shape = [shape, shapeonset, ones(1,length(shapeonset))*30];
            passive = [passive, passiveonset, ones(1,length(passiveonset))*30];
            face = [face, faceonset, ones(1,length(faceonset))*30];
            
            % ensure all events are of the same length.
            if length(orientation{4}) < 3
                orientation{4} = [orientation{4},NaN]; orientation{5} = [orientation{4}, NaN];
            elseif length(colour{4}) < 3
                colour{4} = [colour{4},NaN]; colour{5} = [colour{5}, NaN];
            elseif length(shape{4}) < 3
                shape{4} = [shape{4}, NaN]; shape{5} = [shape{5}, NaN];
            end
            
            passive{4} = [passive{4}, NaN, NaN]; passive{5} = [passive{5}, NaN, NaN];
            face{4} = [face{4}, NaN, NaN]; face{5} = [face{5}, NaN, NaN];
            
            data = [headers; orientation; colour; shape; passive; face];
            
        elseif isequal(condition, '3x3feature')
            
            [horizontal, vertical, diagonal, red, green, blue, circular,...
                square, triangular,passive, face] = deal({ppnum, run});
            [horizontalonset, verticalonset, diagonalonset, redonset, greenonset, blueonset,...
                circularonset, squareonset, triangularonset, passiveonset, faceonset] = deal([]);
            
            for event = 1:length(parfile)
                if parfile(event,2) == 3
                    horizontalonset = [horizontalonset, parfile(event,1)];
                elseif parfile(event,2) == 2
                    verticalonset = [verticalonset, parfile(event,1)];
                elseif parfile(event,2) == 4
                    diagonalonset = [diagonalonset, parfile(event,1)];
                elseif parfile(event,2) == 5
                    redonset = [redonset, parfile(event,1)];
                elseif parfile(event,2) == 6
                    greenonset = [greenonset, parfile(event,1)];
                elseif parfile(event,2) == 7
                    blueonset = [blueonset, parfile(event,1)];
                elseif parfile(event,2) == 8
                    circularonset = [circularonset, parfile(event,1)];
                elseif parfile(event,2) == 9
                    squareonset = [squareonset, parfile(event,1)];
                elseif parfile(event,2) == 10
                    triangularonset = [triangularonset, parfile(event,1)];
                elseif parfile(event,2) == 1
                    passiveonset = [passiveonset, parfile(event,1)];
                elseif parfile(event,2) == 11
                    faceonset = [faceonset, parfile(event,1)];
                end
            end
            
            horizontal = ['HORIZONTAL', horizontal];
            vertical = ['VERTICAL', vertical];
            diagonal = ['DIAGONAL', diagonal];
            red = ['RED', red];
            green = ['GREEN', green];
            blue = ['BLUE', blue];
            circular = ['CIRCULAR', circular];
            square = ['SQUARE', square];
            triangular = ['TRIANGULAR', triangular];
            passive = ['PASSIVE', passive];
            face = ['FACE', face];
            
            horizontal = [horizontal, horizontalonset, ones(1,length(horizontalonset))*30];
            vertical = [vertical, verticalonset, ones(1,length(verticalonset))*30];
            diagonal = [diagonal, diagonalonset, ones(1,length(diagonalonset))*30];
            red = [red, redonset, ones(1,length(redonset))*30];
            green = [green, greenonset, ones(1,length(greenonset))*30];
            blue = [blue, blueonset, ones(1,length(blueonset))*30];
            circular = [circular, circularonset, ones(1,length(circularonset))*30];
            square = [square, squareonset, ones(1,length(squareonset))*30];
            triangular = [triangular, triangularonset, ones(1,length(triangularonset))*30];
            passive = [passive, passiveonset, ones(1,length(passiveonset))*30];
            face = [face, faceonset, ones(1,length(faceonset))*30];
            
            % ensure all events are of the same length.
            if length(horizontal) < 5
                horizontal{4} = NaN; horizontal{5} = NaN;
            elseif length(vertical) < 5
                vertical{4} = NaN; vertical{5} = NaN;
            elseif length(diagonal) < 5
                diagonal{4} = NaN; diagonal{5} = NaN;
            elseif length(red) < 5
                red{4} = NaN; red{5} = NaN;
            elseif length(green) < 5
                green{4} = NaN; green{5} = NaN;
            elseif length(blue) < 5
                blue{4} = NaN; blue{5} = NaN;
            elseif length(circular) < 5
                circular{4} = NaN; circular{5} = NaN;
            elseif length(square) < 5
                square{4} = NaN; square{5} = NaN;
            elseif length(triangular) < 5
                triangular{4} = NaN; triangular{5} = NaN;
            end
            
            data = [headers; horizontal; vertical; diagonal; red; green; blue;...
                circular; square; triangular; passive; face];
        end
        % write the reformatted timing information to a run-specific .csv
        % file.
        cell2csv(strcat(outputdir, sprintf('conn_pp%d_run%d.csv', new_ppnum(i), run)), data);
    end
end

%% --------------------------------------------------------------------- %%