%% Loading of Calcium Imaging Data
% Aniek van Hoek
%
% Goal:     This script combines Inscopix, MED-PC and Neural data to get
%           the indices of the cue, water bouts and ethanol bouts
%
% Input:    Change the date in the script below
%           Loading of data goes per session for all animals at once
% Output:   -   Saved indices per type_idx
%           -   Cell array with trial outcome
%

clc; clearvars;
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));

%%% Retrieving of Data %%%
date = "01082022";
mouse_full = ["4.2", "4.3", "5.1","5.2","5.3"];
mouse_full_save = ["M4-2", "M4-3", "M5-1","M5-2","M5-3"];

for x = 1:5 % For all mice
mouse_load = mouse_full(x);
mouse_save = mouse_full_save(x);

%% Getting Neural data
name_neuron = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/neuron/%s_%s_neuron.mat', date, mouse_save);
load(name_neuron'); 
allNeurons = neuron.C;

%% Getting Inscopix CSV file
name_gpio = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Preprocessing/gpio/%s_%s.csv', date, mouse_load);
raw_lick_data = readtable(name_gpio); %table of licks from INSCOPIX
raw_lick_data.Properties.VariableNames = {'Time' 'Channel' 'Voltage'}; %add new column names to the raw_lick_data table

Inscopix_GPIO_LEDon = raw_lick_data(strcmp(raw_lick_data.Channel,'EX-LED'),:);%INSCOPIX LEDon time array

start_time = Inscopix_GPIO_LEDon{2, 1}; %calculate delay before camera onset for real start time of recording
end_time = Inscopix_GPIO_LEDon{end-1,1}; %start_time + 3935; %Inscopix_GPIO_LEDon{end-1,1}; %

LED_recording_time = end_time - start_time; %calculate total time of LED recording
LED_total_frames = length(allNeurons); %calculate the total number of frames in the LED recording

time = linspace(start_time, end_time, LED_total_frames);


%% Getting MEDPC file
% Setup the Import Options for MedPC files
opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["File", "CMEDPC", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"];
opts.SelectedVariableNames = ["File", "CMEDPC"];
opts.VariableTypes = ["string", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "TrimNonNumeric", true);
opts = setvaropts(opts, 2, "ThousandsSeparator", ",");
opts = setvaropts(opts, [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
name_medpc = sprintf("Reesha/MATLAB/cohort_9_calcium_imaging/Preprocessing/medpc/names/%s_%s", date, mouse_load);
MedPC_Table = readtable(name_medpc, opts);

% Clear temporary variables
clear opts


%% Getting indices for cue, ethanol and water licks
indices = table2array(MedPC_Table(:,1));  %extract indices from table
values = table2array(MedPC_Table(:,2));   %extract values from table
bout_window = 3;

%% for cue light
cue_start = (find(contains(indices,'S'))+1);  
cue_start = cue_start(cue_start>25); 
cue_end = (find(contains(indices, 'V'))-1);  
cue_start = values(cue_start:cue_end,1);  

% Tranlate to Inscopix time so that it corresponds to time vector created
% above
cue_start = cue_start + start_time;

% Find closest frame in time vector for cue light in recording

for i = 1:size(cue_start,1) %for each cue light
    idx = interp1(time,1:length(time),cue_start(i), 'nearest'); %find index in time for that bout
    cue_idx(i) = idx; %save value of index (if not found will be NaN)
end

% Saving cue_idx
name_cueidx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_cue_idx', date, mouse_save);
save(name_cueidx, 'cue_idx');

%% Find Timestamps for licks in medPC
%for water licks
q_start = (find(contains(indices,'Q'))+1);  %find the start of indices for water saved at Q in MedPC
q_end = (find(contains(indices, 'R'))-1);   %find the last indices for water
waterLick_start = values(q_start:q_end,1);  %save the corresponding values of water lick start indices

r_start = (find(contains(indices,'R'))+1);
r_end = (find(contains(indices, 'S'))-1);
r_end = r_end(r_end>25);                    % Make sure no information in first 25 rows is used
waterLick_end = values(r_start:r_end,1);

%for ethanol licks
p_start = (find(contains(indices, 'P'))+1);
p_end = (find(contains(indices, 'Q'))-1);
ethanolLick_start = values(p_start:p_end);

n_start = (find(contains(indices, 'N:'))+1);
n_start = n_start(n_start>25);
n_end = (find(contains(indices, 'P'))-1);
ethanolLick_end = values(n_start:n_end,1);


%% %% Find timestamps for the beginning of a water or ethanol lick bout for neural activity alignment
% beginning of lick bout is defined as any lick after 3 seconds of no licks
% timestamps are adjusted to recording start time delay

%for water
water_interlickinterval = diff(waterLick_start); %calculate the inter-Lick interval
waterBouts_indices = (find(water_interlickinterval>=bout_window))+1;
waterBouts_start = waterLick_start(waterBouts_indices,1);
firstWaterLick_start = waterLick_start(1,1);
allWaterBouts_start = vertcat(firstWaterLick_start, waterBouts_start);    %timestamp for first lick in all water bouts
allWaterBouts_start_corr = allWaterBouts_start(:,1)+ start_time;    %timestamp for first lick in all water bouts corrected for time delay in start of recording
%waterBouts_start_corr = allWaterBouts_start(allWaterBouts_start > (start_time + window_before),1);

%for ethanol
ethanol_interlickinterval = diff(ethanolLick_start); %calculate the inter-Lick interval
ethanolBouts_indices = (find(ethanol_interlickinterval>=bout_window))+1;
ethanolBouts_start = ethanolLick_start(ethanolBouts_indices,1);
firstEthanolLick_start = ethanolLick_start(1,1);
allEthanolBouts_start = vertcat(firstEthanolLick_start, ethanolBouts_start);    %timestamp for first lick in all water bouts
allEthanolBouts_start_corr = allEthanolBouts_start(:,1) + start_time;    %timestamp for first lick in all water bouts corrected for time delay in start of recording
%ethanolBouts_start_corr = allEthanolBouts_start(allEthanolBouts_start > (start_time + window_before),1);

%% find indices in time for start of water/ethanol bouts in recording
water_idx_wide = NaN(3,length(allWaterBouts_start_corr));
ethanol_idx_wide = NaN(3,length(allEthanolBouts_start_corr));

for i = 1: length(allWaterBouts_start_corr) %for each water bout
    idx = interp1(time,1:length(time),allWaterBouts_start_corr(i), 'nearest'); %find index in time for that bout
    water_idx_wide(2,i) = idx; %save value of index (if not found will be NaN)
end
water_idx = water_idx_wide(2,:);

for ii = 1: length(allEthanolBouts_start_corr) %for each ethanol bout
    idx = interp1(time,1:length(time),allEthanolBouts_start_corr(ii), 'nearest'); %find index in time for that bout
    ethanol_idx_wide(2,ii) = idx; %save value of index (if not found will be NaN)
end
ethanol_idx = ethanol_idx_wide(2,:);

% Saving water and ethanol indices
name_wateridx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_water_idx', date, mouse_save);
save(name_wateridx, 'water_idx');

name_etohidx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_ethanol_idx', date, mouse_save);
save(name_idx, 'ethanol_idx');



%% Getting trial outcome
trial_A = [2,10,17,19,26];
trial_B = [5,8,13,23,28];
trial_C = [1,3,4,6,7,9,11,12,14,15,16,18,20,21,22,24,25,27,29,30];

% Find which trial each bout belongs too

% Water
for i = 1:size(water_idx_wide,2)
    for j = 1:size(cue_idx,2)
        if j < size(cue_idx,2)
            if (water_idx_wide(2,i) > cue_idx(j)) && (water_idx_wide(2,i)< cue_idx(j+1))
               water_idx_wide(1,i) = j; 
               water_idx_wide(3,i) = 200;
               break
            end
        else
            if (water_idx_wide(2,i) > cue_idx(30))
                water_idx_wide(1,i) = 30;
                water_idx_wide(3,i) = 200;
            end
        end
    end
end

% Ethanol
for i = 1:size(ethanol_idx_wide,2)
    for j = 1:size(cue_idx,2)
        if j < size(cue_idx,2)
            if (ethanol_idx_wide(2,i) > cue_idx(j)) && (ethanol_idx_wide(2,i)< cue_idx(j+1))
               ethanol_idx_wide(1,i) = j; 
               ethanol_idx_wide(3,i) = 300;
               break
            end
        else
            if (ethanol_idx_wide(2,i) > cue_idx(30))
                ethanol_idx_wide(1,i) = 30;
                ethanol_idx_wide(3,i) = 300;
            end
        end
    end
end

% No drink
none_idx = [];
for i = 1:20
    choice_trial = trial_C(i);
    if ~ismember(choice_trial, ethanol_idx_wide(1,:)) && ~ismember(choice_trial, water_idx_wide(1,:))
        none_idx = [none_idx, choice_trial];
    end
end

for i =1:size(none_idx)
    none_idx(3,:) = 100;
    none_idx(2,:) = NaN;
end

% Transpose
water_idx_wide = water_idx_wide.';
ethanol_idx_wide = ethanol_idx_wide.';
none_idx = none_idx';

% Final table
info = vertcat(water_idx_wide, ethanol_idx_wide, none_idx);
info = info(any((info(:,1)== trial_C),2),:);
info = sortrows(info,[1,2]);

%% Divide choice trials by first lick
included = [];
trial_C_water = [];
trial_C_etoh = [];
trial_C_none = [];

for i = 1:size(info,1)
    if ismember(info(i,1),included)
        continue
    elseif info(i,3) == 100
        trial_C_none = [trial_C_none, info(i,1)];
    elseif info(i,3) == 200
        trial_C_water = [trial_C_water, info(i,1)];
        included = [included, info(i,1)];
    elseif info(i,3) == 300
        trial_C_etoh = [trial_C_etoh, info(i,1)];    
        included = [included, info(i,1)];
    end
end

% Saving trial outcome
all_trials = {trial_A, trial_B, trial_C_none, trial_C_water, trial_C_etoh};
name_trialoutcome = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/trial_outcome/%s_%s_all_trials', date, mouse_save);
save(name_trialoutcome, 'all_trials');

end