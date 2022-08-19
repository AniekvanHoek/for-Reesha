function drinking_data = retrieving_drinking_data(date,mousecage, mousenum,cue)
% Aniek van Hoek

% Goal: script retrieves information (just like loading_CI_data) to make
%       drinking matrix
%       
%       Drinking matrix is used for decoding analysis
%% Loading data
name_neuron = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/neuron/%s_M%d-%d_neuron.mat', date, mousecage, mousenum);
load(name_neuron) %C matrix from CNMFE analysis
allNeurons = neuron.C;

name_gpio = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Preprocessing/gpio/%s_%d.%d.csv', date, mousecage, mousenum);
raw_lick_data = readtable(name_gpio);
raw_lick_data.Properties.VariableNames = {'Time' 'Channel' 'Voltage'}; %add new column names to the raw_lick_data table


%% Setup the Import Options for MedPC files
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
name_medPC = sprintf("Reesha/MATLAB/cohort_9_calcium_imaging/Preprocessing/medpc/names/%s_%d.%d", date, mousecage, mousenum);
MedPC_Table = readtable(name_medPC, opts);

% Clear temporary variables
clear opts

%% set parameters
bout_window = 3; %set the window of no licks necessary to identify a new bout in seconds
window_before_cue = 3; %time to take before bout start in seconds
window_after_cue = 3; %time to take after first lick in bout in seconds
window_before_bout = 6; %time to take before bout start in seconds
window_after_bout= 10; %time to take after first lick in bout in seconds

LED_framerate = 10.06; %recorded at 20hz but temporal downsized by 2
z_score_window = 30; %number of frames in the extracted neural actiity window to use for zscore normalization
poststimulusWindow = 30; %set number to frames postlick to use to determine if the cell was responsive
responsive = 1.6;

shortest_frames = 1130; % Set to number of frames - 30 in trial with shortest ITI
number_trials = 30; % Set to number of cue lights/trials
num_type_trials = 5; % Set to number of different type trials

%% Prepare
trial_A = [2,10,17,19,26];
trial_B = [5,8,13,23,28];
trial_C = [1,3,4,6,7,9,11,12,14,15,16,18,20,21,22,24,25,27,29,30];
%trial_info = {trial_A, trial_B, trial_C};

%% Prepare final matrix
%drinking_data = NaN(5,4000);
drinking_data = NaN(4,length(allNeurons));
all_idx = 1:length(drinking_data);
%% Calculate Time Vector for Calcium Recording in Inscopix

Inscopix_GPIO_LEDon = raw_lick_data(strcmp(raw_lick_data.Channel,'EX-LED'),:);%INSCOPIX LEDon time array

start_time = Inscopix_GPIO_LEDon{2, 1}; %calculate delay before camera onset for real start time of recording
end_time = Inscopix_GPIO_LEDon{end-1,1};

LED_recording_time = end_time - start_time; %calculate total time of LED recording
LED_total_frames = length(allNeurons); %calculate the total number of frames in the LED recording

time = linspace(start_time, end_time, LED_total_frames);
drinking_data(1,:) = time;

%% Find Timestamps for Cue light in medPC
indices = table2array(MedPC_Table(:,1));  %extract indices from table
values = table2array(MedPC_Table(:,2));   %extract values from table

if cue 
    %for cue light
    cue_start = (find(contains(indices,'S'))+1);  %find the start of indices for cue saved at S in MedPC
    cue_start = cue_start(cue_start>25); % Keep only the indices that contain values 
    cue_end = (find(contains(indices, 'V'))-1);   %find the last indices for water
    cue_start = values(cue_start:cue_end,1);  %save the corresponding values of water lick start indices

    % Tranlate to Inscopix time
    cue_start = cue_start + start_time;


    %% Find indices in time for cue light in recording

    for i = 1:size(cue_start,1) %for each cue light
        idx = interp1(time,1:length(time),cue_start(i), 'nearest'); %find index in time for that bout
        cue_idx(i) = idx; %save value of index (if not found will be NaN)
    end

    drinking_data(2,:) = ismember(all_idx, cue_idx);
else
    drinking_data(2,:) = NaN;
end
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

% Convert to inscopix times
%waterLick_start_corr = waterLick_start + start_time;
waterLick_end_corr = waterLick_end + start_time;
%ethanolLick_start_corr = ethanolLick_start + start_time;
ethanolLick_end_corr = ethanolLick_end+ start_time;


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


%for ethanol
ethanol_interlickinterval = diff(ethanolLick_start); %calculate the inter-Lick interval
ethanolBouts_indices = (find(ethanol_interlickinterval>=bout_window))+1;
ethanolBouts_start = ethanolLick_start(ethanolBouts_indices,1);
firstEthanolLick_start = ethanolLick_start(1,1);
allEthanolBouts_start = vertcat(firstEthanolLick_start, ethanolBouts_start);    %timestamp for first lick in all water bouts
allEthanolBouts_start_corr = allEthanolBouts_start(:,1) + start_time;    %timestamp for first lick in all water bouts corrected for time delay in start of recording


%% Find end of bout times
% Water
for i = 1:length(allWaterBouts_start_corr)
    for j = 1:length(waterLick_end_corr)
        if i ~= length(allWaterBouts_start_corr)
            if waterLick_end_corr(j) < (allWaterBouts_start_corr(i+1))
                water_bout_end(i) = waterLick_end_corr(j);
            else
                break
            end
        else
            water_bout_end(i) = waterLick_end_corr(end);
        end
    end
end

% Ethanol
for i = 1:length(allEthanolBouts_start_corr)
    for j = 1:length(ethanolLick_end_corr)
        if i ~= length(allEthanolBouts_start_corr)
            if ethanolLick_end_corr(j) < (allEthanolBouts_start_corr(i+1))
                etoh_bout_end(i) = ethanolLick_end_corr(j);
            else
                break
            end
        else
            etoh_bout_end(i) = ethanolLick_end_corr(end);
        end
    end
end

%% find indices in time for start of water/ethanol bouts in recording
for i = 1: length(allWaterBouts_start_corr) %for each water bout
    idx = interp1(time,1:length(time),allWaterBouts_start_corr(i), 'nearest'); %find index in time for that bout
    water_start_idx(i) = idx; %save value of index (if not found will be NaN)
end

for i = 1: length(water_bout_end) %for each water bout
    idx = interp1(time,1:length(time),water_bout_end(i), 'nearest'); %find index in time for that bout
    water_end_idx(i) = idx; %save value of index (if not found will be NaN)
end

water_bout = [];
counter = 0;
for j = 1:length(water_start_idx)
    temp = water_start_idx(j):water_end_idx(j);
    water_bout((counter+1):(counter +length(temp)),1) = temp;
    counter = counter + length(temp);
    clear temp
end

for ii = 1: length(allEthanolBouts_start_corr) %for each ethanol bout
    idx = interp1(time,1:length(time),allEthanolBouts_start_corr(ii), 'nearest'); %find index in time for that bout
    ethanol_start_idx(ii) = idx; %save value of index (if not found will be NaN)
end

for ii = 1: length(etoh_bout_end) %for each ethanol bout
    idx = interp1(time,1:length(time),etoh_bout_end(ii), 'nearest'); %find index in time for that bout
    ethanol_end_idx(ii) = idx; %save value of index (if not found will be NaN)
end

ethanol_bout = [];
counter = 0;
for j = 1:length(ethanol_start_idx)
    temp = ethanol_start_idx(j):ethanol_end_idx(j);
    ethanol_bout((counter+1):(counter +length(temp)),1) = temp;
    counter = counter + length(temp);
    clear temp
end
drinking_data(3,:) = ismember(all_idx, water_bout);
drinking_data(4,:) = ismember(all_idx, ethanol_bout);

% Print information on licks to see if numbers are legit
fprintf("sum cue = %d\n", sum(drinking_data(2,:)));
fprintf("sum water = %d\n", sum(drinking_data(3,:)));
fprintf("sum etoh = %d\n", sum(drinking_data(4,:)));


end
