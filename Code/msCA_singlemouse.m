%% master script Calcium Analysis per individual mouse
% Aniek van Hoek
%
% Goal:     - This script uses indices of events to perform z-scoring for 
%           all mice per session.
%           - Processing of data goes per session for all animals at once

% Input:    Change parameters below         
%           Change the names of files and folders in name_zscored below
%                (within function)
%
%           if US and CS are both set to false only neuron_event_zscored is
%           saved (no plots)
%           
%
% Output:   Neural data matrix per trial (zscored)
%           
% Plots:    mean traces per trial over different x axis (time) to see if
%           zscoring/trial separateing was performed well and traces are valid.
%
clearvars; clc; close all;
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Change parameters 
%%% Retrieving of Data %%%
date = "02112022";
mouse_full = ["M4-2", "M4-3", "M5-1","M5-2","M5-3"]';
CS = false; 
US = false;
zscore = true; 
ind_trial = true; % Perform zscoring with baseline of individual trials
trial_avg = false; % Used for zscoring US within zscore function
trial_avg_CS = false; % if zscoring using trial average of all baselines for CS, probably set to false.

%% Analysis
for x = 1:5 % For every mouse  
mouse = mouse_full(x);
name_zscored = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/zscored/%s_%s_zscored_etoh', date, mouse);
name_zscored_CS = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/zscored/%s_%s_zscored_trialmatched', date, mouse);
       
%% Loading of neuron C matrix, inscopix time vector and medPC lick information
% Neural data: C & S  matrix from CNMFE analysis
name_neuron = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/neuron/%s_%s_neuron.mat', date, mouse);
load(name_neuron'); 
allNeurons = neuron.C;
allSpikes = neuron.S;

% Trial information
name_trial_info = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/trial_outcome/%s_%s_all_trials.mat', date, mouse);
load(name_trial_info); 

% Indices
name_cueidx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_cue_idx.mat', date, mouse);
load(name_cueidx); 

name_wateridx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_water_idx.mat', date, mouse);
load(name_wateridx); 

name_etohidx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_ethanol_idx.mat', date, mouse);
load(name_etohidx); 

clear name_neuron; clear name_trial_info; clear name_cueidx; 
clear name_wateridx; clear name_ethanolidx; 


%% set parameters
pre_event = 60; %frame window used pre event zscoring 
post_event = 100; % frame window post event z scoring --> can be -30 too
z_score_window = 30; %number of frames in the extracted neural activity window to use for zscore normalization
type_idx = ethanol_idx; %cue_idx, water_idx


%% Perform zscoring
% Check how zscoring function is currently set!
[neuron_event_zscored] = AH_zscore(allNeurons, type_idx, pre_event, post_event, z_score_window, zscore, ind_trial, trial_avg); % Use this one when output of zscoring

% Saving
save(name_zscored, 'neuron_event_zscored');

%% Analyses CS %%
if CS
%% Calculating where bouts happen within trials (CS defined) in the new trial structured neuron matrix 
    % Change number in all_trials{} to choose a trial type

    % Ethanol
    for j = 1:size(ethanol_idx,2)
        for i = 1:size(all_trials{5},2)
            idx = all_trials{5}(i);
            if ethanol_idx(j) > trial_starts(idx) && ethanol_idx(j) < trial_ends(idx)
                ethanol_idx_trial(j) = ethanol_idx(j) - trial_starts(idx);
            else
                continue
            end
        end
    end
    ethanol_idx_trial(ethanol_idx_trial == 0) = []; 

    % Water
    for j = 1:size(water_idx,2)
        for i = 1:size(all_trials{5},2)
            idx = all_trials{5}(i);
            if water_idx(j) > trial_starts(idx) && water_idx(j) < trial_ends(idx)
                water_idx_trial(j) = water_idx(j) - trial_starts(idx);
            else
                continue
            end
        end
    end
    water_idx_trial(water_idx_trial == 0) = []; 
    
 %% Trial averaging and then zscoring   
    if trial_avg_CS
        % Trial averaging all trials
        for i = 1:size(neuron_event_zscored,1) % Run through all neurons
            for j = 1:size(neuron_event_zscored,3)
                for k = [5] % For ethanol / sucrose drinking trials
                   type = all_trials{k};
                   mean_session1(i,j,k) = mean(neuron_event_zscored(i,type,j),2); %(1:trial_match) use if trial matching
                end
                for k = [1,2,3,4]  % For all other trials 
                    type = all_trials{k};
                    mean_session1(i,j,k) = mean(neuron_event_zscored(i,type,j),2);
                end
            end
        end

        % Z scoring using trial average of all baselines
        for i = 1:size(mean_session1,1) %iterate through each neuron
            for j= 1:size(mean_session1,3) %iterate through each trial

                    baseline_mean = mean(mean_session1(i, 1:z_score_window,j));   %find mean of baseline window 3 sec before cue
                    baseline_std = std(mean_session1(i,  1:z_score_window,j));  %find std of baseline window

                    neuron_event_zscored(i,:,j) = (mean_session1(i,:,j)-baseline_mean)./baseline_std;

            end
        end

        % Saving neuron_event_zscored of CS analysis
         save(name_zscored_CS, 'neuron_event_zscored');
    else
        
        %% Trial averaged neuron-averaged trace (first zscored, then trial averaged)
        trial_match = 5;    % Lowest number with trial outcome ethanol for all mice and all sessions

        % If trying a trial match lower than 5, change for loop k = [1...]
        for i = 1:size(neuron_event_zscored,1) % Run through all neurons
            for j = 1:num_frames
                for k = [5] 
                   type = all_trials{k};
                   mean_session1(i,j,k) = mean(neuron_event_zscored(i,type,j)); %(1:trial_match) use if trial matching
                end
                for k = [1,2,3,4]  
                    type = all_trials{k};
                    mean_session1(i,j,k) = mean(neuron_event_zscored(i,type,j));
                end
            end
        end
    end

    %% Plot  (used to see if traces seem valid and analysis was performed correctly)
    % Plots mean of all trials and all neurons
    x=60;
    figure
    suc = plot(mean(mean_session1(:,:,5)));
    err_1 = std(mean_session1(:,:,5))/sqrt(size(mean_session1(:,:,5),1));
    [hl, he] = errorbar_pn(1:size(mean_session1(:,:,5), 2), mean(mean_session1(:,:,5)), err_1, 111, 0.4);
    xlim([0 x])
    xline(z_score_window,'k','Linestyle','--','LineWidth',1);


    %% Plotting average neuron traces for each trial separately (etoh trial selected)
    mean_data = squeeze(mean(neuron_event_zscored,1));
    
    figure
    for i = 1:size(all_trials{5},2) % Change 5 to any number 1-5 to select trial outcome (5 is choice ethanol)
        plot(mean_data(i,:))
        hold on;
    end
    xline(z_score_window,'k','Linestyle','--','LineWidth',1);


    % Extraction of mean value per trial
    mean_trial = mean(neuron_event_zscored(:,:,51:664+50),3);
    mean_ITI = mean(neuron_event_zscored(:,:,664+51:1161),3);

    trial_mean = mean(mean_trial,1);
    ITI_mean = mean(mean_ITI,1);

    for i = 1:size(allNeurons,1)
        cue_start = 51; %31
        for j = 1:30
                cue_end = cue_start + 664;
                mean_neuron_trace(i,(j*2-1)) = mean(neuron_event_zscored(i,j, cue_start:cue_end)); 
                mean_neuron_trace(i, j*2) = mean(neuron_event_zscored(i,j,cue_end:end)); 
        end
    end

    % Plot full trace with means
    figure
    means = plot(mean(mean_neuron_trace,1));
    err_1 = std(mean_neuron_trace)/sqrt(size(mean_neuron_trace,1));
    [hl, he] = errorbar_pn(1:size(mean_neuron_trace, 2), mean(mean_neuron_trace), err_1, 111, 0.4);

end


%% Analyses US %%
if US
    mean_traces_etoh = neuron_event_zscored;
    mean_data = squeeze(mean(neuron_event_zscored,1));

     figure
     for i = 1:size(mean_data,2)
         plot(mean_data(i,:))
         hold on;
     end
     xline(z_score_window,'k','Linestyle','--','LineWidth',1);
     xline(pre_event,'r','Linestyle','--','LineWidth',1);


    %Average trace over all etoh licks
    trial_match = 10;
    for i = 1:size(neuron_event_zscored,1) % Run through all neurons
        for k = 1:size(neuron_event_zscored, 3) % Run through all timepoints
            mean_traces_etoh(i,k) = mean(neuron_event_zscored(i,1:trial_match,k)); 
        end
    end



    %% Plot (used to see if trace looks valid and analysis is performed correctly)
    % Plot first 90 seconds of trial
    time=90;
    figure
    suc = plot(mean(mean_traces_etoh(:,:)), 'r');
    err_1 = std(mean_traces_etoh(:,:))/sqrt(size(mean_traces_etoh(:,:),1));
    [hl, he] = errorbar_pn(1:size(mean_traces_etoh(:,:), 2), mean(mean_traces_etoh(:,:)), err_1, 222, 0.4);
    xline(z_score_window,'k','Linestyle','--','LineWidth',1);
    xline(pre_event,'k','Linestyle','--','LineWidth',1);
    xlim([0 time])

    % Plot full trial after bout
    figure
    suc = plot(mean(mean_traces_etoh(:,:)), 'r');
    err_1 = std(mean_traces_etoh(:,:))/sqrt(size(mean_traces_etoh(:,:),1));
    [hl, he] = errorbar_pn(1:size(mean_traces_etoh(:,:), 2), mean(mean_traces_etoh(:,:)), err_1, 222, 0.4);
    xline(z_score_window,'k','Linestyle','--','LineWidth',1);
    xline(pre_event,'k','Linestyle','--','LineWidth',1);
    xlim([1, pre_event+post_event])
end
end


