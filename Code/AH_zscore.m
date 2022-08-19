function [neuron_event_zscored] = AH_zscore(allNeurons, event_idx, pre_event, post_event, z_score_window, zscore, ind_trial, trial_avg)
% Aniek van Hoek
%
% Goal: z-scoring of 1-photon calcium imaging according to defined strategy
% (zscore, ind_trial, trial_avg)
%
% Input
%       -   allNeurons:     C-matrix with raw neural data (# neurons x #
%                           time frames
%       -   event_idx:      Index start of event to define trials 
%                           (options: cue_idx, water_idx, ethanol_idx)
%       -   pre_event:      # timesframes to add to trial before event
%                           start
%       -   post_event:     # timeframes to add to trial after event start
%       -   z_score_window: # timeframes to use for z-scoring
%       -   zscore:         Whether to perform zscoring, binary (true or false)
%       -   ind_trial:      Whether to zscore per trial or with mean/std of all
%                           baselines combined, binary (true or false)
%       -   trial_avg:      Whether to zscore with mean/std of trial
%                           averaged neuronal data, binary (true or false)
%
%
% Output
%       -   neuron_event_zscored:   Neural data matrix 
%                                   (# neurons x # trials x # timepoints)

if zscore
    if ind_trial == false
        
        %% Z scoring by average of all baselines
        for i = 1:size(allNeurons,1) %iterate through each neuron
            neuronEventData_zscore = []; 
            clear baseline_mean; clear baseline_std;

            for j= 1:size(event_idx,2) %iterate through each trial

                    neuronEventData_zscore(j,:) = allNeurons(i,event_idx(j) - pre_event : event_idx(j) + post_event); 

                    baseline_mean(j) = mean(neuronEventData_zscore(j, 1:z_score_window));   %find mean of baseline window x sec before cue
                    baseline_std(j) = std(neuronEventData_zscore(j, 1:z_score_window));  

                    trial_starts(j) =  event_idx(j) - pre_event;
                    trial_ends(j) = event_idx(j) + post_event;

            end

            full_mean = mean(baseline_mean);
            full_std = mean(baseline_std);

            %neuron_trace_zscored(i,:) = (allNeurons(i,:) - full_mean) ./ full_std;
            neuron_event_zscored(i,:,:) = (neuronEventData_zscore - full_mean) ./ full_std;

            clear full_mean; clear full_std;
        end
    else
        %% Z scoring per trial baseline
        for i = 1:size(allNeurons,1) %iterate through each neuron
            neuronEventData_zscore = []; 

            for j= 1:size(event_idx,2) %iterate through each trial


                    neuronEventData_zscore(j,:) = allNeurons(i,event_idx(j) - pre_event : ...
                        event_idx(j)+ post_event);

                    baseline_mean = mean(neuronEventData_zscore(j, 1:z_score_window));   %find mean of baseline window 3 sec before cue
                    baseline_std = std(neuronEventData_zscore(j,  1:z_score_window));  %find std of baseline window

                    neuronEventData_zscore(j,:) = (neuronEventData_zscore(j,:)-baseline_mean)./baseline_std;

            end

            neuron_event_zscored(i,:,:) = neuronEventData_zscore; %save data for that neuron in a cell array

            clear neuronEventData_zscore;

        end
    end
else   
    %% Splitting in trials without zscoring 
    for i = 1:size(allNeurons,1) %iterate through each neuron
        neuronEventData_zscore = []; 

        for j= 1:size(event_idx,2) %iterate through each trial


                neuronEventData_zscore(j,:) = allNeurons(i,event_idx(j) - pre_event : ...
                    event_idx(j)+ post_event);
        end

        neuron_event_zscored(i,:,:) = neuronEventData_zscore; %save data for that neuron in a cell array

        clear neuronEventData_zscore;

    end    
    
    %% If zscoring using the trial averaged baseline (US only, for CS this is done in the msCA_singlemouse.m function since trial outcome information is necessary for averaging)
    if trial_avg

        for i = 1:size(neuron_event_zscored,1) % Run through all neurons
            for j = 1:size(neuron_event_zscored,3)
                   mean_session1(i,j) = mean(neuron_event_zscored(i,:,j),2); %(1:trial_match) use if trial matching
            end

            baseline_mean = mean(mean_session1(i, 1:z_score_window));   %find mean of baseline window 3 sec before cue
            baseline_std = std(mean_session1(i,  1:z_score_window));  %find std of baseline window

            neuron_event_zscored(i,:) = (mean_session1(i,:)-baseline_mean)./baseline_std;
        end
        
    end
end
    

