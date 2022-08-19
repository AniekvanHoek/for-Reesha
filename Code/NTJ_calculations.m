%% Calculations on neural trajectories
% Aniek van Hoek
%
% Goal:     Calculates length and distance of trial-averaged water/alcohol
%           bouts per session, for individual mice
%           Plots trial-averaged individual mice

% Input:    Change parameters below         
%           
%
% Output:  length_alcohol
%          lengt_water
%          distance_alcwater:   distance per timeframe between alcohol and
%                               water trial
%          num_pcs:             number of PC's needed to reach 90% variance
%                               explained

clearvars; clc
close all;

cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Parameters
T = 90;
timepoints = [1:90];
sucrose = true;
all = false; %Plots all sessions in 1 matrix to make between sessions comparisons

if sucrose
    dates = ["01082022", "02112022", "03012022"];
else
    dates = ["01232022", "02082022", "02282022"];
end

mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];

%% Loading Coregistered cells
reg_idx = cell(5,1);

if sucrose
    for i = 1:5
        filename = sprintf("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/celleg/Suc/%s_task_suc_cellreg.mat", mice(i));
        temp = load(filename);
        reg_idx{i} =  temp.cell_registered_struct.cell_to_index_map;
        reg_idx{i}(any(reg_idx{i}==0,2),:) = [];
    end
else
    for i = 1:5
        filename = sprintf("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/celleg/Etoh/%s_task_etoh_cellreg.mat", mice(i));
        temp = load(filename);
        reg_idx{i} =  temp.cell_registered_struct.cell_to_index_map;
        reg_idx{i}(any(reg_idx{i}==0,2),:) = [];
    end
end

%% Loading neurons of all mice

% Making two cell arrays containing necessary information of all sessions
% and all trials
for mouse=1:5
    for session = 1:3
        name_zscores = sprintf('nonzscored/%s_%s_etoh.mat', dates(session), mice(mouse));
        temp = load(fullfile(cd, name_zscores));
        neural_activity{mouse,session}   = temp.neuron_event;

        name_zscores_water = sprintf('nonzscored/%s_%s_water.mat', dates(session), mice(mouse));
        temp = load(fullfile(cd, name_zscores_water));
        neural_activity{mouse,session+3}   = temp.neuron_event;
    end
end

  

%% Getting neural data per trial
if all                                      % All sessions
    traces_per_trial = cell(size(mice,2),5); 
    for m = 1:size(neural_activity,1)
        for k =1:6
            registration = reg_idx{m};
            neuron_event_zscored = neural_activity{m,k};

            for i = 1:size(registration,1) % Run through all registrated neurons
                for j = 1:T  % Get timepoints used in decoding
                    if k >= 4
                        traces_per_trial{m,k}(i,j,:) = neuron_event_zscored(registration(i,k-3),:,timepoints(j));
                    else
                        traces_per_trial{m,k}(i,j,:) = neuron_event_zscored(registration(i,k),:,timepoints(j));
                    end

                end
            end
        end
    end 
else
    traces_per_trial = cell(size(mice,2),5);
    for m = 1:size(neural_activity,1)
        for k =1:6
            registration = reg_idx{m};
            neuron_event_zscored = neural_activity{m,k};

            for i = 1:size(neuron_event_zscored,1) % Run through all registrated neurons
                for j = 1:T  % Get timepoints used in decoding
                    traces_per_trial{m,k}(i,j,:) = neuron_event_zscored(i,:,timepoints(j));
                end
            end
        end
    end 
end

    

%% PCA all sessions
if all
for mouse = 1:5
    data_matrix = [];
    for session = 1:3    
        % PCA matrix
        data_matrix = horzcat(data_matrix, [mean(traces_per_trial{mouse,session+3}(:,31:90,:),3),mean(traces_per_trial{mouse,session}(:,31:90,:),3)]);  %31:90 because getting rid of z-scored baseline
    end   

    % Mean center and scale
    for i = 1:size(data_matrix,1)
        mean_neuron = mean(data_matrix(i, :));
        std_neuron = std(data_matrix(i,:));

        data_matrix(i,:) = (data_matrix(i,:) - mean_neuron)/std_neuron;
    end

    % Perform PCA
    [coef, ~, ~, ~, explained] = pca(data_matrix', 'Economy', false);
    num_pcs = find(cumsum(explained)>=90, 1); 

    % Project data onto PC's
    data_l(:, :) = coef(:, 1:num_pcs)'*data_matrix(:,:);
    data_water_ses1 = data_l(:,1:60); % Change numbers if changing timepoints
    data_alc_ses1  = data_l(:,61:120);

    data_water_ses2 = data_l(:,121:180);
    data_alc_ses2  = data_l(:,181:240);

    data_water_ses3 = data_l(:,241:300);
    data_alc_ses3  = data_l(:,301:360);


    distance_alcwater(mouse,1) = sum(sqrt(sum(((data_alc_ses1 - data_water_ses1).^2),1))); 
    distance_alcwater(mouse,2) = sum(sqrt(sum(((data_alc_ses2 - data_water_ses2).^2),1)));
    distance_alcwater(mouse,3) = sum(sqrt(sum(((data_alc_ses3 - data_water_ses3).^2),1)));
    
    clear coef explained data_l data_water_ses1 data_alc_ses1 data_water_ses2 data_alc_ses2 data_water_ses3 data_alc_ses3
end

else

%% Per session PCA
for session = 1:3    
for mouse = 1:5
    
        % PCA matrix
        data_matrix = [];
        data_matrix = vertcat(data_matrix, [mean(traces_per_trial{mouse,session+3}(:,31:90,:),3),mean(traces_per_trial{mouse,session}(:,31:90,:),3)]);
        

        % Mean center and scale
        for i = 1:size(data_matrix,1)
            mean_neuron = mean(data_matrix(i, :));
            std_neuron = std(data_matrix(i,:));

            data_matrix(i,:) = (data_matrix(i,:) - mean_neuron)/std_neuron;
        end
        
        % Perform PCA
        [coef, ~, ~, ~, explained] = pca(data_matrix', 'Economy', false);
        num_pcs(mouse,session) = find(cumsum(explained)>=90, 1); 
        
        % Project data onto PC's
        data_l(:, :) = coef(:, 1:num_pcs)'*data_matrix(:,:);
        data_water = data_l(:,1:60);
        data_alc = data_l(:,61:120);
     
        length_alcohol(mouse,session) = sum(sqrt(sum((diff(data_alc,1,2).^2),1)));
        lengt_water(mouse,session) = sum(sqrt(sum((diff(data_water,1,2).^2),1)));
        distance_alcwater (mouse,session) = sum(sqrt(sum(((data_alc - data_water).^2),1))); 
        
        %% Plot trajectories
        % Get data
        counter = 0;
        data = {};
        T=60;
        for i = 1:2
            data = [data, data_l(1:3, counter+1:counter+T)];
            counter = counter + T;
        end
        
        % parameters
        smooth_window = [10 0];
        event_markers = [1, 30, 60];
        event_marker_size = 60;
        line_opacity = 1;
        line_width = 3;
        errorbar_opacity = 1;
        point_size = 5;
        point_freq = 5;

        if sucrose
            line_plot_colors = {[0, 0, 255]/255, [5,153,72]/255};
        else
            line_plot_colors = {[0, 0, 255]/255, [255,0,0]/255};
        end
        scatter_plot_colors = {[0, 0, 0]/255, [255,255,255]/255};
        

        % Plot
        figure
        plot2traj_3Dgaussian_AH(data, smooth_window, {"Water", "Alcohol"}, event_markers, event_marker_size, ... 
        {"Start", "Bout", "End"}, line_plot_colors, line_opacity, line_width, scatter_plot_colors, point_size, point_freq)
        xlabel(sprintf("PC1: %0.2f%%", round(explained(1),2)))
        ylabel(sprintf("PC2: %0.2f%%", round(explained(2),2)))
        zlabel(sprintf("PC3: %0.2f%%", round(explained(3),2)))
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
        set(gcf,'Position',[20 20 1000 1000])
             
             
        clear coef explained data_l data_water data_alc
    end
end
end

