%% Dimensionality reduction using different methods
% Aniek van Hoek
%
% Goal:         Dimensionality reduction using different techniques
%
% Input:        Choose which technique and other settings in section parameters 
%
% Ouput:        Plots depicting neural data of alcohol/water bouts in reduced
%               state space

clearvars;
clc;
close all;

cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Parameters
T = 90;
timepoints = [1:90];
sucrose = false;
%
if sucrose
    dates = ["01082022", "02112022", "03012022"];
else
    dates = ["01232022", "02082022", "02282022"];
end

mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];

%% Settings
sessions_used = 1;
mice_used = 1:5;   
bouts_used = 1:10; % Always 1:4, unless doing PCA earlylate: 1:2

trial_avg = false; % Using trial averaged data
timematch = false; % Time matching bouts for individual trials
US = true; % For PCA analysis
CS = false; % For PCA analysis
MCS = true; % Mean centre and scale data
MC = false; % Mean centre data
earlylate = false; % Taking the earliest and latest trials within session

%% Technique
PLSR = false; % Use 1:10 for bouts_used, 1 session
PCA = false; % Individual vs all mice, US per trial (timematched or earlylate with 2 bouts) vs trial averaged
PCDA = false; % Use timematched, 1:4 for bouts_used, 1 session
LDA = true; % Individual bouts, timematched or earlylate
UMAP = false; % Individual bouts, timematched or earlylate


%% Coregistered cells
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
if sessions_used == 1:3
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


    
%% Find bouts which are closests to each other, divided over the full trial

if timematch
for session =1:3
for mouse = 1:5
% Find closest in time
name_wateridx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_water_idx.mat', dates(session), mice(mouse));
load(name_wateridx); 

name_etohidx = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/indices/%s_%s_ethanol_idx.mat', dates(session), mice(mouse));
load(name_etohidx); 

% Find parameters for dividing bouts over 5 equal sessions
highest_idx = 39364; %max([ethanol_idx, water_idx]); %% FOR SOME REASON DIDNT WORK WHEN USING MAX. FIND OUT LATER WHY.
steps = floor(highest_idx/5);
lb = 1 : steps : highest_idx-steps;
ub = highest_idx/5 : steps : highest_idx;


% Find bouts closests in times in water and ethanol trial 
[minValue_etoh, closestIndex_etoh] = min(abs(water_idx - ethanol_idx.')); % Change order of trial arrays depending on which array is smallest
closestValue_etoh = ethanol_idx(closestIndex_etoh) ;

[minValue_water, closestIndex_water] = min(abs(ethanol_idx - water_idx.')); % Change order of trial arrays depending on which array is smallest
closestValue_water = water_idx(closestIndex_water);


per_cat_etoh = cell(5,1);
for i = 1:size(closestValue_etoh,2)
    for j = 1:5
        if closestValue_etoh(i) > lb(j) && closestValue_etoh(i) < ub(j)
            temp_array_etoh = per_cat_etoh{j};
            new_array_etoh = [temp_array_etoh, closestIndex_etoh(i)];
            per_cat_etoh{j} = new_array_etoh;
        end
    end
end

for j = 1:5
    try
        temp_idx_list = find(ismember(closestIndex_etoh, per_cat_etoh{j}));
        [~, temp_best_idx] = min(minValue_etoh(temp_idx_list));
        best_etoh_idx(j) = closestIndex_etoh(temp_idx_list(temp_best_idx));
    catch ME
        fprintf('No indices in category %d, mouse %d, session %d.\n Error message: %s\n', j, mouse, session, ME.message);
        best_water_idx(j) = NaN;
        continue
    end
end




per_cat_water = cell(5,1);
for i = 1:size(closestValue_water,2)
    for j = 1:5
        if closestValue_water(i) > lb(j) && closestValue_water(i) < ub(j)
            temp_array_water = per_cat_water{j};
            new_array_water = [temp_array_water, closestIndex_water(i)];
            per_cat_water{j} = new_array_water;
        end
    end
end

for j = 1:5
    try
        temp_idx_list_water = find(ismember(closestIndex_water, per_cat_water{j}));
        [~, temp_best_idx_water] = min(minValue_water(temp_idx_list_water));
        best_water_idx(j) = closestIndex_water(temp_idx_list_water(temp_best_idx_water));
    catch ME
        fprintf('No indices in category %d, mouse %d, session %d.\n Error message: %s\n', j, mouse, session, ME.message);
        best_water_idx(j) = NaN;
        continue
    end
end

try
total_indices = horzcat(best_water_idx.', best_etoh_idx.');
total_indices(any(isnan(total_indices), 2), :) = [];
catch ME
    fprintf("Error Message: %s\n", ME.message);
    best_water_idx = best_water_idx(~isnan(best_water_idx));
    best_etoh_idx = best_etoh_idx(~isnan(best_etoh_idx));
    total_indices = horzcat(best_water_idx.', best_etoh_idx.');
    total_indices(any(isnan(total_indices), 2), :) = [];
end

traces_per_trial{mouse, session} = traces_per_trial{mouse,session}(:,:,total_indices(:,2));
traces_per_trial{mouse, session+3} = traces_per_trial{mouse,session+3}(:,:,total_indices(:,1));

end
end
end


%% Data matrix of data used
trial_match = 4;

%Non trial averaged data
data_matrix_w_all = []; 
data_matrix_e_all = [];
data_matrix_ses = {};
data_matrix = [];

if trial_avg == false
if earlylate 
    for m = mice_used
        data_matrix_both = [];


        for session = sessions_used
                data_matrix_w_early = [];
                data_matrix_e_early = [];
                data_matrix_w_late = [];
                data_matrix_e_late = [];
                
                % Make PCA matrix
                for t = bouts_used
                    data_matrix_w_early = cat(2, data_matrix_w_early, traces_per_trial{m,session+3}(:,31:90,t));
                    data_matrix_w_late = cat(2, data_matrix_w_late, traces_per_trial{m,session+3}(:,31:90,end+1-t));
                    data_matrix_e_early= cat(2, data_matrix_e_early, traces_per_trial{m,session}(:,31:90,t));
                    data_matrix_e_late = cat(2, data_matrix_e_late, traces_per_trial{m,session}(:,31:90,end+1-t));
                end

        data_matrix_both = horzcat(data_matrix_both, data_matrix_w_early, data_matrix_w_late, data_matrix_e_early, data_matrix_e_late);
        end

        data_matrix = vertcat(data_matrix, data_matrix_both);
    end
else
    for m = mice_used
        data_matrix_both = [];


        for session = sessions_used
                data_matrix_w = [];
                data_matrix_e =[];

                % Make data matrix
                for t = bouts_used
                    data_matrix_w = cat(2, data_matrix_w, traces_per_trial{m,session+3}(:,31:90,t));                 
                    data_matrix_e= cat(2, data_matrix_e, traces_per_trial{m,session}(:,31:90,t));
                    
                end

        data_matrix_both = horzcat(data_matrix_both, data_matrix_w, data_matrix_e);
        end

        data_matrix = vertcat(data_matrix, data_matrix_both);
    end

end
else % Trial averaged data
    for i = mice_used
        data_matrix = vertcat(data_matrix, [mean(traces_per_trial{i,sessions_used+3}(:,31:90,:),3),mean(traces_per_trial{i,sessions_used}(:,31:90,:),3)]);
    end
end

 
%% Mean center and scale

if MCS
    for i = 1:size(data_matrix,1)
        mean_neuron = mean(data_matrix(i, :));
        std_neuron = std(data_matrix(i,:));

        data_matrix(i,:) = (data_matrix(i,:) - mean_neuron)/std_neuron;
    end
end

if MC
    for i = 1:size(data_matrix,1)
    mean_neuron = mean(data_matrix(i, :));
    data_matrix(i,:) = (data_matrix(i,:) - mean_neuron);
    end
end

%% BELOW Different techniques start 
%% PLSR (Partial least squares regression)
if PLSR
    if earlylate == false 
       classification = [ones(600,1);zeros(600,1)]; %10 bouts of water and 10 bouts of alcohol 
    
       [Xloadings, Yloadings, Xscores, Yscores, betaPLS10, PLSPctVar] = plsregress(data_matrix.', classification, 3);
        
        Y2 = data_matrix.'*Xloadings; % Project data onto new components
       
        % Plotting
        line_plot_colors = 'bbbbbbbbbbrrrrrrrrrr' % 10 trials water, 10 trials alcohol
        T=60;
        counter = 0;
        
        figure
        for i = 1:20
            waitforbuttonpress()
            scatter(Y2(counter+1:counter+T,1), Y2(counter+1:counter+T,2), [], line_plot_colors(i))
            hold on
            counter = counter + T;
        end
        xlabel(sprintf("PC1: %0.2f%%", round(PLSPctVar(1,1),2)))
        ylabel(sprintf("PC2: %0.2f%%", round(PLSPctVar(1,2),2)))      
        title(sprintf("PLSR session %d", sessions_used)); 
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
        
    end
end

%% PCDA (PCA and then LDA)
% Use time matched bouts
if PCDA
    [coef, scores, ~, ~, explained] = pca(data_matrix', 'Economy', false);
     num_pcs = find(cumsum(explained)>=90, 1); 
     new_data = coef(:, 1:num_pcs)'*data_matrix(:,:);
     classification = [ones(240,1);zeros(240,1)]; % 4 bouts water, 4 bouts alcohol
     
        
    Mdl = fitcdiscr(new_data.', classification);
    [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order! 
    lambda = diag(LAMBDA);
    [lambda, SortOrder] = sort(lambda, 'descend');
    W = W(:, SortOrder);
    Y = new_data.'*W;
    
    % Plotting 4 trials of water, 4 of alcohol
    line_plot_colors = [[185 235 255]/255; [145 224 255]/255; [0 171 240]/255; [0 71 100]/255; [255 204 204]/255; [255 102 102]/255;  [255 0 0]/255; [204 0 0]/255];
    T=60;
    counter = 0;
    figure
    for i = 1:8
        scatter3(Y(counter+1:counter+T,1), Y(counter+1:counter+T,2), Y(counter+1:counter+T,3), [], line_plot_colors(i,:))
        hold on
        counter = counter + T;
    end
end

%% LDA
if LDA
    if earlylate
        % Classify per bout
        classification = [repmat([1],240,1); repmat([2],240,1); repmat([3],240,1); repmat([4],240,1)]; % 4 early water, 4 late water, same for alcohol

        Mdl = fitcdiscr(data_matrix.', classification, 'gamma',1 );
        [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order! 
        lambda = diag(LAMBDA);
        [lambda, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        Y = data_matrix.'*W;

        % Plotting
        line_plot_colors = [[185 235 255]/255; [0 71 100]/255; [255 204 204]/255; [204 0 0]/255];
        T=240;
        counter = 0;
        figure
        for i = 1:4
            scatter3(Y(counter+1:counter+T,1), Y(counter+1:counter+T,2), Y(counter+1:counter+T,3), [], line_plot_colors(i,:))
            hold on
            counter = counter + T;
        end
    else
        classification = [ones(240,1); zeros(240,1)]; %4 bouts water, 4 bouts alcohol
        Mdl = fitcdiscr(data_matrix.', classification); 
        [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order! 
        lambda = diag(LAMBDA);
        [lambda, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        Y = data_matrix.'*W; 
       
       % Plotting
       line_plot_colors = [[185 235 255]/255; [145 224 255]/255; [0 171 240]/255; [0 71 100]/255; [255 204 204]/255; [255 102 102]/255;  [255 0 0]/255; [204 0 0]/255];
        T=60;
        counter = 0;
        figure

        for i = 1:20
            waitforbuttonpress();
            scatter3(Y(counter+1:counter+T,1), Y(counter+1:counter+T,2),Y(counter+1:counter+T,3), [], line_plot_colors(i))
            hold on
            counter = counter + T;
        end
        
    end
    
 
        
end
%% UMAP
if UMAP
    [reduction, umap, clusterIdentifiers, extras] = run_umap(data_matrix.', 'cluster_output', 'graphic', ...
        'cluster_detail', 'dbscan arguments',...
        'epsilon', .6, ...
        'minpts', 5);


    cluster_UMAP = kmeans(reduction,2);
    
    % Clustered figure
    figure
    scatter(reduction(cluster_UMAP==1,1), reduction(cluster_UMAP==1,2),'s')
    hold on
    scatter(reduction(cluster_UMAP==2,1), reduction(cluster_UMAP==2,2),'*')
    
    if earlylate
        T = 60;
        line_plot_colors = [repmat([185 235 255]/255, 4,1); repmat([0 71 100]/255, 4,1); repmat([255 204 204]/255, 4,1); repmat([204 0 0]/255, 4,1);];
        counter = 0;
        figure
        for i = 1:16
            scatter(reduction(counter+1:counter+T,1), reduction(counter+1:counter+T,2),[], line_plot_colors(i,:))
            hold on
            counter = counter + T;
        end
    else      
        % Figure colored by bout number
        T= 60;
        line_plot_colors = [[185 235 255]/255; [145 224 255]/255; [0 171 240]/255; [0 71 100]/255; [255 204 204]/255; [255 102 102]/255;  [255 0 0]/255; [204 0 0]/255];
        counter = 0;
        figure
        for i = 1:8
            waitforbuttonpress;
            scatter(reduction(counter+1:counter+T,1), reduction(counter+1:counter+T,2),[], line_plot_colors(i,:))
            hold on
            counter = counter + T;
        end
    end
end
%% PCA
if PCA
    [coef, ~, ~, ~, explained] = pca(data_matrix', 'Economy', false);
    num_pcs = find(cumsum(explained)>=90, 1); 
    
    if mice_used == 1:5

        %% Getting animal_ID per cell
        cell_ID = [];
        for i = 1:size(traces_per_trial,1)
            cell_num = size(traces_per_trial{i,1}, 1);
            cell_ID = vertcat(cell_ID, i.*ones(cell_num,1));
        end

        
        % PCA all mice, CS, trial_averaged
        if CS
        %% Getting LOO data CS
        water_forced = pca2LOO(coef, data_matrix(:,1:60), cell_ID, num_pcs);
        etoh_forced = pca2LOO(coef, data_matrix(:,61:120), cell_ID, num_pcs);
        choice_none = pca2LOO(coef, data_matrix(:,121:180), cell_ID, num_pcs);
        choice_water = pca2LOO(coef, data_matrix(:,181:240), cell_ID, num_pcs);
        choice_etoh = pca2LOO(coef, data_matrix(:,241:300), cell_ID, num_pcs);

        plot_data = {water_forced([1:3],:,:), etoh_forced([1:3],:,:), choice_none([1:3],:,:),...
            choice_water([1:3],:,:), choice_etoh([1:3],:,:)}; %choice_none([4,5],:,:)
        end

        
        % PCA all mice, US
        if US
        %% LOO data US
        
        % All sessions
        if session == 1:3
            water_US_ses1 = pca2LOO(coef, data_matrix(:,1:240), cell_ID, num_pcs);
            etoh_US_ses1 = pca2LOO(coef, data_matrix(:,241:480), cell_ID, num_pcs);
            water_US_ses2 = pca2LOO(coef, data_matrix(:,481:720), cell_ID, num_pcs);
            etoh_US_ses2 = pca2LOO(coef, data_matrix(:,721:960), cell_ID, num_pcs);
            water_US_ses3 = pca2LOO(coef, data_matrix(:,961:1200), cell_ID, num_pcs);
            etoh_US_ses3 = pca2LOO(coef, data_matrix(:,1201:1440), cell_ID, num_pcs);

            data_l = horzcat(water_US_ses1, etoh_US_ses1, water_US_ses2, etoh_US_ses2, water_US_ses3, etoh_US_ses3); %, water_US_ses2, etoh_US_ses2, water_US_ses3, etoh_US_ses3);
        
        % One session    
        else
            water_US_ses1 = pca2LOO(coef, data_matrix(:,1:240), cell_ID, num_pcs);
            etoh_US_ses1 = pca2LOO(coef, data_matrix(:,241:480), cell_ID, num_pcs);
            data_l = horzcat(water_US_ses1, etoh_US_ses1);
        end
        
        % Data per trial
        counter = 0;
        data = {};
        T=60;
        
        for i = 1:2*trial_match %Change depending on trial averaged or not *3 if doing 3 sessions!
            data = [data, data_l(1:3, counter+1:counter+T,:)];
            counter = counter + T;
        end
        end
    
    % Individual mice, US    
    else
        % Per trial individual mice  
        if trial_avg == false
            data_1(:, :) = coef(:, 1:num_pcs)'*data_matrix(:,:);

            counter = 0;
            data = {};
            T=60;
            for i = 1:2*trial_match %Change depending on trial averaged or not
                data = [data, data_1(1:3, counter+1:counter+T)];
                counter = counter + T;
            end
            
             %% Loadings plot to check PCA performance
             vbls = string(1:size(data_matrix,1))';
             figure
             biplot(coef(:,1:3),'VarLabels',vbls);


            %% Neural trajectory plotting
            % Parameters
            smooth_window = [10 0];
            event_markers = [1, 30, 60];
            event_marker_size = 60;
            line_opacity = 1;
            line_width = 3;
            errorbar_opacity = 1;
            point_size = 5;
            point_freq = 5;

            scatter_plot_colors = {[0, 0, 0]/255, [0, 0, 0]/255, [0, 0, 0]/255, [0, 0, 0]/255, [255, 255, 255]/255, [255, 255, 255]/255, [255, 255, 255]/255, [255, 255, 255]/255}; % 5trajectories [255, 255, 255]/255,
            if sucrose 
                line_plot_colors = {[185 235 255]/255, [145 224 255]/255, [0 171 240]/255, [0 71 100]/255, [206 250 208]/255, [108 187 60]/255, [0 128 0]/255, [37 65 23]/255};
            else
                line_plot_colors = {[185 235 255]/255, [145 224 255]/255, [0 171 240]/255, [0 71 100]/255, [255 204 204]/255, [255 102 102]/255,  [255 0 0]/255, [204 0 0]/255};
            end    


            %% Plotting neural trajectory
            figure
            plot8traj_3Dgaussian_AH(data, smooth_window, {"Water 1", "Water 2", "Water 3", "Water 4",  ...
                "Alcohol 1", "Alcohol 2", "Alcohol 3", "Alcohol 4"}, event_markers, event_marker_size, ... 
                {"Start", "Bout", "End"}, line_plot_colors, line_opacity, line_width, scatter_plot_colors, point_size, point_freq)
            xlabel(sprintf("PC1: %0.2f%%", round(explained(1),2)))
            ylabel(sprintf("PC2: %0.2f%%", round(explained(2),2)))
            zlabel(sprintf("PC3: %0.2f%%", round(explained(3),2)))
           % legend();
            % %title(sprintf("Neural trajectory of %s - %s", mice(mouse), dates(session)));
             set(findall(gcf,'-property','FontSize'),'FontSize',14)
             set(gcf,'Position',[20 20 1000 1000])
            
        % Trial averaged individual mice    
        else
            data_1(:, :) = coef(:, 1:num_pcs)'*data_matrix(:,:);
            
            counter = 0;
            data = {};
            T=60;
            for i = 1:2
                data = [data, data_1(1:3, counter+1:counter+T)];
                counter = counter + T;
            end
            
            %% Loadings plot to check PCA performance
             vbls = string(1:size(data_matrix,1))';
             figure
             biplot(coef(:,1:3),'VarLabels',vbls);


            %% Neural trajectory plotting of trial averaged data
            % Parameters
            smooth_window = [10 0];
            event_markers = [1, 30, 60];
            event_marker_size = 60;
            line_opacity = 1;
            line_width = 3;
            errorbar_opacity = 1;
            point_size = 5;
            point_freq = 5;

            scatter_plot_colors = {[0, 0, 0]/255, [255,255,255]/255};
            line_plot_colors = {[0, 0, 255]/255, [255,0,0]/255};
            
            % Trial averaged individual mice
            figure
            plot2traj_3Dgaussian_AH(data, smooth_window, {"Water", "Alcohol"}, event_markers, event_marker_size, ... 
                {"Start", "Bout", "End"}, line_plot_colors, line_opacity, line_width, scatter_plot_colors, point_size, point_freq)
            xlabel(sprintf("PC1: %0.2f%%", round(explained(1),2)))
            ylabel(sprintf("PC2: %0.2f%%", round(explained(2),2)))
            zlabel(sprintf("PC3: %0.2f%%", round(explained(3),2)))
            %legend();
            % %title(sprintf("Neural trajectory of %s - %s", mice(mouse), dates(session)));
             set(findall(gcf,'-property','FontSize'),'FontSize',14)
             set(gcf,'Position',[20 20 1000 1000])
        end
    end 

end

