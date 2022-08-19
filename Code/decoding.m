%% Decoding 
% Aniek van Hoek
%
% Goal:         Decoding of ethanol vs water trials based on cue response
%
% Description:  Raw neural data is PCA transformed and then used to decode
%               if a mouse drinks alcohol or water first using a linear SVM
%               and random forest per timeframe 

% Input:    Change parameters below         
%           
% Output:   Plots decoding AUC of both models

clearvars; clc
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Parameters
T = 90;
timepoints = [1:90];
session = 1; %1 = pre isolation, 2 = post isolation, 3 = regroup
dates = ["01232022", "02082022", "02282022"];
mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];
num_shuffles = 10;
K = 5;

%% Loading neurons of all mice
cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 


% Making two cell arrays containing necessary information of all sessions
% and all trials
for k=1:5
    for j =1:3
        
        name_zscores = sprintf('nonzscored/%s_%s_cue.mat', dates(j), mice(k));
          temp = load(fullfile(cd, name_zscores));
          neural_activity{k,j}   = temp.neuron_event;

          name_trials= sprintf('trial_outcome/%s_%s_all_trials.mat',dates(j), mice(k));
          temp = load(fullfile(cd, name_trials)); 
          trial_ses{k,j}       = temp.all_trials;
 
    end
end



%% Making cell array per trial type (containing all trials) for 1 chosen session

% Finding lowest number of trials in water/ethanol choice
trial_match = 30;
for i = 1:size(mice,2)
    for j = [4,5]
        if size(trial_ses{i,session}{j},2) < trial_match
            trial_match = size(trial_ses{i,session}{j},2);
        end
    end
end

  

traces_per_trial = cell(size(mice,2),5);
for m = 1:size(neural_activity,1)
    neuron_event_zscored = neural_activity{m,session};
    
    for i = 1:size(neuron_event_zscored,1) % Run through all neurons
        for j = 1:T  % Get timepoints used in decoding
            for k = 1:5 	% Run through all possible trial types
                trials = trial_ses{m,session};
                type = trials{k};
                traces_per_trial{m,k}(i,j,:) = neuron_event_zscored(i,type,timepoints(j));
            end
        end
    end
end

%% PCA
% Include in PCA matrix which trial types need to be decoded at X for
% traces_per_trial{i,X}
pca_matrix = [];
for i = 1:size(traces_per_trial,1)
    pca_matrix = vertcat(pca_matrix, [mean(traces_per_trial{i,4}(:,:,1:trial_match),3), mean(traces_per_trial{i,5}(:,:,1:trial_match),3)]); %mean(traces_per_trial{i,3}(:,:,1:trial_match),3), 
end

[coef, ~, ~, ~, explained] = pca(pca_matrix', 'Economy', false);
num_pcs = find(cumsum(explained)>=90, 1);


%% Data for decoding
data_all = [];
labels_all = []; % 3 = choice none, 4 = choice water, 5 = choice ethanol
counter = 0;


for i = 1:size(traces_per_trial,1)
    cells_per_mice = 1:size(traces_per_trial{i,1},1);
    for j = 4:5 % for every trial type used in decoding
        for k = 1:size(traces_per_trial{i,j},3) % For every trial in that type
            data_all = cat(3, data_all, coef(counter+cells_per_mice, 1:num_pcs)'*traces_per_trial{i,j}(:,:,k));
        end
        labels_all = vertcat(labels_all, ones(size(traces_per_trial{i,j},3),1)*j);
    end
    counter = counter + size(traces_per_trial{i,j},1);
end

data_all = permute(data_all, [3, 2, 1]);
shuffle = randperm(size(data_all, 1));
data_all = data_all(shuffle, :, :);
labels_all = labels_all(shuffle, :);

%% Decoding water vs ethanol
auc_watervsetoh_glm = zeros(T, K);
auc_watervsetoh_RF= zeros(T, K);

data_water = data_all(labels_all==4, :, :);
data_etoh = data_all(labels_all==5, :, :);
N_trials_water = sum(labels_all==4);
N_trials_etoh = sum(labels_all==5);

for k = 1: K
    parfor t = 1: T
        data_test_water= squeeze(data_water(floor((k-1)*N_trials_water/K) + 1: floor(k*N_trials_water/K), t, :));      
        data_test_etoh= squeeze(data_etoh(floor((k-1)*N_trials_etoh/K) + 1: floor(k*N_trials_etoh/K), t, :));
        data_train_water = squeeze(data_water(setdiff(1: N_trials_water, floor((k-1)*N_trials_water/K) + 1: floor(k*N_trials_water/K)), t, :));
        data_train_etoh = squeeze(data_etoh(setdiff(1: N_trials_etoh, floor((k-1)*N_trials_etoh/K) + 1: floor(k*N_trials_etoh/K)), t, :));
        data_train = cat(1, data_train_water, data_train_etoh);
        labels_train = cat(1, ones(size(data_train_water, 1), 1), zeros(size(data_train_etoh, 1), 1)); % Here labels are set to 1 and 0 again
        data_test = cat(1, data_test_water, data_test_etoh);
        labels_test = cat(1, ones(size(data_test_water, 1), 1), zeros(size(data_test_etoh, 1), 1));

        model = fitcensemble(data_train, labels_train, 'Method', 'Bag', 'Learners', 'tree');
        [~, scores] = predict(model, data_test);
        [~, ~, ~, auc_local] = perfcurve(labels_test, scores(:, 2)', 1);
        auc_watervsetoh_RF(t, k) = auc_local;
        
        model = mnrfit(data_train, labels_train+1); % specialized function that only does logistic regression. GLM is linear combination, the linear combi is pushed into link function. If its logistic or sigmoids its logistic regression. If poission its a poisson process.
        score = mnrval(model, data_test);
        [~, ~, ~, auc_local] = perfcurve(labels_test, score(:, 2)', 1);
        auc_watervsetoh_glm(t, k) = auc_local;
    end
end

clearvars data_water data_etoh N_trials_water N_trials_etoh k t data_test_water
clearvars data_test_etoh data_train_water data_train_etoh data_train
clearvars labels_train data_test labels_test model score auc_local

%% Decoding shuffled control
auc_watervsetoh_glm_control = cell(1, num_shuffles);
auc_watervsetoh_RF_control = cell(1, num_shuffles);

data_water = data_all(labels_all==4, :, :);
data_etoh = data_all(labels_all==5, :, :);
N_trials_water = sum(labels_all==4);
N_trials_etoh = sum(labels_all==5);

for s = 1: num_shuffles
    auc_watervsetoh_glm_control_local = zeros(T, K);
    auc_watervsetoh_RF_control_local = zeros(T, K);
    for k = 1: K
        parfor t = 1: T
            data_test_water = squeeze(data_water(floor((k-1)*N_trials_water/K) + 1: floor(k*N_trials_water/K), t, :));      
            data_test_etoh= squeeze(data_etoh(floor((k-1)*N_trials_etoh/K) + 1: floor(k*N_trials_etoh/K), t, :));
            data_train_water = squeeze(data_water(setdiff(1: N_trials_water, floor((k-1)*N_trials_water/K) + 1: floor(k*N_trials_water/K)), t, :));
            data_train_etoh = squeeze(data_etoh(setdiff(1: N_trials_etoh, floor((k-1)*N_trials_etoh/K) + 1: floor(k*N_trials_etoh/K)), t, :));
            data_train = cat(1, data_train_water, data_train_etoh);
            labels_train = cat(1, ones(size(data_train_water, 1), 1), zeros(size(data_train_etoh, 1), 1));
            data_test = cat(1, data_test_water, data_test_etoh);
            labels_test = cat(1, ones(size(data_test_water, 1), 1), zeros(size(data_test_etoh, 1), 1));
            
            labels_train = labels_train(randperm(length(labels_train)), 1);

            model = fitcensemble(data_train, labels_train, 'Method', 'Bag', 'Learners', 'tree');
            [~, score] = predict(model, data_test);
            [~, ~, ~, auc_local] = perfcurve(labels_test, score(:, 2)', 1);
            auc_watervsetoh_RF_control_local(t, k) = auc_local;

            model = mnrfit(data_train, labels_train+1);
            score = mnrval(model, data_test);
            [~, ~, ~, auc_local] = perfcurve(labels_test, score(:, 2)', 1);
            auc_watervsetoh_glm_control_local(t, k) = auc_local;
        end
    end
    auc_watervsetoh_RF_control{1, s} = auc_watervsetoh_RF_control_local;
    auc_watervsetoh_glm_control{1, s} = auc_watervsetoh_glm_control_local;
end

clearvars num_shuffles K data_water data_etoh N_trials_water N_trials_etoh s
clearvars auc_watervsetoh_glm_control_local
clearvars auc_watervsetoh_svmgauss_control_local k t data_test_water
clearvars data_test_etoh data_train_water data_train_etoh data_train
clearvars labels_train data_test labels_test model score auc_local

%% Plot
plot_model_performance_Kanha(T, {auc_watervsetoh_RF, cat(2, auc_watervsetoh_RF_control{:})}, ...
   {"blue", "black"}, {"Test Data AUC", "Test Data AUC: Shuffled"}, ...
   "Time", "Score", "RF", "on");

plot_model_performance_Kanha(T, {auc_watervsetoh_glm, cat(2, auc_watervsetoh_glm_control{:})}, ...
   {"blue", "black"}, {"Test Data AUC", "Test Data AUC: Shuffled"}, ...
   "Time", "Score", "GLM", "on");

%% Statistics
% Average
[h_mean,p_mean,ks2stat_mean] = kstest2(mean(auc_watervsetoh_RF), mean(cat(2, auc_watervsetoh_RF_control{:})));

h = [];
p = [];
ks2stat = [];
control = cat(2, auc_watervsetoh_RF_control{:});

% Per time point 
for i = 1:size(auc_watervsetoh_glm,1)
    [h(i),p(i), ks2stat(i)] = kstest2(auc_watervsetoh_RF(i,:), control(i,:));
end
