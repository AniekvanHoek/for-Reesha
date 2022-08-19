%% Pseudopopulation decoding all neurons (water vs substance)
% Aniek van Hoek
%
% Goal:         Script builds pseudopopulation to decode water vs alcohol
%           
% Description:  Other scripts for sucrose exist in code/decoding
%               1) Randomly select neurons to minimum number found over all
%               sessions (425)
%               2) Get drinking data, select bouts that are far enough apart
%               3) Build pseudopopulation for train and test set separately
%               4) Decode water vs substance (including shuffled control)
%
% Input:        See parameters below
%
% Output:       AUC, Accuracy for linear SVM, LDA, random forest, 
%               logistic regression
%               List of important neurons (beta coefficients) for LDA and linear SVM

clearvars; clc
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 

%% Parameters
dates = ["01232022", "02082022", "02282022"];  % Alcohol
mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];
min_distance = 20;
n_timepoints = 20000;
n_reps = 20;
n_shuffles = 1;
n_bins = 1;
n_neurons = 425;


%% Loading neurons of all mice

neural_activity = cell(5,3);
for session =1:3
    for mouse=1:5
              name_neuraldata = sprintf('C_matrices/%s_%s_neuron.mat', dates(session), mice(mouse));
              temp = load(fullfile(cd, name_neuraldata));
              neural_activity{mouse,session}   = temp.C_matrix;
              n_neurons_mouse(mouse,session) = size(neural_activity{mouse,session},1);

              name_drinking_data = sprintf('drinking_data/%s_%s_drinking_data.mat', dates(session), mice(mouse));
              temp= load(fullfile(cd, name_drinking_data));
              drinking_data{mouse,session} = temp.drinking_data;
    end
end

%% Start repetitions
for session = 1:3
for rep = 1:n_reps
    %clearvars -except aucresult_watervsetoh result_watervsetoh_shuffle session rep min_distance n_timepoints n_reps session n_shuffles n_bins n_neurons neural_activity drinking_data n_neurons_mouse
    clear neural_activity_mouse_rep neurons_keep n_neurons_delete_mouse n_neurons_keep
    clear neural_activity_curr
    clear water_bouts_test water_bouts_train etoh_bouts_test etoh_bouts_train water_neural_test water_neural_train etoh_neural_test etoh_neural_train
    fprintf("Repetition %d\n", rep); 

    %% Take random samples of neurons for each mouse 
    n_neurons_delete = sum(n_neurons_mouse(:,session)) - n_neurons;
    for mouse =1:5
        neural_activity_mouse_rep{mouse,1} = neural_activity{mouse,session};
        n_neurons_delete_mouse(mouse) = round(n_neurons_mouse(mouse,session)/sum(n_neurons_mouse(:,session),1) * n_neurons_delete);
        n_neurons_keep(mouse) = n_neurons_mouse(mouse,session) - n_neurons_delete_mouse(mouse);
        neurons_keep{mouse} = randi(n_neurons_mouse(mouse,session), n_neurons_keep(mouse),1);

        neural_activity_mouse_rep{mouse,1} = neural_activity_mouse_rep{mouse}(neurons_keep{mouse},:);
    end
    
    n_neurons = sum(n_neurons_keep);
    %% Normalize data
    neural_activity_norm = cell(5,3);

    for mouse=1:5
            neural_activity_norm{mouse} = normalize(neural_activity_mouse_rep{mouse},2);
    end

    %% Rebin neural data and behavior
    neural_norm_rebin = cell(5,1);
    for mouse=1:5
            neural_norm_rebin_temp = NaN(size(neural_activity_norm{mouse},1), floor((size(neural_activity_norm{mouse},2)/n_bins) ));
            for i = 1:size(neural_norm_rebin_temp,2)
                neural_norm_rebin{mouse}(:,i) = mean(neural_activity_norm{mouse}(:, (i-1)*n_bins+1: (i)*n_bins),2);
            end
    end
    %% Get epochs and divide into train and test set
    water_bouts = cell(5,1);
    %neural_data_water = cell(5,1);
    etoh_bouts = cell(5,1);
    %neural_data_etoh = cell(5,1);

    for mouse = 1:5
            drinking_data_curr = drinking_data{mouse,session};

            drinking_data_rebin_temp = NaN(size(drinking_data_curr,1), floor((size(drinking_data_curr,2)/n_bins) ));
            for i = 1:size(drinking_data_rebin_temp,2)
                drinking_data_rebin_temp(:,i) = max(drinking_data_curr(:, (i-1)*n_bins+1: (i)*n_bins),[],2);
            end
            
            neural_activity_curr = neural_norm_rebin{mouse};
            [water_bouts{mouse}, ~, etoh_bouts{mouse}, ~] = get_epochs(drinking_data_rebin_temp, neural_activity_curr, min_distance/n_bins);

    end

    %% Start of making pseudopopluations 
    
    %% Split in train and test set and concatenate all the bouts
    
    for mouse = 1:5
        sample_water = randsample(size(water_bouts{mouse},2), round(size(water_bouts{mouse},2)/3)); % Sample bouts for test and train set (keeping licks in one bout together)
        water_bouts_test{mouse} =  water_bouts{mouse}(sample_water); % Concatenating all bouts together
        water_bouts_test{mouse} = horzcat(water_bouts_test{mouse}{:});
        water_bouts_train{mouse} = water_bouts{mouse}(setdiff(1:size(water_bouts{mouse},2), sample_water));
        water_bouts_train{mouse} = horzcat(water_bouts_train{mouse}{:});

        water_neural_test{mouse} = neural_norm_rebin{mouse}(:,water_bouts_test{mouse});
        water_neural_train{mouse} = neural_norm_rebin{mouse}(:,water_bouts_train{mouse});
        
        sample_etoh = randsample(size(etoh_bouts{mouse},2), round(size(etoh_bouts{mouse},2)/3));
        etoh_bouts_test{mouse} =  etoh_bouts{mouse}(sample_etoh);
        etoh_bouts_test{mouse} = horzcat(etoh_bouts_test{mouse}{:});
        etoh_bouts_train{mouse} = etoh_bouts{mouse}(setdiff(1:size(etoh_bouts{mouse},2), sample_etoh));
        etoh_bouts_train{mouse} = horzcat(etoh_bouts_train{mouse}{:});

        etoh_neural_test{mouse} = neural_norm_rebin{mouse}(:,etoh_bouts_test{mouse});
        etoh_neural_train{mouse} = neural_norm_rebin{mouse}(:,etoh_bouts_train{mouse});
    end


%% Create pseudopopulation
    % Test pseudopopulation
    behavior_test{rep,session} = NaN(n_timepoints,1);
    new_data_test{rep,session} = NaN(n_neurons, n_timepoints);


    for t = 1:n_timepoints
        beh_label = randi(2)-1;
        behavior_test{rep,session}(t) = beh_label;

        for mouse = 1:5
            potential_water_tp = find(water_bouts_test{mouse});
            potential_etoh_tp = find(etoh_bouts_test{mouse});

            if beh_label == 0
                rand_idx = randi(size(potential_water_tp,2));
                timepoint{mouse} = water_neural_test{mouse}(:,rand_idx);

            elseif beh_label == 1
                rand_idx = randi(size(potential_etoh_tp,2));
                timepoint{mouse} = etoh_neural_test{mouse}(:,rand_idx);
            end

        end

        new_data_test{rep,session}(:,t) = vertcat(timepoint{:});
    end
    
    clear timepoint rand_idx potential_etoh_tp potential_water_tp beh_label
    
    
    % Train pseudopopulation
    behavior_train{rep,session} = NaN(n_timepoints,1);
    new_data_train{rep,session} = NaN(n_neurons, n_timepoints);


    for t = 1:n_timepoints
        beh_label = randi(2)-1;
        behavior_train{rep,session}(t) = beh_label;

        for mouse = 1:5
            potential_water_tp = find(water_bouts_train{mouse});
            potential_etoh_tp = find(etoh_bouts_train{mouse});

            if beh_label == 0
                rand_idx = randi(size(potential_water_tp,2));
                timepoint{mouse} = water_neural_train{mouse}(:,rand_idx);

            elseif beh_label == 1
                rand_idx = randi(size(potential_etoh_tp,2));
                timepoint{mouse} = etoh_neural_train{mouse}(:,rand_idx);
            end

        end

        new_data_train{rep,session}(:,t) = vertcat(timepoint{:});
    end
    
    clear timepoint rand_idx potential_etoh_tp potential_water_tp beh_label
end
end
  

    %% Decoding
for session = 1:3
    parfor rep=1:n_reps
        fprintf("Decoding %d\n", rep); 
        
%         model = fitcsvm(new_data_train{rep,session}.', behavior_train{rep,session}, 'KernelFunction', 'Gaussian');
%         [preds, scores] = predict(model, new_data_test{rep,session}.');
%         [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
%         auc_watervsetoh_svm_gau(rep) = auc_local;
%         acc_watervsetoh_svm_gau(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

        model = fitcsvm(new_data_train{rep,session}.', behavior_train{rep,session}, 'KernelFunction', 'Linear');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_svm_lin(rep) = auc_local;
        acc_watervsetoh_svm_lin(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

        model = fitcensemble(new_data_train{rep,session}.', behavior_train{rep,session}, 'Method', 'Bag', 'Learners', 'tree');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_RF(rep) = auc_local;
        acc_watervsetoh_RF(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

%         model = fitclinear(new_data_train{rep,session}.', behavior_train{rep,session}, 'Learner', 'logistic');
%         [preds, scores] = predict(model, new_data_test{rep,session}.');
%         [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
%         auc_watervsetoh_lm(rep) = auc_local; 
%         acc_watervsetoh_lm(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);



        %% Shuffled control
        fprintf("Shuffle %d\n", rep); 

        % Shuffle labels
        behavior_train_par = behavior_train{rep,session}(randperm(length(behavior_train{rep,session})), 1);

        % Run models
%         model = fitcsvm(new_data_train{rep,session}.', behavior_train_par, 'KernelFunction', 'Gaussian');
%         [preds, scores] = predict(model, new_data_test{rep,session}.');
%         [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
%         auc_watervsetoh_svm_gau_shuffle(rep) = auc_local;
%         acc_watervsetoh_svm_gau_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

         model = fitcsvm(new_data_train{rep,session}.', behavior_train_par, 'KernelFunction', 'Linear');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_svm_lin_shuffle(rep) = auc_local;
        acc_watervsetoh_svm_lin_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

        model = fitcensemble(new_data_train{rep,session}.', behavior_train_par, 'Method', 'Bag', 'Learners', 'tree');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_RF_shuffle(rep) = auc_local;
        acc_watervsetoh_RF_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

%         model = fitclinear(new_data_train{rep,session}.', behavior_train_par, 'Learner', 'logistic');
%         [preds, scores] = predict(model, new_data_test{rep,session}.');
%         [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
%         auc_watervsetoh_lm_shuffle(rep) = auc_local; 
%         acc_watervsetoh_lm_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

    end
    
   % result_watervsetoh{1,session} = [auc_watervsetoh_svm_gau];
   % result_watervsetoh{2,session} = [acc_watervsetoh_svm_gau];
    result_watervsetoh{1,session} = [auc_watervsetoh_svm_lin];
    result_watervsetoh{2,session} = [acc_watervsetoh_svm_lin];
    result_watervsetoh{3,session} = [auc_watervsetoh_RF];
    result_watervsetoh{4,session} = [acc_watervsetoh_RF];
   % result_watervsetoh{7,session} = [auc_watervsetoh_lm];
   % result_watervsetoh{8,session} = [acc_watervsetoh_lm];

    save(sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/%s_result_watervsetoh',datestr(now, 'mmddyyyy')), 'result_watervsetoh')

    
    
   % result_watervsetoh_shuffle{1,session} = auc_watervsetoh_svm_gau_shuffle;
   % result_watervsetoh_shuffle{2,session} = acc_watervsetoh_svm_gau_shuffle;
    result_watervsetoh_shuffle{1,session} = auc_watervsetoh_svm_lin_shuffle;
    result_watervsetoh_shuffle{2,session} = acc_watervsetoh_svm_lin_shuffle;
    result_watervsetoh_shuffle{3,session} = auc_watervsetoh_RF_shuffle;
    result_watervsetoh_shuffle{4,session} = acc_watervsetoh_RF_shuffle;
    %result_watervsetoh_shuffle{7,session} = auc_watervsetoh_lm_shuffle;
   % result_watervsetoh_shuffle{8,session} = acc_watervsetoh_lm_shuffle;

    save(sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/%s_result_watervsetoh_shuffle',datestr(now, 'mmddyyyy')), 'result_watervsetoh_shuffle')
end    

