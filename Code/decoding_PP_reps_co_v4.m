%% Pseudopopulation decoding co-registered neurons (water vs substance)
% Aniek van Hoek
%
% Goal:         Script builds pseudopopulation to decode water vs alcohol
%           
% Description:  Other scripts for sucrose exist in code/decoding
%               1) Get drinking data, select bouts that are far enough apart
%               2) Build pseudopopulation for train and test set separately
%               3) Decode water vs substance (including shuffled control)
%
% Input:        See parameters below
%
% Output:       AUC, Accuracy for linear SVM, LDA random forest, 
%               logistic regression
%               List of important neurons (beta coefficients) for LDA and linear SVM

clearvars; clc
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 

%% Parameters
dates = ["01232022", "02082022", "02282022"]; % Alcohol, different script for sucrose
mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];
min_distance = 20; % How far bouts have to be apart (20 frames = 2 seocnds)
n_timepoints = 20000; % Size of new data matrix
n_reps = 20; % Number of repetitions, comparable to 20 cross validations
n_shuffles = 1;
n_bins = 1; % Data can be rebinned, more course data (higher n_bins) makes decoding easier

%% Loading co registration of all mice
reg_idx = cell(5,1);
for i = 1:5
    filename = sprintf("celleg/Etoh/%s_task_etoh_cellreg.mat", mice(i));
    temp = load(fullfile(cd,filename));
    reg_idx{i} =  temp.cell_registered_struct.cell_to_index_map;
    reg_idx{i}(any(reg_idx{i}==0,2),:) = [];
    n_neurons_mouse(i) = size(reg_idx{i}(:,1),1);
end

n_neurons = sum(n_neurons_mouse);


%% Loading neurons of all mice

neural_activity = cell(5,3);
for session =1:3
    for mouse=1:5
        registration = reg_idx{mouse};
              name_neuraldata = sprintf('C_matrices/%s_%s_neuron.mat', dates(session), mice(mouse));
              temp = load(fullfile(cd, name_neuraldata));
              temp_neural_activity{mouse,session}   = temp.C_matrix;
              
              % Keep only cells coregistered over all 3  sessions
              for i =1:size(registration,1) 
                    neural_activity{mouse,session}(i,:) = temp_neural_activity{mouse,session}(registration(i,session),:);
              end           
              neural_activity{mouse,session}(all(~neural_activity{mouse,session},2), : ) = [];
              
              % Drinking data
              name_drinking_data = sprintf('drinking_data/%s_%s_drinking_data.mat', dates(session), mice(mouse));
              temp= load(fullfile(cd, name_drinking_data));
              drinking_data{mouse,session} = temp.drinking_data;
    end
end


%% Start repetitions
for session = 1:3
for rep = 1:n_reps
    clear neural_activity_curr
    clear water_bouts_test water_bouts_train etoh_bouts_test etoh_bouts_train water_neural_test water_neural_train etoh_neural_test etoh_neural_train
    fprintf("Repetition %d\n", rep); 

    %% Normalize data
    neural_activity_norm = cell(5,1);

    for mouse=1:5
        neural_activity_norm{mouse} = normalize(neural_activity{mouse,session},2);
    end

    %% Rebin neural data and behavior
    neural_norm_rebin = cell(5,1);
    for mouse=1:5
            neural_norm_rebin_temp = NaN(size(neural_activity_norm{mouse},1), floor((size(neural_activity_norm{mouse},2)/n_bins) ));
            for i = 1:size(neural_norm_rebin_temp,2)
                neural_norm_rebin{mouse}(:,i) = mean(neural_activity_norm{mouse}(:, (i-1)*n_bins+1: (i)*n_bins),2);
            end
    end
    %% Get epochs (bouts) and divide into train and test set
    water_bouts = cell(5,1);
    etoh_bouts = cell(5,1);

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
    
    %% Split in train and test set and concatenate all the bouts into one array
    
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


    for t = 1:n_timepoints                              % For every timepoints
        beh_label = randi(2)-1;                         % Choose random behavior label
        behavior_test{rep,session}(t) = beh_label;

        for mouse = 1:5
            potential_water_tp = find(water_bouts_test{mouse}); % Get array of bouts for every mouse
            potential_etoh_tp = find(etoh_bouts_test{mouse});

            if beh_label == 0                                   % Choose random neural timepoint for behavior label for each mouse
                rand_idx = randi(size(potential_water_tp,2));
                timepoint{mouse} = water_neural_test{mouse}(:,rand_idx);

            elseif beh_label == 1
                rand_idx = randi(size(potential_etoh_tp,2));
                timepoint{mouse} = etoh_neural_test{mouse}(:,rand_idx);
            end

        end

        new_data_test{rep,session}(:,t) = vertcat(timepoint{:}); % Vertical concatenation of all neurons (of all mice) at that timepoints
    end
    
    clear timepoint rand_idx potential_etoh_tp potential_water_tp beh_label
    
    
    % Train pseudopopulation (same as above)
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
        
        % Gaussian SVM
%         model = fitcsvm(new_data_train{rep,session}.', behavior_train{rep,session}, 'KernelFunction', 'Gaussian');
%         [preds, scores] = predict(model, new_data_test{rep,session}.');
%         [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
%         auc_watervsetoh_svm_gau(rep) = auc_local;
%         acc_watervsetoh_svm_gau(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

        % Linear SVM
        model = fitcsvm(new_data_train{rep,session}.', behavior_train{rep,session}, 'KernelFunction', 'Linear');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_svm_lin(rep) = auc_local;
        acc_watervsetoh_svm_lin(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);
        parameters_SVM(:,rep) = model.Beta; % Save weights of variables (neurons)
        
        % Linear Discriminant Classifier (LDA)
        Mdl = fitcdiscr(new_data_train{rep,session}.', behavior_train{rep,session});
        [preds, scores] = predict(Mdl, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_LDA(rep) = auc_local;
        acc_watervsetoh_LDA(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);
        % Get out weights of variables (neurons)        
        [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order! 
        lambda = diag(LAMBDA);
        [lambda, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        parameters_LDA(:,rep) = W(:,1);    

        % Random Forest    
        model = fitcensemble(new_data_train{rep,session}.', behavior_train{rep,session}, 'Method', 'Bag', 'Learners', 'tree');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_RF(rep) = auc_local;
        acc_watervsetoh_RF(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

        % Logistic regression    
        model = fitclinear(new_data_train{rep,session}.', behavior_train{rep,session}, 'Learner', 'logistic');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_lm(rep) = auc_local; 
        acc_watervsetoh_lm(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);



        %% Shuffled control
        fprintf("Shuffle %d\n", rep); 

        % Shuffle labels
        behavior_train_par = behavior_train{rep,session}(randperm(length(behavior_train{rep,session})), 1);

        % Run models
        % Gaussian SVM
%         model = fitcsvm(new_data_train{rep,session}.', behavior_train_par, 'KernelFunction', 'Gaussian');
%         [preds, scores] = predict(model, new_data_test{rep,session}.');
%         [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
%         auc_watervsetoh_svm_gau_shuffle(rep) = auc_local;
%         acc_watervsetoh_svm_gau_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);
        
        % Linear SVM
        model = fitcsvm(new_data_train{rep,session}.', behavior_train_par, 'KernelFunction', 'Linear');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_svm_lin_shuffle(rep) = auc_local;
        acc_watervsetoh_svm_lin_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);
        parameters_shuffle_SVM(:,rep) = model.Beta;
        
        % Linear discriminant analysis (LDA)
        Mdl = fitcdiscr(new_data_train{rep,session}.', behavior_train_par);
        [preds, scores] = predict(Mdl, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_LDA_shuffle(rep) = auc_local;
        acc_watervsetoh_LDA_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);
        [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); %Must be in the right order! 
        lambda = diag(LAMBDA);
        [lambda, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        parameters_shuffle_LDA(:,rep) = W(:,1);

        % Random Forest
        model = fitcensemble(new_data_train{rep,session}.', behavior_train_par, 'Method', 'Bag', 'Learners', 'tree');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_RF_shuffle(rep) = auc_local;
        acc_watervsetoh_RF_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);
        
        % Logistic regression
        model = fitclinear(new_data_train{rep,session}.', behavior_train_par, 'Learner', 'logistic');
        [preds, scores] = predict(model, new_data_test{rep,session}.');
        [~, ~, ~, auc_local] = perfcurve(behavior_test{rep,session}, scores(:, 2)', 1);
        auc_watervsetoh_lm_shuffle(rep) = auc_local; 
        acc_watervsetoh_lm_shuffle(rep) = sum((preds == behavior_test{rep,session}))/size(behavior_test{rep,session},1);

    end
    
    result_watervsetoh_co{1,session} = [auc_watervsetoh_svm_lin];
    result_watervsetoh_co{2,session} = [acc_watervsetoh_svm_lin];
    result_watervsetoh_co{3,session} = parameters_SVM;
    result_watervsetoh_co{4,session} = [auc_watervsetoh_LDA];
    result_watervsetoh_co{5,session} = [acc_watervsetoh_LDA];
    result_watervsetoh_co{6,session} = parameters_LDA;
    result_watervsetoh_co{7,session} = [auc_watervsetoh_RF];
    result_watervsetoh_co{8,session} = [acc_watervsetoh_RF];
    result_watervsetoh_co{9,session} = [auc_watervsetoh_lm];
    result_watervsetoh_co{10,session} = [acc_watervsetoh_lm];

    save(sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/%s_result_watervsetoh_co',datestr(now, 'mmddyyyy')), 'result_watervsetoh_co')

   
    result_watervsetoh_shuffle_co{1,session} = [auc_watervsetoh_svm_lin_shuffle];
    result_watervsetoh_shuffle_co{2,session} = [acc_watervsetoh_svm_lin_shuffle];
    result_watervsetoh_shuffle_co{3,session} = parameters_shuffle_SVM;
    result_watervsetoh_shuffle_co{4,session} = [auc_watervsetoh_LDA_shuffle];
    result_watervsetoh_shuffle_co{5,session} = [acc_watervsetoh_LDA_shuffle];
    result_watervsetoh_shuffle_co{6,session} = parameters_shuffle_LDA;
    result_watervsetoh_shuffle_co{7,session} = [auc_watervsetoh_RF_shuffle];
    result_watervsetoh_shuffle_co{8,session} = [acc_watervsetoh_RF_shuffle];
    result_watervsetoh_shuffle_co{9,session} = [auc_watervsetoh_lm_shuffle];
    result_watervsetoh_shuffle_co{10,session} = [acc_watervsetoh_lm_shuffle];
       
    save(sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/%s_result_watervsetoh_shuffle_co',datestr(now, 'mmddyyyy')), 'result_watervsetoh_shuffle_co')
end    

