%% PSTH analyses on all mice
%  Name:            Aniek van Hoek
%  Date:            02/03/2022
%  Goal:     Analyses on z-scored neural data of all mice in cue
%                   task experiments
%             
%  Trial options:
%                   - A: only water
%                   - B: only ethanol
%                   - C: both
%
%  Trial order:     C, A, C, C, B, C, C, B, C, A, C, C, B, C, C, C, A, C, A, C, C, C, B, C, C, A, C, B, C, C
%     
%  Input:   set parameters below
%           
%
%  Output:  Heatmaps of neuronal responses
%           Responsive neurons
%           Mean traces of coregistered neurons
%           Hierarchical clustering of coregistered neurons
%           Correlation between AUC of traces and rank

clearvars; close all;
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 

%% set parameters
sucrose = false;
pre_event = 30; %time to take before bout start in seconds
post_event = 1130; %set number to frames postlick to use to determine if the cell was responsive
z_score_window = 30; %number of frames in the extracted neural actiity window to use for zscore normalization
mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];
trial_match = 5;

%% Loading co registration of all mice
if sucrose
    dates = ["01082022", "02112022", "03012022"];
    reg_idx = cell(5,1);
    for mouse = 1:5
        filename = sprintf("celleg/Suc/%s_task_suc_cellreg.mat", mice(mouse));
        temp = load(fullfile(cd,filename));
        reg_idx{mouse} =  temp.cell_registered_struct.cell_to_index_map;
        reg_idx{mouse}(any(reg_idx{mouse}==0,2),:) = [];
        n_neurons_mouse(mouse) = size(reg_idx{mouse}(:,1),1);
    end
else
    dates = ["01232022", "02082022", "02282022"];
    reg_idx = cell(5,1);
    for mouse = 1:5
        filename = sprintf("celleg/Etoh/%s_task_etoh_cellreg.mat", mice(mouse));
        temp = load(fullfile(cd,filename));
        reg_idx{mouse} =  temp.cell_registered_struct.cell_to_index_map;
        reg_idx{mouse}(any(reg_idx{mouse}==0,2),:) = [];
        n_neurons_mouse(mouse) = size(reg_idx{mouse}(:,1),1);
    end
end

%% Loading neurons of all mice
neural_ses = cell(5,2);
trial_ses = cell(5,2);

for k=1:5 % all mice
    for j =1:3 % all sessions
        
        name_zscores = sprintf('zscored/%s_%s_zscored_cue.mat', dates(j), mice(k));
          temp = load(fullfile(cd, name_zscores));
          neural_ses{k,j}   = temp.neuron_event_zscored;

          name_trials= sprintf('trial_outcome/%s_%s_all_trials.mat',dates(j), mice(k));
          temp = load(fullfile(cd, name_trials)); 
          trial_ses{k,j}       = temp.all_trials;
    end
end


%% Trial average
for m = 1:size(neural_ses,1) % for all mice
    for n = 1:size(neural_ses,2) % for all sessions

        neuron_event_zscored = neural_ses{m,n};
    
    % If trying a trial match lower than 5, change for loop k = [1...]
        for i = 1:size(neuron_event_zscored,1) % Run through all neurons
            for j = 1:size(neuron_event_zscored,3)
                for k = [5] 
                    trials = trial_ses{m,n};
                    type = trials{k};
                    mean_traces_temp(i,j,k) = mean(neuron_event_zscored(i,type(1:trial_match),j));
                end

                for k = [1,2,3,4]  
                    trials = trial_ses{m,n};
                    type = trials{k};
                    mean_traces_temp(i,j,k) = mean(neuron_event_zscored(i,type,j));
                end
            end
        end

        mean_traces{m,n} = mean_traces_temp;
        clear mean_traces_temp;
    end
    
end

%% Concatenating the mean traces per neuron of all mice in a matrix

for k=1:size(mean_traces,2)
    counter= 0;
    for i = 1:size(mean_traces,1)
        temp = mean_traces{i,k};
        full_traces{k}(counter+1:counter+size(temp,1),:,:) = temp;
        counter = counter +size(temp,1);
    end
end

%% Heatmaps of ALL neurons of ALL mice for choice water and choice ethanol

Sessions = ["Pre isolation", "Post Isolation", "Regroup"];
figure
for i=1:size(Sessions,2)
subplot(1,3,i)
imagesc(full_traces{i}(:,1:60,5))
xline(z_score_window,'w','Linestyle','-','LineWidth',2)
colormap(redbluecmap)
title(Sessions(i));
caxis([-5 5])
colorbar
end

%% Responsive cells

mov = 10;
pre_iso =movmean(full_traces{1}(:,31:60,5),mov, 2); %% ALL NEURONS
[highestvalue, loc_highest] = (max(abs(pre_iso), [],2));

for i = 1:length(loc_highest)
    pre_iso_final(i,1) = pre_iso(i,loc_highest(i));
end

threshold = 1.96;    

pre_iso_all = abs(pre_iso_final)>threshold;
pre_iso_e = pre_iso_final>=threshold;
pre_iso_i = pre_iso_final<=threshold/-1;

results = [];
results.pre_iso_all_percent = sum(pre_iso_all)/size(pre_iso_all,1);
results.pre_iso_e_percent = sum(pre_iso_e)/size(pre_iso_e,1);
results.pre_iso_i_percent = sum(pre_iso_i)/size(pre_iso_i,1);

% Post
post_iso=movmean(full_traces{2}(:,31:60,5),mov, 2);
%post_iso = movmean(full_ethanol(:,31:60), mov,2);
[highestvalue_post, loc_highest_post] = (max(abs(post_iso), [],2));


for i = 1:length(loc_highest_post)
    post_iso_final(i,1) = post_iso(i,loc_highest_post(i));
end

post_iso_all = abs(post_iso_final)>threshold;
post_iso_e = post_iso_final>=threshold;
post_iso_i = post_iso_final<=threshold/-1;

results.post_iso_all_percent = sum(post_iso_all)/size(post_iso_all,1);
results.post_iso_e_percent = sum(post_iso_e)/size(post_iso_e,1);
results.post_iso_i_percent = sum(post_iso_i)/size(post_iso_i,1);

% Rehousing

reh_iso=movmean(full_traces{3}(:,31:60,5),mov, 2);
[highestvalue_reh, loc_highest_reh] = (max(abs(reh_iso), [],2));


for i = 1:length(loc_highest_reh)
    reh_iso_final(i,1) = reh_iso(i,loc_highest_reh(i));
end

reh_iso_all = abs(reh_iso_final)>threshold;
reh_iso_e = reh_iso_final>=threshold;
reh_iso_i = reh_iso_final<=threshold/-1;

results.reh_iso_all_percent = sum(reh_iso_all)/size(reh_iso_all,1);
results.reh_iso_e_percent = sum(reh_iso_e)/size(reh_iso_e,1);
results.reh_iso_i_percent = sum(reh_iso_i)/size(reh_iso_i,1);

%% Plots of mean and STD of excited and inhibited cells (responsive)
resp_pre_exc = full_traces{1}(:,1:60,5); % all neurons
resp_pre_exc(~pre_iso_e,:) = [];

resp_pre_inh = full_traces{1}(:,1:60,5);
resp_pre_inh(~pre_iso_i,:) = [];

resp_post_exc = full_traces{2}(:,1:60,5); 
resp_post_exc(~post_iso_e,:) = [];

resp_post_inh = full_traces{2}(:,1:60,5);
resp_post_inh(~post_iso_i,:) = [];

resp_reh_exc = full_traces{3}(:,1:60,5);
resp_reh_exc(~reh_iso_e,:) = [];

resp_reh_inh = full_traces{3}(:,1:60,5);
resp_reh_inh(~reh_iso_i,:) = [];


figure
pre_exc = plot(mean(resp_pre_exc), 'k');
err_1 = std(resp_pre_exc)/sqrt(size(resp_pre_exc,1));
[hl, he] = errorbar_pn(1:size(resp_pre_exc, 2), mean(resp_pre_exc), err_1, 111, 0.4);
hold on

post_exc = plot(mean(resp_post_exc), 'k');
err_1 = std(resp_post_exc)/sqrt(size(resp_post_exc,1));
[hl, he] = errorbar_pn(1:size(resp_post_exc, 2), mean(resp_post_exc), err_1, 222, 0.4);
hold on

reh_exc = plot(mean(resp_reh_exc), 'k');
err_1 = std(resp_reh_exc)/sqrt(size(resp_reh_exc,1));
[hl, he] = errorbar_pn(1:size(resp_reh_exc, 2), mean(resp_reh_exc), err_1, 333, 0.4);
%ylim([-0.25 0.35])
hold off


figure
pre_inh= plot(mean(resp_pre_inh), 'k');
err_1 = std(resp_pre_inh)/sqrt(size(resp_pre_inh,1));
[hl, he] = errorbar_pn(1:size(resp_pre_inh, 2), mean(resp_pre_inh), err_1, 111, 0.4);
hold on

post_inh = plot(mean(resp_post_inh), 'k');
err_1 = std(resp_post_inh)/sqrt(size(resp_post_inh,1));
[hl, he] = errorbar_pn(1:size(resp_post_inh, 2), mean(resp_post_inh), err_1, 222, 0.4);
hold on

reh_inh = plot(mean(resp_reh_inh), 'k');
err_1 = std(resp_reh_inh)/sqrt(size(resp_reh_inh,1));
[hl, he] = errorbar_pn(1:size(resp_reh_inh, 2), mean(resp_reh_inh), err_1, 333, 0.4);
%ylim([-0.25 0.35])
hold off



%% Getting neuron traces of co registered mice for HC
mouse = cell(5,1);
for i = 1:length(mice) % all mice
    for j = 1:size(mean_traces,2) % all sessions
        
        temp_trace = cell2mat(mean_traces(i,j));
        registration = reg_idx{i};
        mouse_number = i.*ones(length(registration),1);
        
        for k = 1:size(registration,1)
            for l = 1:size(temp_trace,3)
                temp_trace_trial{l}(k,:) = temp_trace(registration(k,j),1:100,l);
            end
        end
        temp_trace_co{j} = horzcat(temp_trace_trial{:});
        clear temp_trace_trial;
    end
    mouse{i} = horzcat(temp_trace_co{:}, mouse_number);
    clear temp_trace_co;
end
full_trace = vertcat(mouse{:});

         
%% Plot mean of coregistered cells per session

% Session 1
figure
suc = plot(mean(full_trace(:,401:480))); 
err_1 = std(full_trace(:,401:480))/sqrt(size(full_trace(:,401:480),1));
[hl, he] = errorbar_pn(1:size(full_trace(:,401:480), 2), mean(full_trace(:,401:480)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);

% Session 2
figure
suc = plot(mean(full_trace(:,901:980))); 
err_1 = std(full_trace(:,901:980))/sqrt(size(full_trace(:,901:980),1));
[hl, he] = errorbar_pn(1:size(full_trace(:,901:980), 2), mean(full_trace(:,901:980)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);

% Session 3
figure
suc = plot(mean(full_trace(:,1401:1500))); 
err_1 = std(full_trace(:,1401:1480))/sqrt(size(full_trace(:,1401:1480),1));
[hl, he] = errorbar_pn(1:size(full_trace(:,1401:1480), 2), mean(full_trace(:,1401:1480)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);


%% Hierarchical Clustering 
HC = full_trace(:,[1:60, 501:560, 1001:1060]); %(:,[90:180, 585:675, 1080:1170]); %[435:495, 930:990, 1425:1485]) ;%(:,[401:460]); %, 901:960, 1401:1460]); %concatenate matrices for ethanol and sucrose events
ttl = 'Trial averaged Z-scored Calcium Trace during Ethanol pre and post iso';
cutoff = .50;   %adjust maybe .2 or .3
[Z, T, C, I , f1, outperm, cluster_order] = HAC(HC, 'ward', 'euclidean', cutoff, ttl);


%% Mean traces per cluster
% If error, look into clusters with only 1 cell


label = [];
data = [];

figure
for i = 1:length(cluster_order)
    
    idx = find(T==i);
       
    for j = 1:length(idx);       
        data(j,:) = HC(idx(j),:);   
    end
    
    err = std(data)/sqrt(size(data,1));
   
    subplot(length(cluster_order),1, i);
    
    [hl, he] = errorbar_pn(1:size(HC,2), mean(data), err, 123, 0.4);

    hold on
    cluster = num2str(i);
    title(cluster);
      xline(30,'k','Linestyle','--','LineWidth',1);
      xline(90,'k','Linestyle','--','LineWidth',1);
      xline(150, 'k', 'Linestyle', '--', 'LineWidth', 1);
%     
     xline(60,'k','Linestyle','-','LineWidth',1);
     xline(120,'k','Linestyle','-','LineWidth',1);
%     
    ylim([-1 6.5]);
    %axis tight
   % xlim([0 270]);
    clear('data')
end

%% Find which mouse is which cluster
cellID_cluster = horzcat(full_trace(:,1501), T);




%% Correlation to social rank of AUC of mean of coregistered cells
for m = 1:size(mice,2)

    for n=1:size(mean_traces,2)

        temp = mean_traces{m,n};
        registration = reg_idx{m};
    
        for j = 1:size(temp,3)

                for l = 1:size(registration,1)
                    AUC{n}(m,l,j) = trapz(temp(registration(l,n),1:664,j));
                end

        end
    end
end

for m = 1:size(mice,2)
    for n = 1:size(mean_traces,2)
        for j = 1:size(temp,3)
            gem_AUC{n}(m,j) = mean(AUC{n}(m,:,j),2);
        end
    end
end

AUC_ses1(:,6) = [3,2,1,1,2,3];
AUC_ses2(:,6) = [3,2,1,1,2,3];

figure
scatter(AUC_ses1(:,6), AUC_ses1(:,5), 50, 'filled')
xlim([0.9 3]);
xticks([1 2 3]);
xlabel("Rank")
ylabel("AUC sucrose")

figure
scatter(AUC_ses2(:,6), AUC_ses2(:,5), 50, 'filled')
xlim([0.9 3]);
xticks([1 2 3]);
xlabel("Rank")
ylabel("AUC ethanol")
