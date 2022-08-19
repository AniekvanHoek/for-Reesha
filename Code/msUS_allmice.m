%% PSTH analyses on all mice
%  Name:            Aniek van Hoek
%  Date:            02/03/2022
%  Goal:     Analyses on z-scored neural data of all mice in cued 2BC for
%               US
%
%  Input:   set parameters below
%           
%
%  Output:  Heatmaps of neuronal responses (unsorted and sorted)
%           Mean traces of all neurons
%           Responsive neurons within all neurons
%           Mean traces of coregistered neurons
%           Hierarchical clustering of coregistered neurons

close all; clearvars;

addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 
%% set parameters
pre_event = 60; %time to take before bout start in seconds
post_event = 100; %set number to frames postlick to use to determine if the cell was responsive
z_score_window = 30; %number of frames in the extracted neural actiity window to use for zscore normalization
mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];
trial_match = 10;   
sucrose = false;

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
ses = cell(5,3);

for k=1:5 % all mice
    for j = 1:3 % all sessions
      name_zscores = sprintf('zscored/%s_%s_zscored_etoh.mat', dates(j), mice(k));
      temp = load(fullfile(cd, name_zscores));
      ses{k,j}   = temp.neuron_event_zscored;
    end
end


%% Find average trace per trial type

for m = 1:size(ses,1)
    for n = 1:size(ses,2)
        neuron_event_zscored = ses{m,n};
    
        for i = 1:size(neuron_event_zscored,1) % Run through all neurons
            for k = 1:size(neuron_event_zscored, 3) % Run through all timepoints
            	mean_traces_temp(i,k) = mean(neuron_event_zscored(i,1:trial_match,k)); 
            end
        end

        mean_traces{m,n} = mean_traces_temp;
        clear mean_traces_temp;
    end
end


%% Concatenating the mean traces per neuron of all mice in a matrix
for k = 1:size(mean_traces,2)
    counter= 0;
    for i = 1:size(mean_traces,1) 
        temp = mean_traces{i,k};
        full_traces{k}(counter+1:counter+size(temp,1),:,:) = temp;
        counter = counter +size(temp,1);
    end
end

%% Heatmaps of ALL NEURON for ALL MICE

figure
Sessions = ["Pre-isolation", "Isolation", "Post-isolation"];

for i = 1:3
    subplot(1,3,i)
    imagesc(full_traces{i}(:,1:90))
    xline(pre_event,'w','Linestyle','-','LineWidth',2)
    xline(z_score_window,'w','Linestyle','-','LineWidth',2)
    colormap(redbluecmap)
    title(Sessions(i));
    caxis([-5 5])
    colorbar
end

%% First cluster and sort and then show all responses
cluster_pre = kmeans(full_traces{1}(:,[60:90]),8, 'Replicates', 20); % Set replicates to make sure it seems more deterministic (e.g. not other clusters everytime)
cluster_iso = kmeans(full_traces{2}(:,[60:90]),8,'Replicates', 20);
cluster_post = kmeans(full_traces{3}(:,[60:90]),8,'Replicates', 20);

[clusters__pre_sort, idx_pre] = sort(cluster_pre);
[clusters__iso_sort, idx_iso] = sort(cluster_iso);
[clusters__post_sort, idx_post] = sort(cluster_post);

full_traces2{1} = full_traces{1}(idx_pre,:);
full_traces2{2} = full_traces{2}(idx_iso,:);
full_traces2{3} = full_traces{3}(idx_post,:);

% Change numbers below depending on which cluster order you want 
%Pre iso
clusters_new_pre(cluster_pre==1) = 2;
clusters_new_pre(cluster_pre==2) = 6;
clusters_new_pre(cluster_pre==3) = 1;
clusters_new_pre(cluster_pre==4) = 7;
clusters_new_pre(cluster_pre==5) = 8;
clusters_new_pre(cluster_pre==6) = 5;
clusters_new_pre(cluster_pre==7) = 3;
clusters_new_pre(cluster_pre==8) = 4;

[clusters_sort_pre_new, idx_pre_new] = sort(clusters_new_pre);


% Post iso
clusters_new_iso(cluster_iso==1) = 4;
clusters_new_iso(cluster_iso==2) = 8;
clusters_new_iso(cluster_iso==3) = 6;
clusters_new_iso(cluster_iso==4) = 3;
clusters_new_iso(cluster_iso==5) = 5;
clusters_new_iso(cluster_iso==6) = 7;
clusters_new_iso(cluster_iso==7) = 1;
clusters_new_iso(cluster_iso==8) = 2;

 [clusters_sort_iso_new, idx_iso_new] = sort(clusters_new_iso);

% Rehousing
clusters_new_post(cluster_post==1) = 4;
clusters_new_post(cluster_post==2) = 2;
clusters_new_post(cluster_post==3) = 6;
clusters_new_post(cluster_post==4) = 5;
clusters_new_post(cluster_post==5) = 7;
clusters_new_post(cluster_post==6) = 8;
clusters_new_post(cluster_post==7) = 3;
clusters_new_post(cluster_post==8) = 1;

[clusters_sort_post_new, idx_post_new] = sort(clusters_new_post);

% Sorted responses
full_traces_new{1} =full_traces{1}(idx_pre_new,:);
full_traces_new{2} =full_traces{2}(idx_iso_new,:);
full_traces_new{3} =full_traces{3}(idx_post_new,:);

    figure
    Sessions = ["Pre isolation", "Isolation", "Post-isolation"];
    for session = 1:3
        subplot(1,3,session)
        imagesc(full_traces_new{session}(:,1:90))
        xline(pre_event,'w','Linestyle','-','LineWidth',2)
        xline(z_score_window,'w','Linestyle','-','LineWidth',2)
        colormap(redbluecmap)
        title(Sessions(session));
        caxis([-5 5])
        colorbar
    end
    
%% Plot mean of all cells per session

% Session 1
figure
suc = plot(mean(full_traces{1}(:,1:90))); 
err_1 = std(full_traces{1}(:,1:90))/sqrt(size(full_traces{1}(:,1:90),1));
[hl, he] = errorbar_pn(1:size(full_traces{1}(:,1:90), 2), mean(full_traces{1}(:,1:90)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);
xline(60,'k','Linestyle','--','LineWidth',1);
ylim([-0.15 0.45]);
yline(0);
xticks({});

% Session 2
figure
suc = plot(mean(full_traces{2}(:,1:90))); 
err_1 = std(full_traces{2}(:,1:90))/sqrt(size(full_traces{2}(:,1:90),1));
[hl, he] = errorbar_pn(1:size(full_traces{2}(:,1:90), 2), mean(full_traces{2}(:,1:90)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);
xline(60,'k','Linestyle','--','LineWidth',1);
xticks({});
ylim([-0.15 0.45]);
yline(0);

% Session 3
figure
suc = plot(mean(full_traces{3}(:,1:90))); 
err_1 = std(full_traces{3}(:,1:90))/sqrt(size(full_traces{3}(:,1:90),1));
[hl, he] = errorbar_pn(1:size(full_traces{3}(:,1:90), 2), mean(full_traces{3}(:,1:90)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);
xline(60,'k','Linestyle','--','LineWidth',1);
xticks({});
ylim([-0.15 0.45]);
yline(0);    

%% Get data for prism significance tests ALL
pre_mean = mean(full_traces{1}(:,1:90));
pre_err= std(full_traces{1}(:,1:90))/sqrt(size(full_traces{1}(:,1:90),1));

post_mean = mean(full_traces{2}(:,1:90));
post_err = std(full_traces{2}(:,1:90))/sqrt(size(full_traces{2}(:,1:90),1));

reh_mean= mean(full_traces{3}(:,1:90));
reh_err = std(full_traces{3}(:,1:90))/sqrt(size(full_traces{3}(:,1:90),1));


%% Responsive cells ALL neurons 
mov = 10;
pre_iso =movmean(full_traces{1}(:,31:90),mov, 2); %% ALL NEURONS
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
post_iso=movmean(full_traces{2}(:,31:90),mov, 2);
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

reh_iso=movmean(full_traces{3}(:,31:90),mov, 2);
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
resp_pre_exc = full_traces{1}(:,1:90); % all neurons
resp_pre_exc(~pre_iso_e,:) = [];

resp_pre_inh = full_traces{1}(:,1:90);
resp_pre_inh(~pre_iso_i,:) = [];

resp_post_exc = full_traces{2}(:,1:90); 
resp_post_exc(~post_iso_e,:) = [];

resp_post_inh = full_traces{2}(:,1:90);
resp_post_inh(~post_iso_i,:) = [];

resp_reh_exc = full_traces{3}(:,1:90);
resp_reh_exc(~reh_iso_e,:) = [];

resp_reh_inh = full_traces{3}(:,1:90);
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


%% Get data for prism significance tests
pre_mean_exc = mean(resp_pre_exc);
pre_err_exc = std(resp_pre_exc)/sqrt(size(resp_pre_exc,1));

post_mean_exc = mean(resp_post_exc);
post_err_exc = std(resp_post_exc)/sqrt(size(resp_post_exc,1));

reh_mean_exc = mean(resp_reh_exc);
reh_err_exc = std(resp_reh_exc)/sqrt(size(resp_reh_exc,1));

pre_mean_inh = mean(resp_pre_inh);
pre_err_inh = std(resp_pre_inh)/sqrt(size(resp_pre_inh,1));

post_mean_inh = mean(resp_post_inh);
post_err_inh = std(resp_post_inh)/sqrt(size(resp_post_inh,1));

reh_mean_inh = mean(resp_reh_inh);
reh_err_inh = std(resp_reh_inh)/sqrt(size(resp_reh_inh,1));

%% Combining neuron traces of co registered mice

for i =1:length(mice)
    for j = 1:size(mean_traces,2)
        
        temp_trace = cell2mat(mean_traces(i,j));
        registration = reg_idx{i};
        mouse_number = i.*ones(length(registration),1);
    
        for k = 1:length(registration)  
            temp_trace_co{j}(k,:) = temp_trace(registration(k,j),:);

        end
         
    end
    mouse{i} = horzcat(temp_trace_co{:}, mouse_number);
    clear temp_trace_co;
end
full_trace = vertcat(mouse{:});

%% Delete neuron 39 and 5 (ORDER IMPORTANT! only neuron 5 when 39 is already deleted!)
%full_trace(39,:)= [];
%full_trace(5,:)= [];


%% Plot mean and SEM of coregistered neurons.
% Session 1
figure
suc = plot(mean(full_trace(:,1:90))); 
err_1 = std(full_trace(:,1:90))/sqrt(size(full_trace(:,1:90),1));
[hl, he] = errorbar_pn(1:size(full_trace(:,1:90), 2), mean(full_trace(:,1:90)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);
xline(60,'k','Linestyle','--','LineWidth',1);
xticks({});
ylim([-0.45, 0.5]);

% Session 2
figure
suc = plot(mean(full_trace(:,162:251))); 
err_1 = std(full_trace(:,162:251))/sqrt(size(full_trace(:,162:251),1));
[hl, he] = errorbar_pn(1:size(full_trace(:,162:251), 2), mean(full_trace(:,162:251)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);
xline(60,'k','Linestyle','--','LineWidth',1);
xticks({});
ylim([-0.45, 0.5]);

% Session 3
figure
suc = plot(mean(full_trace(:,323:412))); 
err_1 = std(full_trace(:,323:412))/sqrt(size(full_trace(:,323:413),1));
[hl, he] = errorbar_pn(1:size(full_trace(:,323:412), 2), mean(full_trace(:,323:412)), err_1, 222, 0.4);
xline(30,'k','Linestyle','--','LineWidth',1);
xline(60,'k','Linestyle','--','LineWidth',1);
xticks({});
ylim([-0.45, 0.5]);

%% Get data for prism significance tests
pre_mean = mean(full_trace(:,1:90));
pre_err= std(full_trace(:,1:90))/sqrt(size(full_trace(:,1:90),1));

post_mean = mean(full_trace(:,162:251));
post_err = std(full_trace(:,162:251))/sqrt(size(full_trace(:,162:251),1));

reh_mean= mean(full_trace(:,323:412));
reh_err = std(full_trace(:,323:412))/sqrt(size(full_trace(:,323:412),1));


%% Hierarchical Clustering 
HC = full_trace(:,[1:90, 162:251, 323:413]); %(:,[162:251, 484:573]); %162:251, 485:574]%concatenate matrices for ethanol and sucrose events
ttl = 'Trial averaged Z-scored Calcium Trace during Sucrose and Ethanol';
cutoff = .40;   %adjust maybe .2 or .3
[Z, T, C, I , f1, outperm, cluster_order] = HAC(HC, 'ward', 'euclidean', cutoff, ttl);

%% K means clustering concatenated
clusters = kmeans(full_trace(:,[61:90, 221:251, 383:413]),4, 'Replicates', 20);
[clusters_sort, idx] = sort(clusters);
full_traces2 = full_trace(idx,[1:90, 162:251, 323:413]);
figure
imagesc(full_traces2)
colormap(redbluecmap)
caxis([-5 5])
colorbar
    
%% Find cluster numbers per mice
cellID_cluster = horzcat(full_trace(:,end), T);

cluster_per_mouse = zeros(5,6);
for mouse = 1:size(mice,2)
    for cluster = 1:size(cluster_order,1)
        for i = 1:size(cellID_cluster,1)
            if cellID_cluster(i,1) == mouse &&  cellID_cluster(i,2) == cluster
                cluster_per_mouse(mouse,cluster) = cluster_per_mouse(mouse,cluster) + 1;
            end
        end
    end
end

%% Find cluster size 
for i = 1:size(cluster_order,1)
    clusters_size(i,1) = sum(cellID_cluster(:,2)==i);
    clusters_size(i,2) = sum(cellID_cluster(:,2)==i)/119*100;
end


%% Averages plot

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
    xline(60, 'r','Linestyle','--','LineWidth',1);
   % xline(91,'k', 'Linestyle','-','LineWidth',1);
    
    
    xline(120, 'k', 'Linestyle', '--', 'LineWidth', 1);
    xline(150, 'r','Linestyle', '--', 'LineWidth', 1);
   % xline(181,'k','Linestyle','-','LineWidth',1);
    
    xline(210,'k','Linestyle','--','LineWidth',1);
    xline(240, 'r', 'Linestyle', '--', 'LineWidth', 1);
   
    xlim([0 271]);
    xticks({});
   % ylim([-2.9 5]);
% axis tight
    clear('data')
end
